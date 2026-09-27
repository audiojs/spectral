// Phase-vocoder engine — the spectral bin-scatter/phase-lock machinery behind
// pitch shifting (@audio/shift-pvoc*), time stretching (@audio/stretch-pvoc*), and
// spectral freezing. Operates on STFT magnitude/phase frames; pair with
// fourier-transform/stft (ana/syn hop) or @audio/stft (AMS).
//
//   - scatterGated: Bernsee/SMB peak-gated bin scatter (pitch shift by ratio)
//   - scatterLocked: Laroche-Dolson rigid-ROI peak-locked scatter (phase coherent)
//   - lockPhase: lock non-peak bin phases to their nearest peak's rotation (time stretch)
//   - lockMap: which peak each bin locks to, the same regions, for complex-domain locking
//   - lockState / lockAdvance: the phase-locked vocoder step on complex bins (no polar round trip)
//   - findPeaks / nearestPeak / makeFrameRatio / wrapPhase: the supporting primitives

export const PI2 = Math.PI * 2

export function wrapPhase(p) {
  return p - Math.floor(p / PI2 + 0.5) * PI2  // floor(x + 0.5) rounds like Math.round here, several times faster in V8 (hot: per bin per frame)
}

// First-order local magnitude peaks above a fraction of the frame's peak.
// ±1 comparison keeps closely-spaced chord partials whose mainlobes overlap. `>=` on the
// left / `>` on the right reports the trailing edge of an exact-magnitude plateau exactly
// once, instead of a strict `>` on both sides missing the whole plateau.
export function findPeaks(mag, half) {
  let maxM = 0
  for (let k = 0; k <= half; k++) if (mag[k] > maxM) maxM = mag[k]
  let floor = Math.max(1e-8, maxM * 0.005)
  let peaks = []
  for (let k = 1; k < half; k++) {
    let v = mag[k]
    if (v < floor) continue
    if (v >= mag[k - 1] && v > mag[k + 1]) peaks.push(k)
  }
  return peaks
}

// Binary-search nearest peak index for bin k.
export function nearestPeak(peaks, k) {
  if (!peaks.length) return -1
  let lo = 0, hi = peaks.length - 1
  while (lo < hi) {
    let mid = (lo + hi) >> 1
    if (peaks[mid] < k) lo = mid + 1
    else hi = mid
  }
  if (lo > 0 && Math.abs(peaks[lo - 1] - k) <= Math.abs(peaks[lo] - k)) return lo - 1
  return lo
}

// Peak-gated bin scatter (Bernsee/SMB scheme): every analysis bin at or adjacent (±1) to a
// local magnitude peak advances phase at its own instantaneous frequency, scales it by
// `ratio`, and deposits into the destination bin that frequency implies. Colliding bins
// accumulate in the energy domain (Σmag², √ at the end) — synthesis treats each bin as an
// independent oscillator, so energies add where magnitude-summing overshoots (+4.3 dB for
// a Hann mainlobe's own ±1 bins landing together) and last-writer-wins discards every
// other contributor. The frequency written to a bin is its loudest contributor's (a
// quieter contributor's frequency estimate is masked anyway).
// The gate keeps only mainlobe cores; skirt bins it drops carry real energy belonging to
// the same partials. The frame is renormalized so kept-bin energy matches the input
// frame's — minus content whose destination fell outside Nyquist, which is legitimately
// lost — times WIN_GAIN: concentrating a windowed mainlobe into one bin makes the ISTFT
// frame a pure sinusoid where the analysis frame was a windowed one, and through the
// engine's w·(·)/Σw² overlap-add that costs exactly mean(w)/rms(w). Per-frame and causal,
// so batch and stream reconstruct at identical loudness with no whole-signal correction.
// `newMag`/`newFreq`/`peakMag` are caller-owned scratch sized `half+1`, zero-filled by the
// caller before the call. `prevPhase` is the previous frame's unwrapped phase, or `null`
// on the first frame.
export function scatterGated(mag, phase, prevPhase, ratio, ctx, newMag, newFreq, peakMag) {
  let { half, hop, freqPerBin } = ctx
  let maxM = 0
  for (let k = 0; k <= half; k++) if (mag[k] > maxM) maxM = mag[k]
  let floor = Math.max(1e-8, maxM * 0.005)
  let eIn = 0, eOut = 0
  for (let k = 0; k <= half; k++) {
    let e = mag[k] * mag[k]
    eIn += e
    let eligible = false
    for (let d = -1; d <= 1; d++) {
      let j = k + d
      if (j <= 0 || j >= half) continue
      if (mag[j] >= floor && mag[j] >= mag[j - 1] && mag[j] > mag[j + 1]) { eligible = true; break }
    }
    if (!eligible) continue
    let trueFreq
    if (!prevPhase) trueFreq = k * freqPerBin
    else {
      let dp = wrapPhase(phase[k] - prevPhase[k] - k * freqPerBin * hop)
      trueFreq = k * freqPerBin + dp / hop
    }
    let shifted = trueFreq * ratio
    let destBin = Math.round(shifted / freqPerBin)
    if (destBin < 0 || destBin > half) { eIn -= e; continue }
    eOut += e
    let r = lobeGain(shifted - destBin * freqPerBin, ctx.N)
    newMag[destBin] += e / (r * r)
    if (mag[k] > peakMag[destBin]) { peakMag[destBin] = mag[k]; newFreq[destBin] = shifted }
  }
  let g = (eOut > 1e-24 && eIn > 1e-24 ? Math.sqrt(eIn / eOut) : 1) * WIN_GAIN
  for (let k = 0; k <= half; k++) if (newMag[k]) newMag[k] = Math.sqrt(newMag[k]) * g
}

// rms(w)/mean(w) for the engine's periodic Hann: sqrt(3/8)/(1/2) = sqrt(3/2).
export const WIN_GAIN = Math.sqrt(1.5)

// Mean overlap-add amplitude of a partial whose intra-frame (bin-grid) and inter-frame
// (true) frequencies differ by `dw` rad/sample: |W(dw)|/W(0) for the engine's periodic
// Hann of length N. The synthesized bin oscillates on the bin grid inside each frame
// while its phase steps at the true frequency across frames, so overlapping frames sum
// slightly incoherently — the classic vocoder scalloping loss, deterministic per bin.
function lobeGain(dw, N) {
  if (!dw) return 1
  let d = (t) => {
    let s = Math.sin(t / 2)
    return Math.abs(s) < 1e-12 ? N : Math.sin(N * t / 2) / s
  }
  let b = PI2 / N
  return Math.abs(0.5 * d(dw) + 0.25 * d(dw - b) + 0.25 * d(dw + b)) / (0.5 * N)
}

// Peak-locked rigid-ROI bin scatter (Laroche-Dolson): each `findPeaks` peak advances its own
// phase at its instantaneous frequency × `ratio`; every other bin rides along rigidly at its
// nearest peak's integer bin-shift, carrying phase relative to that peak (phase coherence
// across the peak's region of influence). Colliding destination bins accumulate in the
// energy domain (Σmag², √ at the end); the phase written is the loudest contributor's —
// same RMS-preserving collision policy as `scatterGated`.
// `reset` skips phase-derivative estimation (first frame, or a caller-detected phase
// discontinuity such as a transient) and uses the analysis phase directly instead of
// integrating. `syn` holds frame-center phases indexed by SOURCE bin, including each
// peak's neighbours so moving peaks inherit coherent history. NaN marks untracked bins.
// `newMag`/`newPhase`/`peakMag` are caller-owned scratch sized `half+1`, zero-filled
// by the caller; `peakDest`/`peakSynPhase` are caller-owned scratch sized `peaks.length`.
export function scatterLocked(mag, phase, prevPhase, reset, peaks, ratio, ctx, syn, newMag, newPhase, peakDest, peakSynPhase, peakMag) {
  let { half, hop, freqPerBin } = ctx
  // Without a resolved partial there is no pitch to move. Preserve DC and
  // short boundary frames rather than dropping their energy; re-seed later peaks.
  if (!peaks.length) {
    for (let k = 0; k <= half; k++) { newMag[k] = mag[k]; newPhase[k] = phase[k]; syn[k] = NaN }
    return
  }
  if (_boost.length <= half) _boost = new Float64Array(half + 1)
  for (let i = 0; i < peaks.length; i++) {
    let k = peaks[i]
    let trueFreq
    if (reset) trueFreq = k * freqPerBin
    else {
      let dp = wrapPhase(phase[k] - prevPhase[k] - k * freqPerBin * hop)
      trueFreq = k * freqPerBin + dp / hop
    }
    let shifted = trueFreq * ratio
    // Shift the lobe by the integer bin count closest to the true frequency delta: the
    // lobe's own frac is preserved, so the intra-/inter-frame frequency mismatch stays
    // within ±half a bin — where the scalloping model below is accurate.
    let destBin = k + Math.round((shifted - trueFreq) / freqPerBin)
    if (destBin < 0 || destBin > half) { peakDest[i] = -1; continue }
    // Center-referenced phase survives destination-bin changes (each bin adds π
    // at N/2). Read all previous phases before updating source-bin history below.
    let newSyn = reset || !Number.isFinite(syn[k]) ? phase[k] + k * Math.PI : wrapPhase(syn[k] + shifted * hop)
    peakDest[i] = destBin
    peakSynPhase[i] = newSyn - destBin * Math.PI
    let r = lobeGain(shifted - (trueFreq + (destBin - k) * freqPerBin), ctx.N)
    _boost[i] = 1 / (r * r)
  }

  // No WIN_GAIN here: the rigid ROI shift carries the whole mainlobe shape to the
  // destination, so the ISTFT frame stays a windowed one — only collision and
  // past-Nyquist bookkeeping need correction. Bins riding a peak whose destination
  // fell outside Nyquist are legitimately lost and excluded from the energy target.
  let eIn = 0, eOut = 0
  for (let k = 0; k <= half; k++) {
    let pi = nearestPeak(peaks, k)
    if (pi < 0) { syn[k] = NaN; continue }
    let destBin = peakDest[pi]
    if (destBin < 0) { syn[k] = NaN; continue }
    let e = mag[k] * mag[k]
    eIn += e
    let pk = peaks[pi]
    syn[k] = wrapPhase(peakSynPhase[pi] + destBin * Math.PI + (phase[k] - phase[pk]) + (k - pk) * Math.PI)
    let dest = destBin + (k - pk)
    if (dest < 0 || dest > half) continue
    eOut += e
    newMag[dest] += e * _boost[pi]
    if (mag[k] > peakMag[dest]) { peakMag[dest] = mag[k]; newPhase[dest] = peakSynPhase[pi] + (phase[k] - phase[pk]) }
  }
  let g = eOut > 1e-24 && eIn > 1e-24 ? Math.sqrt(eIn / eOut) : 1
  for (let k = 0; k <= half; k++) if (newMag[k]) newMag[k] = Math.sqrt(newMag[k]) * g
}

// scatterLocked-internal per-peak scalloping boosts; grown once, reused across frames.
let _boost = new Float64Array(0)

// Variable-ratio resolver for STFT process callbacks. Returns `{ scalar, at }`
// where `scalar` is the ratio at t=0 and `at(frameStart, sampleRate)` resolves
// the ratio for a given frame position.
export function makeFrameRatio(ratio) {
  if (typeof ratio !== 'function') return { scalar: ratio, at: () => ratio }
  let scalar = ratio(0)
  return {
    scalar,
    at(frameStart, sampleRate) {
      let r = ratio(Math.max(0, frameStart) / sampleRate)
      return (!Number.isFinite(r) || r <= 0) ? (scalar || 1) : r
    }
  }
}

// ── Time-stretch phase locking (Laroche & Dolson 1999) ──
// Prominence-based peak mask tuned for stretch: keeps chord partials, rejects
// leakage shoulders. Internal to lockPhase; distinct from the shift-side findPeaks.
let _peakMask = new Uint8Array(0)
function peakMask(mag, half) {
  if (_peakMask.length !== half + 1) _peakMask = new Uint8Array(half + 1)
  let peaks = _peakMask
  peaks.fill(0)
  if (half <= 1) {
    peaks[0] = 1
    if (half === 1) peaks[1] = 1
    return peaks
  }

  let maxMag = 0
  for (let k = 0; k <= half; k++) if (mag[k] > maxMag) maxMag = mag[k]

  let minMag = Math.max(1e-8, maxMag * 0.015)
  let minProm = Math.max(1e-9, maxMag * 0.003)
  let lastPeak = -2, lastPeakMag = 0

  for (let k = 1; k < half; k++) {
    let value = mag[k]
    if (value < minMag || value < mag[k - 1] || value < mag[k + 1]) continue

    let shoulder = Math.max(mag[k - 1], mag[k + 1], k > 1 ? mag[k - 2] : 0, k + 2 <= half ? mag[k + 2] : 0)
    if (value - shoulder < minProm && value < maxMag * 0.1) continue

    if (k - lastPeak <= 1) {
      if (value > lastPeakMag) { peaks[lastPeak] = 0; peaks[k] = 1; lastPeak = k; lastPeakMag = value }
      continue
    }
    peaks[k] = 1; lastPeak = k; lastPeakMag = value
  }

  let found = false
  for (let k = 0; k <= half; k++) if (peaks[k]) { found = true; break }
  if (!found) {
    let best = 0
    for (let k = 1; k <= half; k++) if (mag[k] > mag[best]) best = k
    peaks[best] = 1
  }
  return peaks
}

// Which bins ride which peak: `owner[k]` = the peak whose rotation bin k takes (its region of
// influence reaches halfway to the neighbouring peaks), or -1 for a free bin: a peak itself, or a
// bin under 3% of its peak's magnitude, which keeps its own phase advance. Returns `owner`.
let _peakBins = new Int32Array(0)
export function lockMap(mag, half, owner) {
  let peaks = peakMask(mag, half)
  if (_peakBins.length < half + 1) _peakBins = new Int32Array(half + 1)
  let peakBins = _peakBins, nBins = 0
  for (let k = 0; k <= half; k++) if (peaks[k]) peakBins[nBins++] = k
  owner.fill(-1, 0, half + 1)
  for (let i = 0; i < nBins; i++) {
    let pk = peakBins[i]
    let start = i === 0 ? 0 : Math.floor((peakBins[i - 1] + pk) * 0.5) + 1
    let end = i === nBins - 1 ? half : Math.floor((pk + peakBins[i + 1]) * 0.5)
    let lockFloor = Math.max(1e-10, mag[pk] * 0.03)
    for (let k = start; k <= end; k++) if (k !== pk && !(mag[k] < lockFloor)) owner[k] = pk
  }
  return owner
}

// Lock non-peak bin phases to nearest peak's rotation, in place on `propPhase`.
let _owner = new Int32Array(0)
export function lockPhase(phase, propPhase, mag, half) {
  if (_owner.length < half + 1) _owner = new Int32Array(half + 1)
  let owner = lockMap(mag, half, _owner)
  for (let k = 0; k <= half; k++) {
    let pk = owner[k]
    if (pk >= 0) propPhase[k] = phase[k] + (propPhase[pk] - phase[pk])
  }
}

// ── Complex-domain phase locking ──
// The propagate-then-lockPhase vocoder step on complex bins (pair with fourier-transform/stft's `complex` frames).
// A bin riding a peak takes the peak's rotation R = Y(p)·conj(X(p))/|X(p)|², Y(k) = X(k)·R: one complex multiply,
// the identity locking of Laroche & Dolson. Only free bins (the peaks, and bins too quiet to lock) advance at
// their instantaneous frequency, one atan2 and one sin/cos each. State keeps unit synthesis phasors and the previous
// analysis bins, so no bin needs its phase as an angle.

/** Frame-to-frame state over half+1 bins: `mag` (the caller fills this frame's magnitudes), the lock map, the
 *  previous analysis bins and magnitudes, the unit synthesis phasors. */
export function lockState(half) {
  let n = half + 1, f = () => new Float64Array(n)
  return { mag: f(), owner: new Int32Array(n), xr: f(), xi: f(), xm: f(), ur: f(), ui: f() }
}

/** One vocoder frame in place: `re`/`im` hold the analysis bins, then the synthesis bins. `reset` (first frame, a
 *  transient) restarts each bin's synthesis phase at its analysis phase; otherwise free bins advance by their
 *  instantaneous frequency over the hops and locked bins take their peak's rotation. */
export function lockAdvance(re, im, st, reset, anaHop, synHop, freqPerBin, half) {
  let { mag, owner, xr, xi, xm, ur, ui } = st, n = half + 1
  if (reset) for (let k = 0; k < n; k++) {
    let m = mag[k]
    if (m) { ur[k] = re[k] / m; ui[k] = im[k] / m; continue }
    let a = Math.atan2(im[k], re[k])   // a zero bin: the direction its zeros' signs give atan2
    ur[k] = Math.cos(a); ui[k] = Math.sin(a)
  } else {
    lockMap(mag, half, owner)
    for (let k = 0; k < n; k++) {
      if (owner[k] >= 0) continue
      // phase advance since the previous frame, arg(X·conj(Xprev)); a zero bin has no product angle: its own angles
      let d = mag[k] && xm[k] ? Math.atan2(im[k] * xr[k] - re[k] * xi[k], re[k] * xr[k] + im[k] * xi[k]) : Math.atan2(im[k], re[k]) - Math.atan2(xi[k], xr[k])
      let dp = wrapPhase(d - k * freqPerBin * anaHop)
      let a = (k * freqPerBin + dp / anaHop) * synHop, c = Math.cos(a), s = Math.sin(a), u = ur[k], v = ui[k]
      ur[k] = u * c - v * s
      ui[k] = u * s + v * c
    }
    for (let k = 0, p = -1, rr = 0, ri = 0; k < n; k++) {
      let pk = owner[k]
      if (pk < 0) continue
      if (pk !== p) { p = pk; let m = mag[p]; rr = (ur[p] * re[p] + ui[p] * im[p]) / m; ri = (ui[p] * re[p] - ur[p] * im[p]) / m }
      let q = 1 / mag[k]
      ur[k] = (re[k] * rr - im[k] * ri) * q
      ui[k] = (re[k] * ri + im[k] * rr) * q
    }
  }
  for (let k = 0; k < n; k++) {
    xr[k] = re[k]; xi[k] = im[k]; xm[k] = mag[k]
    re[k] = mag[k] * ur[k]; im[k] = mag[k] * ui[k]
  }
}
