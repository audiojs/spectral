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
// The peaks in order into `list`, their count returned (an adjacent louder one takes the last one's place).
function peakList(mag, half, list) {
  if (half <= 1) {
    list[0] = 0
    if (half === 1) { list[1] = 1; return 2 }
    return 1
  }

  let maxMag = 0
  for (let k = 0; k <= half; k++) if (mag[k] > maxMag) maxMag = mag[k]

  let minMag = Math.max(1e-8, maxMag * 0.015)
  let minProm = Math.max(1e-9, maxMag * 0.003)
  let lastPeak = -2, lastPeakMag = 0, n = 0

  for (let k = 1; k < half; k++) {
    let value = mag[k]
    if (value < minMag || value < mag[k - 1] || value < mag[k + 1]) continue

    let shoulder = Math.max(mag[k - 1], mag[k + 1], k > 1 ? mag[k - 2] : 0, k + 2 <= half ? mag[k + 2] : 0)
    if (value - shoulder < minProm && value < maxMag * 0.1) continue

    if (k - lastPeak <= 1) {
      if (value > lastPeakMag) { list[n - 1] = k; lastPeak = k; lastPeakMag = value }
      continue
    }
    list[n++] = k; lastPeak = k; lastPeakMag = value
  }

  if (!n) {
    let best = 0
    for (let k = 1; k <= half; k++) if (mag[k] > mag[best]) best = k
    list[n++] = best
  }
  return n
}

// Which bins ride which peak: `owner[k]` = the peak whose rotation bin k takes (its region of
// influence reaches halfway to the neighbouring peaks), or -1 for a free bin: a peak itself, or a
// bin under 3% of its peak's magnitude, which keeps its own phase advance. Returns `owner`.
// The regions tile 0…half (each starts where the one before ends), so every bin is set once; the free bins and the
// locked ones are also listed in order (_free, _locked; _nf, _nl of them), for lockAdvance to run through.
let _peakBins = new Int32Array(0), _free = new Int32Array(0), _locked = new Int32Array(0), _nf = 0, _nl = 0
export function lockMap(mag, half, owner) {
  if (_peakBins.length < half + 1) { _peakBins = new Int32Array(half + 1); _free = new Int32Array(half + 1); _locked = new Int32Array(half + 1) }
  let peakBins = _peakBins, nBins = peakList(mag, half, peakBins), free = _free, locked = _locked, nf = 0, nl = 0
  for (let i = 0; i < nBins; i++) {
    let pk = peakBins[i]
    let start = i === 0 ? 0 : Math.floor((peakBins[i - 1] + pk) * 0.5) + 1
    let end = i === nBins - 1 ? half : Math.floor((pk + peakBins[i + 1]) * 0.5)
    let lockFloor = Math.max(1e-10, mag[pk] * 0.03)
    for (let k = start; k <= end; k++) {
      if (k !== pk && !(mag[k] < lockFloor)) { owner[k] = pk; locked[nl++] = k }
      else { owner[k] = -1; free[nf++] = k }
    }
  }
  _nf = nf; _nl = nl
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
// With whole hops (the stft engines' frames start at whole samples), a free bin's advance splits in two: the bin
// centre's, ω_k·hop with ω_k = 2πk/N, a whole number of N-th turns, read from a table (no reduction of an angle of a
// thousand radians, whose rounding alone is 1e-13); and its deviation from it, within a turn: the measured product
// rotated back by the centre's analysis advance gives it as one atan2, already wrapped, and its synthesis share is
// rotated on by sin and cos of under an eighth of a turn, the quarter turns again from the table. The trigonometry is
// fdlibm's polynomials (Sun Microsystems 1993: s_atan.c, k_sin.c, k_cos.c) in basic float ops, within 2 ulp of the
// exact values and the same on every platform (Math's differ by platform in the last bit).

const S1 = -1.66666666666666324348e-01, S2 = 8.33333333332248946124e-03, S3 = -1.98412698298579493134e-04,
  S4 = 2.75573137070700676789e-06, S5 = -2.50507602534068634195e-08, S6 = 1.58969099521155010221e-10
const C1 = 4.16666666666666019037e-02, C2 = -1.38888888888741095749e-03, C3 = 2.48015872894767294178e-05,
  C4 = -2.75573143513906633035e-07, C5 = 2.08757232129817482790e-09, C6 = -1.13596475577881948265e-11
const T0 = 3.33333333333329318027e-01, T1 = -1.99999999998764832476e-01, T2 = 1.42857142725034663711e-01,
  T3 = -1.11111104054623557880e-01, T4 = 9.09088713343650656196e-02, T5 = -7.69187620504482999495e-02,
  T6 = 6.66107313738753120669e-02, T7 = -5.83357013379057348645e-02, T8 = 4.97687799461593236017e-02,
  T9 = -3.65315727442169155270e-02, T10 = 1.62858201153657823623e-02
const PIO2_1 = 1.57079632673412561417e+00, PIO2_1T = 6.07710050650619224932e-11, INV_PIO2 = 6.36619772367581382433e-01
const PIO4_LO = 3.06161699786838301793e-17, PIO2_LO = 6.123233995736766e-17, PI_LO = 1.2246467991473532e-16, TAN8 = Math.SQRT2 - 1

/** atan2(y, x) of a nonzero (x, y): the octant by min/max, atan on |u| ≤ tan(π/8) by fdlibm's polynomial (its range is
 *  |x| < 7/16), the octant's angle added back with its low part */
function atan2(y, x) {
  let ax = Math.abs(x), ay = Math.abs(y), t = Math.min(ax, ay) / Math.max(ax, ay), big = t > TAN8
  let u = big ? (t - 1) / (t + 1) : t, z = u * u, w = z * z
  let p = u * (z * (T0 + w * (T2 + w * (T4 + w * (T6 + w * (T8 + w * T10))))) + w * (T1 + w * (T3 + w * (T5 + w * (T7 + w * T9)))))
  let r = big ? Math.PI / 4 - ((p - PIO4_LO) - u) : u - p
  if (ay > ax) r = Math.PI / 2 - (r - PIO2_LO)
  if (x < 0) r = Math.PI - (r - PI_LO)
  return y < 0 ? -r : r
}
/** sin and cos of |x| ≤ π/4 */
const sink = x => { let z = x * x; return x + z * x * (S1 + z * (S2 + z * (S3 + z * (S4 + z * (S5 + z * S6))))) }
const cosk = x => { let z = x * x, r = z * (C1 + z * (C2 + z * (C3 + z * (C4 + z * (C5 + z * C6))))), h = 0.5 * z, w = 1 - h; return w + (((1 - w) - h) + z * r) }

/** e^(2πi·m/N) for m < N: the first eighth by the kernels, the rest by its symmetries (exact) */
const TURNS = new Map()
function turns(N) {
  let t = TURNS.get(N)
  if (t) return t
  let c = new Float64Array(N), s = new Float64Array(N), q = N >> 2
  let put = (m, x, y) => { for (let j = 0; j < 4; j++, [x, y] = [-y, x]) { c[(m + j * q) % N] = x; s[(m + j * q) % N] = y } }
  for (let m = 0; m <= N >> 3; m++) { let a = 2 * Math.PI * m / N, x = cosk(a), y = sink(a); put(m, x, y); put(q - m, y, x) }
  TURNS.set(N, t = { c, s })
  return t
}

// pvoc.wat, compiled by wat2wasm (3,824 bytes): lockAdvance's frame with whole hops on a state's own memory, two bins to a
// float64x2 where they run alike, the same float ops in the same order as the JS below (Math's atan2, cos and sin
// imported for the bins that take them)
const WASM = new Uint8Array([
  0,97,115,109,1,0,0,0,1,46,7,96,2,124,124,1,124,96,1,124,1,124,96,1,124,1,127,96,1,123,1,123,96,2,123,123,1,123,96,7,
  127,127,127,127,127,127,127,0,96,5,127,127,127,124,127,0,2,27,3,1,109,5,97,116,97,110,50,0,0,1,109,3,99,111,115,0,1,1,
  109,3,115,105,110,0,1,3,9,8,2,1,3,3,3,4,5,6,5,3,1,0,1,6,16,3,127,1,65,0,11,127,1,65,0,11,127,1,65,0,11,7,17,2,6,109,
  101,109,111,114,121,2,0,4,108,111,99,107,0,10,10,227,28,8,20,0,32,0,68,0,0,0,0,0,0,0,0,98,32,0,32,0,97,113,11,38,0,32,
  0,32,0,68,24,45,68,84,251,33,25,64,163,68,0,0,0,0,0,0,224,63,160,156,68,24,45,68,84,251,33,25,64,162,161,11,88,1,1,
  123,32,0,253,117,33,1,32,1,32,0,32,1,253,241,1,253,12,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,253,76,253,12,0,0,0,0,0,0,
  240,63,0,0,0,0,0,0,240,63,253,78,253,240,1,33,1,32,0,32,1,253,12,0,0,0,0,0,0,0,128,0,0,0,0,0,0,0,128,253,82,11,176,1,
  1,1,123,32,0,32,0,253,242,1,33,1,32,0,32,1,32,0,253,242,1,253,12,73,85,85,85,85,85,197,191,73,85,85,85,85,85,197,191,
  32,1,253,12,166,248,16,17,17,17,129,63,166,248,16,17,17,17,129,63,32,1,253,12,213,97,193,25,160,1,42,191,213,97,193,
  25,160,1,42,191,32,1,253,12,125,254,177,87,227,29,199,62,125,254,177,87,227,29,199,62,32,1,253,12,235,156,43,138,230,
  229,90,190,235,156,43,138,230,229,90,190,32,1,253,12,124,213,207,90,58,217,229,61,124,213,207,90,58,217,229,61,253,
  242,1,253,240,1,253,242,1,253,240,1,253,242,1,253,240,1,253,242,1,253,240,1,253,242,1,253,240,1,253,242,1,253,240,1,
  11,133,2,1,4,123,32,0,32,0,253,242,1,33,1,32,1,253,12,76,85,85,85,85,85,165,63,76,85,85,85,85,85,165,63,32,1,253,12,
  119,81,193,22,108,193,86,191,119,81,193,22,108,193,86,191,32,1,253,12,144,21,203,25,160,1,250,62,144,21,203,25,160,1,
  250,62,32,1,253,12,173,82,156,128,79,126,146,190,173,82,156,128,79,126,146,190,32,1,253,12,196,177,180,189,158,238,33,
  62,196,177,180,189,158,238,33,62,32,1,253,12,212,56,136,190,233,250,168,189,212,56,136,190,233,250,168,189,253,242,1,
  253,240,1,253,242,1,253,240,1,253,242,1,253,240,1,253,242,1,253,240,1,253,242,1,253,240,1,253,242,1,33,2,253,12,0,0,0,
  0,0,0,224,63,0,0,0,0,0,0,224,63,32,1,253,242,1,33,3,253,12,0,0,0,0,0,0,240,63,0,0,0,0,0,0,240,63,32,3,253,241,1,33,4,
  32,4,253,12,0,0,0,0,0,0,240,63,0,0,0,0,0,0,240,63,32,4,253,241,1,32,3,253,241,1,32,1,32,2,253,242,1,253,240,1,253,240,
  1,11,135,5,1,9,123,32,1,253,236,1,33,2,32,0,253,236,1,33,3,32,2,32,3,253,244,1,32,2,32,3,253,245,1,253,243,1,33,4,32,
  4,253,12,52,239,252,153,121,130,218,63,52,239,252,153,121,130,218,63,253,74,33,5,32,4,253,12,0,0,0,0,0,0,240,63,0,0,0,
  0,0,0,240,63,253,241,1,32,4,253,12,0,0,0,0,0,0,240,63,0,0,0,0,0,0,240,63,253,240,1,253,243,1,32,4,32,5,253,82,33,6,32,
  6,32,6,253,242,1,33,7,32,7,32,7,253,242,1,33,8,32,6,32,7,253,12,13,85,85,85,85,85,213,63,13,85,85,85,85,85,213,63,32,
  8,253,12,255,131,0,146,36,73,194,63,255,131,0,146,36,73,194,63,32,8,253,12,110,32,76,197,205,69,183,63,110,32,76,197,
  205,69,183,63,32,8,253,12,81,61,208,160,102,13,177,63,81,61,208,160,102,13,177,63,32,8,253,12,235,13,118,36,75,123,
  169,63,235,13,118,36,75,123,169,63,32,8,253,12,17,218,34,227,58,173,144,63,17,218,34,227,58,173,144,63,253,242,1,253,
  240,1,253,242,1,253,240,1,253,242,1,253,240,1,253,242,1,253,240,1,253,242,1,253,240,1,253,242,1,32,8,253,12,196,235,
  152,153,153,153,201,191,196,235,152,153,153,153,201,191,32,8,253,12,113,22,35,254,198,113,188,191,113,22,35,254,198,
  113,188,191,32,8,253,12,109,154,116,175,242,176,179,191,109,154,116,175,242,176,179,191,32,8,253,12,154,253,222,82,45,
  222,173,191,154,253,222,82,45,222,173,191,32,8,253,12,47,108,106,44,68,180,162,191,47,108,106,44,68,180,162,191,253,
  242,1,253,240,1,253,242,1,253,240,1,253,242,1,253,240,1,253,242,1,253,240,1,253,242,1,253,240,1,253,242,1,33,9,253,12,
  24,45,68,84,251,33,233,63,24,45,68,84,251,33,233,63,32,9,253,12,7,92,20,51,38,166,129,60,7,92,20,51,38,166,129,60,253,
  241,1,32,6,253,241,1,253,241,1,32,6,32,9,253,241,1,32,5,253,82,33,10,253,12,24,45,68,84,251,33,249,63,24,45,68,84,251,
  33,249,63,32,10,253,12,7,92,20,51,38,166,145,60,7,92,20,51,38,166,145,60,253,241,1,253,241,1,32,10,32,3,32,2,253,74,
  253,82,33,10,253,12,24,45,68,84,251,33,9,64,24,45,68,84,251,33,9,64,32,10,253,12,7,92,20,51,38,166,161,60,7,92,20,51,
  38,166,161,60,253,241,1,253,241,1,32,10,32,1,253,12,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,253,73,253,82,33,10,32,10,253,237,
  1,32,10,32,0,253,12,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,253,73,253,82,11,230,5,2,13,127,8,124,65,0,33,7,3,64,32,1,32,7,65,
  3,116,106,43,3,0,33,20,32,20,32,22,100,4,64,32,20,33,22,11,32,7,65,1,106,34,7,32,0,76,13,0,11,68,58,140,48,226,142,
  121,69,62,32,22,68,184,30,133,235,81,184,142,63,162,165,33,23,68,149,214,38,232,11,46,17,62,32,22,68,250,126,106,188,
  116,147,104,63,162,165,33,24,65,126,33,13,65,1,33,7,2,64,3,64,32,7,32,0,78,13,1,2,64,32,1,32,7,65,3,116,106,33,17,32,
  17,43,3,0,33,20,32,20,32,23,99,13,0,32,20,32,17,65,8,107,43,3,0,99,13,0,32,20,32,17,43,3,8,99,13,0,32,17,65,8,107,43,
  3,0,32,17,43,3,8,165,32,7,65,1,74,4,124,32,17,65,16,107,43,3,0,5,68,0,0,0,0,0,0,0,0,11,32,7,65,2,106,32,0,76,4,124,32,
  17,43,3,16,5,68,0,0,0,0,0,0,0,0,11,165,165,33,26,32,20,32,26,161,32,24,99,32,20,32,22,68,154,153,153,153,153,153,185,
  63,162,99,113,13,0,32,7,32,13,107,65,1,76,4,64,32,20,32,25,100,4,64,32,6,32,8,65,1,107,65,2,116,106,32,7,54,2,0,32,7,
  33,13,32,20,33,25,11,12,1,11,32,6,32,8,65,2,116,106,32,7,54,2,0,32,8,65,1,106,33,8,32,7,33,13,32,20,33,25,11,32,7,65,
  1,106,33,7,12,0,11,11,32,8,69,4,64,65,1,33,7,3,64,32,1,32,7,65,3,116,106,43,3,0,32,1,32,16,65,3,116,106,43,3,0,100,4,
  64,32,7,33,16,11,32,7,65,1,106,34,7,32,0,76,13,0,11,32,6,32,16,54,2,0,65,1,33,8,11,65,0,33,9,3,64,32,6,32,9,65,2,116,
  106,40,2,0,33,10,32,9,4,127,32,6,32,9,65,1,107,65,2,116,106,40,2,0,32,10,106,65,1,117,65,1,106,5,65,0,11,33,11,32,9,
  32,8,65,1,107,70,4,127,32,0,5,32,10,32,6,32,9,65,1,106,65,2,116,106,40,2,0,106,65,1,117,11,33,12,68,187,189,215,217,
  223,124,219,61,32,1,32,10,65,3,116,106,43,3,0,68,184,30,133,235,81,184,158,63,162,165,33,27,32,11,33,7,2,64,3,64,32,7,
  32,12,74,13,1,32,7,65,3,116,33,17,32,1,32,17,106,43,3,0,33,20,32,7,32,10,71,32,20,32,27,99,69,113,33,18,32,3,32,7,65,
  2,116,106,32,10,65,127,32,18,27,54,2,0,32,7,65,0,71,32,7,32,0,72,113,32,20,68,0,0,0,0,0,0,0,0,98,32,20,32,20,97,113,
  32,2,32,17,106,43,3,0,34,21,68,0,0,0,0,0,0,0,0,98,32,21,32,21,97,113,113,113,33,19,32,4,32,14,65,2,116,106,32,7,54,2,
  0,32,5,32,15,65,2,116,106,32,7,54,2,0,32,14,32,18,69,32,19,113,106,33,14,32,15,32,18,69,32,19,69,113,106,33,15,32,7,
  65,1,106,33,7,12,0,11,11,32,9,65,1,106,34,9,32,8,72,13,0,11,32,14,36,0,32,15,36,1,32,8,36,2,11,161,13,3,37,127,24,124,
  14,123,32,0,65,3,106,65,126,113,33,5,32,5,65,3,116,33,6,32,0,65,1,116,33,7,32,7,65,1,107,33,8,32,6,33,11,32,6,65,2,
  108,33,12,32,6,65,3,108,33,13,32,6,65,4,108,33,14,32,6,65,5,108,33,15,32,6,65,6,108,33,16,32,6,65,7,108,33,17,32,6,65,
  8,108,33,18,32,6,65,9,108,33,19,32,6,65,10,108,33,20,32,6,65,11,108,33,21,32,6,65,12,108,33,22,32,6,65,13,108,33,23,
  32,6,65,14,108,33,24,32,6,65,15,108,33,25,32,25,32,7,65,3,116,106,33,26,32,26,32,7,65,3,116,106,33,27,32,27,32,5,65,2,
  116,106,33,28,32,28,32,5,65,2,116,106,33,29,32,29,32,5,65,2,116,106,33,30,32,0,65,1,106,65,3,116,33,9,32,4,4,64,65,0,
  33,33,3,64,32,33,253,0,4,0,33,71,32,11,32,33,106,253,0,4,0,33,72,32,12,32,33,106,32,71,32,71,253,242,1,32,72,32,72,
  253,242,1,253,240,1,253,239,1,253,11,4,0,32,33,65,16,106,34,33,32,9,72,13,0,11,11,32,0,32,12,32,15,32,27,32,28,32,29,
  32,30,16,9,35,0,33,37,35,1,33,38,35,2,33,39,32,1,183,33,42,32,2,183,33,43,65,0,33,31,2,64,3,64,32,31,32,37,78,13,1,32,
  28,32,31,65,2,116,106,40,2,0,33,32,32,32,65,3,116,33,33,32,33,43,3,0,33,48,32,11,32,33,106,43,3,0,33,49,32,13,32,33,
  106,43,3,0,33,50,32,14,32,33,106,43,3,0,33,51,32,48,32,50,162,32,49,32,51,162,160,33,54,32,49,32,50,162,32,48,32,51,
  162,161,33,55,32,32,32,1,108,32,8,113,65,3,116,33,34,32,25,32,34,106,43,3,0,33,56,32,26,32,34,106,43,3,0,33,57,32,31,
  65,3,116,33,33,32,18,32,33,106,32,55,32,56,162,32,54,32,57,162,161,57,3,0,32,19,32,33,106,32,54,32,56,162,32,55,32,57,
  162,160,57,3,0,32,31,65,1,106,33,31,12,0,11,11,32,42,253,20,33,66,32,43,253,20,33,67,65,0,33,33,32,37,65,3,116,33,34,
  2,64,3,64,32,33,32,34,78,13,1,32,18,32,33,106,253,0,4,0,32,19,32,33,106,253,0,4,0,16,8,32,66,253,243,1,32,67,253,242,
  1,33,68,32,68,253,12,131,200,201,109,48,95,228,63,131,200,201,109,48,95,228,63,253,242,1,16,5,33,69,32,68,32,69,253,
  12,0,0,64,84,251,33,249,63,0,0,64,84,251,33,249,63,253,242,1,253,241,1,32,69,253,12,49,99,98,26,97,180,208,61,49,99,
  98,26,97,180,208,61,253,242,1,253,241,1,33,70,32,20,32,33,106,32,69,253,11,4,0,32,21,32,33,106,32,70,16,7,253,11,4,0,
  32,22,32,33,106,32,70,16,6,253,11,4,0,32,33,65,16,106,33,33,12,0,11,11,65,0,33,31,2,64,3,64,32,31,32,37,78,13,1,32,28,
  32,31,65,2,116,106,40,2,0,33,32,32,31,65,3,116,33,33,32,32,32,2,108,32,20,32,33,106,43,3,0,252,2,32,7,65,2,117,108,
  106,32,8,113,65,3,116,33,35,32,21,32,33,106,43,3,0,33,58,32,22,32,33,106,43,3,0,33,59,32,25,32,35,106,43,3,0,33,60,32,
  26,32,35,106,43,3,0,33,61,32,60,32,58,162,32,61,32,59,162,161,33,46,32,60,32,59,162,32,61,32,58,162,160,33,47,32,32,
  65,3,116,33,33,32,16,32,33,106,43,3,0,33,44,32,17,32,33,106,43,3,0,33,45,32,16,32,33,106,32,44,32,46,162,32,45,32,47,
  162,161,57,3,0,32,17,32,33,106,32,44,32,47,162,32,45,32,46,162,160,57,3,0,32,31,65,1,106,33,31,12,0,11,11,65,0,33,31,
  2,64,3,64,32,31,32,38,78,13,1,32,29,32,31,65,2,116,106,40,2,0,33,32,32,32,65,3,116,33,33,32,33,43,3,0,33,48,32,11,32,
  33,106,43,3,0,33,49,32,13,32,33,106,43,3,0,33,50,32,14,32,33,106,43,3,0,33,51,32,12,32,33,106,43,3,0,16,3,32,15,32,33,
  106,43,3,0,16,3,113,4,124,32,49,32,50,162,32,48,32,51,162,161,32,48,32,50,162,32,49,32,51,162,160,16,0,5,32,49,32,48,
  16,0,32,51,32,50,16,0,161,11,33,62,32,62,32,32,183,32,3,162,32,42,162,161,16,4,33,62,32,32,183,32,3,162,32,62,32,42,
  163,160,32,43,162,33,63,32,63,16,1,33,46,32,63,16,2,33,47,32,16,32,33,106,43,3,0,33,44,32,17,32,33,106,43,3,0,33,45,
  32,16,32,33,106,32,44,32,46,162,32,45,32,47,162,161,57,3,0,32,17,32,33,106,32,44,32,47,162,32,45,32,46,162,160,57,3,0,
  32,31,65,1,106,33,31,12,0,11,11,65,0,33,31,3,64,32,30,32,31,65,2,116,106,40,2,0,33,36,32,31,4,127,32,30,32,31,65,1,
  107,65,2,116,106,40,2,0,32,36,106,65,1,117,65,1,106,5,65,0,11,33,40,32,31,32,39,65,1,107,70,4,127,32,0,5,32,36,32,30,
  32,31,65,1,106,65,2,116,106,40,2,0,106,65,1,117,11,33,41,32,36,65,3,116,33,33,32,12,32,33,106,43,3,0,33,52,32,16,32,
  33,106,43,3,0,33,44,32,17,32,33,106,43,3,0,33,45,32,33,43,3,0,33,48,32,11,32,33,106,43,3,0,33,49,32,44,32,48,162,32,
  45,32,49,162,160,32,52,163,33,64,32,45,32,48,162,32,44,32,49,162,161,32,52,163,33,65,32,64,253,20,33,74,32,65,253,20,
  33,75,32,40,65,3,116,33,33,32,41,65,1,106,65,3,116,33,34,3,64,32,23,32,33,106,32,74,253,11,4,0,32,24,32,33,106,32,75,
  253,11,4,0,32,33,65,16,106,34,33,32,34,72,13,0,11,32,31,65,1,106,34,31,32,39,72,13,0,11,65,0,33,33,3,64,32,33,253,0,4,
  0,33,71,32,11,32,33,106,253,0,4,0,33,72,32,12,32,33,106,253,0,4,0,33,73,32,23,32,33,106,253,0,4,0,33,74,32,24,32,33,
  106,253,0,4,0,33,75,32,27,32,33,65,1,118,106,253,93,3,0,253,199,1,253,12,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,253,219,1,33,
  77,253,12,0,0,0,0,0,0,240,63,0,0,0,0,0,0,240,63,32,73,253,243,1,33,76,32,71,32,74,253,242,1,32,72,32,75,253,242,1,253,
  241,1,32,76,253,242,1,32,16,32,33,106,253,0,4,0,32,77,253,82,33,78,32,71,32,75,253,242,1,32,72,32,74,253,242,1,253,
  240,1,32,76,253,242,1,32,17,32,33,106,253,0,4,0,32,77,253,82,33,79,32,16,32,33,106,32,78,253,11,4,0,32,17,32,33,106,
  32,79,253,11,4,0,32,13,32,33,106,32,71,253,11,4,0,32,14,32,33,106,32,72,253,11,4,0,32,15,32,33,106,32,73,253,11,4,0,
  32,33,32,73,32,78,253,242,1,253,11,4,0,32,11,32,33,106,32,73,32,79,253,242,1,253,11,4,0,32,33,65,16,106,34,33,32,9,72,
  13,0,11,11,
])

let mod
/** A state's arrays on a WebAssembly memory of its own, with the kernel that advances them (`k`, not enumerated);
 *  null where WebAssembly does not run */
function wasmState(half) {
  if (mod === undefined) try { mod = new WebAssembly.Module(WASM) } catch { mod = null }
  if (!mod || half < 2) return null
  let n = half + 1, n2 = (n + 2) & ~1, B = 8 * n2, N = 2 * half, need = 15 * B + 16 * N + 16 * n2
  let { exports: x } = new WebAssembly.Instance(mod, { m: { atan2: Math.atan2, cos: Math.cos, sin: Math.sin } })
  if (need > x.memory.buffer.byteLength) x.memory.grow(Math.ceil((need - x.memory.buffer.byteLength) / 65536))
  let buf = x.memory.buffer, f = i => new Float64Array(buf, i * B, n)
  if ((N & (N - 1)) === 0) { let { c, s } = turns(N); new Float64Array(buf, 15 * B, N).set(c); new Float64Array(buf, 15 * B + 8 * N, N).set(s) }
  let st = { mag: f(2), owner: new Int32Array(buf, 15 * B + 16 * N, n), xr: f(3), xi: f(4), xm: f(5), ur: f(6), ui: f(7) }
  return Object.defineProperty(st, 'k', { value: { half, re: f(0), im: f(1), lock: x.lock } })
}

/** Frame-to-frame state over half+1 bins: `mag` (the caller fills this frame's magnitudes), the lock map, the
 *  previous analysis bins and magnitudes, the unit synthesis phasors. */
export function lockState(half) {
  let w = wasmState(half)
  if (w) return w
  let n = half + 1, f = () => new Float64Array(n)
  return { mag: f(), owner: new Int32Array(n), xr: f(), xi: f(), xm: f(), ur: f(), ui: f() }
}

// whole hops on a power-of-two frame: the turns' table serves every free bin (the kernel's frames, and the JS loop's)
const whole = (N, A, S, freqPerBin) => (N & (N - 1)) === 0 && Number.isInteger(A) && Number.isInteger(S) && freqPerBin === PI2 / N
const put = (to, from, n) => to.set(from.length === n ? from : from.subarray(0, n))

/** The frame through the state's kernel, on the bins copied into its memory and back, its magnitudes first with `mags`;
 *  false where it does not run it (no kernel, a reset, other hops, its i32 products k·hop out of range) */
function kernelStep(re, im, st, reset, anaHop, synHop, freqPerBin, half, mags) {
  let w = st.k, n = half + 1
  if (!w || reset || w.half !== half || !whole(2 * half, anaHop, synHop, freqPerBin) || !(Math.abs(half * anaHop) < 2 ** 31 && Math.abs(half * synHop) < 2 ** 31)) return false
  put(w.re, re, n); put(w.im, im, n)
  w.lock(half, anaHop, synHop, freqPerBin, mags)
  re.set(w.re); im.set(w.im)
  return true
}

/** lockAdvance with this frame's magnitudes taken from its bins, √(re² + im²) into `st.mag`, first */
export function lockStep(re, im, st, reset, anaHop, synHop, freqPerBin, half) {
  if (kernelStep(re, im, st, reset, anaHop, synHop, freqPerBin, half, 1)) return
  let { mag } = st
  for (let k = 0; k <= half; k++) mag[k] = Math.sqrt(re[k] * re[k] + im[k] * im[k])
  lockAdvance(re, im, st, reset, anaHop, synHop, freqPerBin, half)
}

/** One vocoder frame in place: `re`/`im` hold the analysis bins, then the synthesis bins. `reset` (first frame, a
 *  transient) restarts each bin's synthesis phase at its analysis phase; otherwise free bins advance by their
 *  instantaneous frequency over the hops and locked bins take their peak's rotation. */
export function lockAdvance(re, im, st, reset, anaHop, synHop, freqPerBin, half) {
  if (kernelStep(re, im, st, reset, anaHop, synHop, freqPerBin, half, 0)) return
  let { mag, owner, xr, xi, xm, ur, ui } = st, n = half + 1
  if (reset) for (let k = 0; k < n; k++) {
    let m = mag[k]
    if (m) { ur[k] = re[k] / m; ui[k] = im[k] / m; continue }
    let a = Math.atan2(im[k], re[k])   // a zero bin: the direction its zeros' signs give atan2
    ur[k] = Math.cos(a); ui[k] = Math.sin(a)
  } else {
    lockMap(mag, half, owner)
    let N = 2 * half, M = N - 1, A = anaHop, S = synHop, free = _free, nf = _nf, locked = _locked, nl = _nl
    let tab = whole(N, A, S, freqPerBin) ? turns(N) : null, tc = tab?.c, ts = tab?.s
    for (let i = 0; i < nf; i++) {
      let k = free[i], c, s, u = ur[k], v = ui[k]
      if (tab && k && k < half && mag[k] && xm[k]) {
        // X·conj(Xprev) rotated back by the centre's ω_k·A: e^(i·dp), the deviation, wrapped; its synthesis share th
        let dr = re[k] * xr[k] + im[k] * xi[k], di = im[k] * xr[k] - re[k] * xi[k], m = (k * A) & M, wc = tc[m], ws = ts[m]
        let th = atan2(di * wc - dr * ws, dr * wc + di * ws) / A * S, q = Math.round(th * INV_PIO2), r = th - q * PIO2_1 - q * PIO2_1T
        // on by the centre's ω_k·S and th: their quarter turns from the table, the rest by the kernels
        let j = (k * S + q * (N >> 2)) & M, cr = cosk(r), sr = sink(r)
        c = tc[j] * cr - ts[j] * sr; s = tc[j] * sr + ts[j] * cr
      } else {
        // phase advance since the previous frame, arg(X·conj(Xprev)); a zero bin has no product angle: its own angles
        let d = mag[k] && xm[k] ? Math.atan2(im[k] * xr[k] - re[k] * xi[k], re[k] * xr[k] + im[k] * xi[k]) : Math.atan2(im[k], re[k]) - Math.atan2(xi[k], xr[k])
        let dp = wrapPhase(d - k * freqPerBin * anaHop), a = (k * freqPerBin + dp / anaHop) * synHop
        c = Math.cos(a); s = Math.sin(a)
      }
      ur[k] = u * c - v * s
      ui[k] = u * s + v * c
    }
    for (let i = 0, p = -1, rr = 0, ri = 0; i < nl; i++) {
      let k = locked[i], pk = owner[k]
      if (pk !== p) { p = pk; let m = mag[p]; rr = (ur[p] * re[p] + ui[p] * im[p]) / m; ri = (ui[p] * re[p] - ur[p] * im[p]) / m }
      let q = 1 / mag[k]
      ur[k] = (re[k] * rr - im[k] * ri) * q
      ui[k] = (re[k] * ri + im[k] * rr) * q
    }
  }
  xr.set(re.subarray(0, n)); xi.set(im.subarray(0, n)); xm.set(mag.subarray(0, n))
  for (let k = 0; k < n; k++) { re[k] = mag[k] * ur[k]; im[k] = mag[k] * ui[k] }
}
