import test, { almost, ok, is } from 'tst'
import { fft } from 'fourier-transform'
import { centroid, spread, flatness, rolloff, flux, slope, crest, mfcc, ltas, edit, zcr, target, deviation, smooth } from './index.js'
import { zcr as zcrStat } from '@audio/spectral-zcr/audio'

const fs = 44100
const N = 4096

function sine (freq, n, sr = fs) {
	let d = new Float32Array(n)
	for (let i = 0; i < n; i++) d[i] = Math.sin(2 * Math.PI * freq * i / sr)
	return d
}
function saw (freq, n, H = 12, sr = fs) {
	let d = new Float32Array(n)
	for (let i = 0; i < n; i++) {
		let s = 0
		for (let h = 1; h <= H; h++) { if (h * freq < sr / 2) s += Math.sin(2 * Math.PI * h * freq * i / sr) / h }
		d[i] = s
	}
	return d
}
function whiteNoise (n, seed = 7) {
	let d = new Float32Array(n), s = seed
	for (let i = 0; i < n; i++) { s = (s * 1103515245 + 12345) & 0x7fffffff; d[i] = s / 0x3fffffff - 1 }
	return d
}
// hann-windowed magnitude spectrum
function mags (data) {
	let buf = new Float64Array(data.length)
	for (let i = 0; i < data.length; i++) buf[i] = data[i] * (0.5 - 0.5 * Math.cos(2 * Math.PI * i / (data.length - 1)))
	let [re, im] = fft(buf)
	let half = data.length / 2
	let m = new Float32Array(half + 1)
	for (let k = 0; k <= half; k++) m[k] = Math.sqrt(re[k] * re[k] + im[k] * im[k])
	return m
}
// mulberry32 (public domain) — seeded PRNG for pinkNoise below. The file's own `whiteNoise`
// LCG is fine for the loose ≈0.5/flatness checks it's used for elsewhere, but LCGs have a known
// lattice/spectral-test weakness (Knuth, TAOCP Vol.2 §3.3.4) that shows up as multi-dB narrowband
// coloration once filtered — too coarse for the ±1.5 dB LTAS-vs-target check below.
function mulberry32 (seed) {
	let s = seed >>> 0
	return function () {
		s = (s + 0x6d2b79f5) | 0
		let t = Math.imul(s ^ (s >>> 15), 1 | s)
		t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t
		return ((t ^ (t >>> 14)) >>> 0) / 4294967296
	}
}
// Paul Kellet's refined pink noise filter (musicdsp.org "pink noise", public domain), −3 dB/oct
// to within ~±0.05 dB over ~9.5 octaves — colors seeded white noise for the spectral-target ×
// spectral-ltas integration test below.
function pinkNoise (n, seed = 11) {
	let rnd = mulberry32(seed)
	let d = new Float32Array(n)
	let b0 = 0, b1 = 0, b2 = 0, b3 = 0, b4 = 0, b5 = 0, b6 = 0
	for (let i = 0; i < n; i++) {
		let x = rnd() * 2 - 1
		b0 = 0.99886 * b0 + x * 0.0555179
		b1 = 0.99332 * b1 + x * 0.0750759
		b2 = 0.96900 * b2 + x * 0.1538520
		b3 = 0.86650 * b3 + x * 0.3104856
		b4 = 0.55000 * b4 + x * 0.5329522
		b5 = -0.7616 * b5 - x * 0.0168980
		let out = b0 + b1 + b2 + b3 + b4 + b5 + b6 + x * 0.5362
		b6 = x * 0.115926
		d[i] = out * 0.11
	}
	return d
}
// single-bin energy (Goertzel)
function energyAt (data, freq, from = 0, to = data.length, sr = fs) {
	let w = 2 * Math.PI * freq / sr, cw = Math.cos(w)
	let s0 = 0, s1 = 0, s2 = 0
	for (let i = from; i < to; i++) { s0 = data[i] + 2 * cw * s1 - s2; s2 = s1; s1 = s0 }
	let n = to - from
	return Math.sqrt(Math.max(0, s1 * s1 + s2 * s2 - 2 * cw * s1 * s2)) / n
}

test('centroid — pure tone sits at its frequency', () => {
	almost(centroid(mags(sine(1000, N)), { fs, n: N }), 1000, 30)
	almost(centroid(mags(sine(5000, N)), { fs, n: N }), 5000, 60)
	is(centroid(new Float32Array(64)), 0, 'silence → 0')
})

test('spread — tone narrow, noise wide', () => {
	let tone = spread(mags(sine(1000, N)), { fs, n: N })
	let noise = spread(mags(whiteNoise(N)), { fs, n: N })
	ok(tone < 300, `tone spread ${tone.toFixed(0)} Hz`)
	ok(noise > 3000, `noise spread ${noise.toFixed(0)} Hz`)
})

test('flatness — noise ≫ tone, bounded 0..1', () => {
	let fTone = flatness(mags(sine(1000, N)))
	let fNoise = flatness(mags(whiteNoise(N)))
	ok(fTone < 0.02, `tone flatness ${fTone.toFixed(4)}`)
	ok(fNoise > 0.2, `noise flatness ${fNoise.toFixed(3)}`)
	ok(fNoise <= 1 && fTone >= 0)
})

test('rolloff — pure tone rolls off at the tone; noise much higher', () => {
	almost(rolloff(mags(sine(500, N)), { fs, n: N }), 500, 60)
	ok(rolloff(mags(whiteNoise(N)), { fs, n: N }) > 5000)
})

test('flux — zero for identical frames, positive for change; rectified ignores decay', () => {
	let a = mags(sine(1000, N))
	is(flux(a, a), 0)
	ok(flux(mags(sine(2000, N)), a) > 0)
	let quiet = a.map(x => x * 0.5)
	is(flux(quiet, a, { rectified: true }), 0, 'pure decay → no rectified flux')
})

test('slope — sign follows spectral tilt', () => {
	ok(slope(Float32Array.from([4, 3, 2, 1]), { fs: 8, n: 6 }) < 0)
	ok(slope(Float32Array.from([1, 2, 3, 4]), { fs: 8, n: 6 }) > 0)
})

test('crest — flat spectrum → 1, single peak → n', () => {
	almost(crest(Float32Array.from([1, 1, 1, 1])), 1, 1e-9)
	almost(crest(Float32Array.from([0, 0, 1, 0])), 4, 1e-9)
	ok(crest(mags(sine(1000, N))) > crest(mags(whiteNoise(N))), 'tone peakier than noise')
})

test('mfcc — 13 coefficients; c1..12 gain-invariant; timbres differ', () => {
	let x = saw(220, 2048)
	let c = mfcc(x, { fs })
	is(c.length, 13)
	let loud = mfcc(x.map(v => v * 4), { fs })
	for (let k = 1; k < 13; k++) almost(loud[k], c[k], 1e-3, `c${k} gain-invariant`)
	ok(Math.abs(loud[0] - c[0]) > 1, 'c0 tracks level')
	let toneC = mfcc(sine(220, 2048), { fs })
	let dist = 0
	for (let k = 1; k < 13; k++) dist += (toneC[k] - c[k]) ** 2
	ok(Math.sqrt(dist) > 1, 'sine vs saw timbres separate')
})

test('mfcc — tables kept by their parameters: every call the dense definition, to the bit', () => {
	// the definition, built whole on every call: a dense triangle per band over all bins, the DCT's cosines in the sum
	const hz2mel = f => 2595 * Math.log10(1 + f / 700), mel2hz = m => 700 * (10 ** (m / 2595) - 1)
	function dense (data, { fs = 44100, bins = 13, nMel = 40, fmin = 0, fmax = fs / 2 } = {}) {
		let n = data.length, half = n / 2, [re, im] = fft(Float64Array.from(data)), lo = hz2mel(fmin), hi = hz2mel(fmax)
		let c = Array.from({ length: nMel + 2 }, (_, i) => mel2hz(lo + (hi - lo) * i / (nMel + 1)) * n / fs), log = []
		for (let m = 0; m < nMel; m++) {
			let e = 0
			for (let k = 0; k <= half; k++) {
				let w = k < Math.ceil(c[m]) || k >= Math.min(c[m + 2], half + 1) ? 0 : k < c[m + 1] ? (k - c[m]) / (c[m + 1] - c[m]) : (c[m + 2] - k) / (c[m + 2] - c[m + 1])
				e += (re[k] * re[k] + im[k] * im[k]) * (w < 0 ? 0 : w)
			}
			log.push(Math.log(e + 1e-10))
		}
		return Float32Array.from({ length: bins }, (_, k) => log.reduce((s, v, m) => s + v * Math.cos(Math.PI * k * (2 * m + 1) / (2 * nMel)), 0))
	}
	// parameter sets taking turns, so a table kept under too few of them would serve the wrong call
	let sets = [{ fs }, { fs: 16000 }, { fs, bins: 20 }, { fs, nMel: 26 }, { fs, fmin: 300, fmax: 8000 }, { fs: 48000, nMel: 128, bins: 40 }]
	for (let round = 0; round < 2; round++) for (let n of [256, 2048]) for (let o of sets) {
		let x = whiteNoise(n, 3 + n), got = mfcc(x, o), want = dense(x, o)
		is(got.length, want.length)
		ok(got.every((v, i) => Object.is(v, want[i])), `n ${n}, ${JSON.stringify(o)}`)
	}
	ok(mfcc(new Float64Array(64), { fs }).every(v => Number.isFinite(v)), 'silence: log of the floor, finite')
})

test('ltas — stationary tone peaks at its bin', () => {
	let m = ltas(sine(1000, fs * 2), { frameSize: 4096 })
	is(m.length, 2049)
	let peak = 0
	for (let k = 1; k < m.length; k++) if (m[k] > m[peak]) peak = k
	almost(peak * fs / 4096, 1000, fs / 4096 + 1, 'peak bin at 1 kHz')
})

test('edit — no regions reconstructs input', () => {
	let d = sine(440, fs)
	let r = edit(d, { fs })
	let err = 0
	for (let i = 4096; i < d.length - 4096; i++) err = Math.max(err, Math.abs(r[i] - d[i]))
	ok(err < 1e-3, `max reconstruction error ${err.toExponential(1)}`)
})

test('edit — band kill removes one tone, keeps the other', () => {
	let n = fs
	let d = new Float32Array(n)
	let a = sine(500, n), b = sine(3000, n)
	for (let i = 0; i < n; i++) d[i] = a[i] + b[i]
	let r = edit(d, { fs, regions: [{ f0: 2500, f1: 3500, gain: 0 }] })
	let e3k = energyAt(r, 3000, 4096, n - 4096) / energyAt(d, 3000, 4096, n - 4096)
	let e500 = energyAt(r, 500, 4096, n - 4096) / energyAt(d, 500, 4096, n - 4096)
	ok(e3k < 0.02, `3 kHz suppressed to ${(20 * Math.log10(e3k)).toFixed(0)} dB`)
	almost(e500, 1, 0.05, '500 Hz preserved')
})

test('edit — time-scoped region only affects its span', () => {
	let n = fs
	let d = sine(3000, n)
	let r = edit(d, { fs, regions: [{ t0: 0.5, f0: 2500, f1: 3500, gain: 0 }] })
	let head = energyAt(r, 3000, 4096, n / 2 - 4096) / energyAt(d, 3000, 4096, n / 2 - 4096)
	let tail = energyAt(r, 3000, n / 2 + 4096, n - 4096) / energyAt(d, 3000, n / 2 + 4096, n - 4096)
	almost(head, 1, 0.05, 'first half intact')
	ok(tail < 0.05, 'second half suppressed')
})

import { freeze, contrast, harmonics, cqt } from './index.js'

test('freeze — sustains the frozen tone for the whole duration', () => {
	let d = sine(500, fs)
	let out = freeze(d, { fs, at: 0.2, duration: 2 })
	is(out.length, 2 * fs)
	ok(out.every(isFinite))
	let e = (x, f, from, to) => energyAt(x, f, from, to)
	ok(e(out, 500, Math.round(1.7 * fs), 2 * fs - 2048) > 0.05 * e(d, 500, 8192, fs - 8192), 'tone persists at the tail')
})

test('contrast — tone contrasty, noise flat', () => {
	let c1 = contrast(mags(sine(1000, N)), { fs, n: N })
	let c2 = contrast(mags(whiteNoise(N)), { fs, n: N })
	ok(Math.max(...c1) > Math.max(...c2) + 10, 'tone band contrast exceeds noise by >10 dB')
})

test('harmonics — saw signatures: T1 from 1/k amps, low inharmonicity; square odd-dominant', () => {
	let h = harmonics(saw(220, 4096), { fs, f0: 220 })
	// saw amps ∝ 1/k → e_k ∝ 1/k²: T1 = 1/Σ(1/k²) ≈ 0.645 for 10 harmonics
	almost(h.tristimulus[0], 0.645, 0.05, 'T1 ' + h.tristimulus[0].toFixed(3))
	ok(h.inharmonicity < 0.01, 'harmonic tone: inharmonicity ' + h.inharmonicity.toFixed(4))
	let sq = new Float32Array(4096)
	for (let i = 0; i < 4096; i++) { let s = 0; for (let k = 1; k <= 9; k += 2) s += Math.sin(2 * Math.PI * k * 220 * i / fs) / k; sq[i] = s }
	let hq = harmonics(sq, { fs, f0: 220 })
	ok(hq.oddEven > 20, 'square odd/even ' + hq.oddEven.toFixed(1))
})

test('cqt — semitone-spaced bins peak at the played notes', () => {
	let d = new Float32Array(fs)
	for (let i = 0; i < fs; i++) d[i] = 0.5 * Math.sin(2 * Math.PI * 220 * i / fs) + 0.5 * Math.sin(2 * Math.PI * 440 * i / fs)
	let { freqs, mag } = cqt(d, { fs, fmin: 55, at: 0.5 })
	let peak = f => {
		let best = 0
		for (let b = 1; b < mag.length; b++) if (Math.abs(freqs[b] - f) < Math.abs(freqs[best] - f)) best = b
		return mag[best]
	}
	let off = peak(311) // G#4/Eb — between the tones, quiet
	ok(peak(220) > off * 8 && peak(440) > off * 8, 'A3 + A4 bins dominate')
	almost(freqs[12] / freqs[0], 2, 1e-6, 'octave spacing exact')
})

// ── @audio/spectral-pvoc ──
import { findPeaks, nearestPeak, lockPhase, lockMap, lockState, lockAdvance, lockStep, jointState, lockJoint, makeFrameRatio, wrapPhase, PI2, scatterLocked } from '@audio/spectral-pvoc'

test('spectral-pvoc — findPeaks locates isolated spectral peaks', () => {
  let half = 512
  let mag = new Float64Array(half + 1)
  for (let k of [50, 100, 200]) { mag[k] = 1; mag[k - 1] = 0.5; mag[k + 1] = 0.5 }
  let peaks = findPeaks(mag, half)
  is(peaks.length, 3)
  is(peaks[0], 50); is(peaks[1], 100); is(peaks[2], 200)
  is(nearestPeak(peaks, 60), 0, 'bin 60 → peak 50')
  is(nearestPeak(peaks, 90), 1, 'bin 90 → peak 100')
})

test('spectral-pvoc — makeFrameRatio resolves scalar and function ratios', () => {
  let s = makeFrameRatio(1.5)
  is(s.scalar, 1.5); is(s.at(44100, 44100), 1.5)
  let v = makeFrameRatio(t => 1 + t)
  is(v.scalar, 1); is(v.at(44100, 44100), 2, 'ratio(1s) = 2')
  ok(Math.abs(wrapPhase(3 * PI2 + 0.5) - 0.5) < 1e-12, 'wrapPhase')
})

test('spectral-pvoc — lockPhase rigidly co-rotates a peak region', () => {
  let half = 64
  let mag = new Float64Array(half + 1)
  let phase = new Float64Array(half + 1)
  let prop = new Float64Array(half + 1)
  // one strong peak at bin 32 with shoulders
  mag[31] = 0.6; mag[32] = 1; mag[33] = 0.6
  for (let k = 0; k <= half; k++) { phase[k] = 0.1 * k; prop[k] = phase[k] }
  prop[32] = phase[32] + 1.0 // peak advanced by 1 rad
  lockPhase(phase, prop, mag, half)
  ok(Math.abs(prop[31] - (phase[31] + 1.0)) < 1e-12, 'shoulder locked to peak rotation')
  ok(Math.abs(prop[33] - (phase[33] + 1.0)) < 1e-12, 'other shoulder locked')
})

test('spectral-pvoc – lockMap: peaks and quiet bins stay free, the rest ride the nearest peak up to halfway', () => {
  let half = 64, mag = new Float64Array(half + 1), owner = new Int32Array(half + 1)
  // peaks at 20 and 40, sloping to a valley at 30; bin 10 under 3% of peak 20
  for (let k = 20; k <= 30; k++) mag[k] = 1 - 0.08 * (k - 20)
  for (let k = 31; k <= 40; k++) mag[k] = 0.8 - 0.06 * (40 - k)
  mag[19] = 0.5; mag[41] = 0.4; mag[10] = 0.01
  lockMap(mag, half, owner)
  is([owner[20], owner[40]], [-1, -1], 'peaks are free')
  is([owner[19], owner[21], owner[30], owner[31], owner[39], owner[41]], [20, 20, 20, 40, 40, 40], 'bins ride their peak, the region ending halfway between')
  is([owner[10], owner[50]], [-1, -1], 'a bin under 3% of its peak, and a silent one, are free')
})

test('spectral-pvoc – lockAdvance: reset passes the spectrum; magnitudes kept; locked bins keep their offset to the peak', () => {
  let half = 64, w = PI2 / (2 * half), st = lockState(half)
  // two partials near bins 20 and 40 with shoulders, phases advancing at their own frequencies; silence elsewhere
  let frame = t => {
    let re = new Float64Array(half + 1), im = new Float64Array(half + 1)
    for (let [c, a, f] of [[20, 1, 20.3], [40, 0.8, 39.8]]) for (let k = c - 3; k <= c + 3; k++) {
      let m = a * (1 - Math.abs(k - c) / 4), ph = f * w * 256 * t + 0.3 * (k - c)
      re[k] = m * Math.cos(ph); im[k] = m * Math.sin(ph)
    }
    return [re, im]
  }
  let mags = (re, im) => re.map((v, k) => Math.hypot(v, im[k]))
  let arg = (re, im, k, p) => Math.atan2(im[k] * re[p] - re[k] * im[p], re[k] * re[p] + im[k] * im[p])  // arg(X_k / X_p)
  let [re, im] = frame(0), x0 = re.slice(), y0 = im.slice()
  st.mag.set(mags(re, im))
  lockAdvance(re, im, st, true, 256, 205, w, half)
  ok(re.every((v, k) => Math.abs(v - x0[k]) < 1e-12 && Math.abs(im[k] - y0[k]) < 1e-12), 'reset: synthesis = analysis')
  for (let t = 1; t < 5; t++) {
    [re, im] = frame(t)
    let ar = re.slice(), ai = im.slice(), m = mags(re, im)
    st.mag.set(m)
    lockAdvance(re, im, st, false, 256, 205, w, half)
    ok(m.every((v, k) => Math.abs(Math.hypot(re[k], im[k]) - v) < 1e-12), `frame ${t}: every bin keeps its magnitude, silent ones stay silent`)
    ok([17, 19, 21, 23, 37, 39, 41, 43].every(k => { let p = k < 30 ? 20 : 40; return Math.abs(arg(re, im, k, p) - arg(ar, ai, k, p)) < 1e-9 }), `frame ${t}: locked bins keep their analysis offset to the peak`)
  }
})

test('spectral-pvoc – lockAdvance: free bins advance at their instantaneous frequency, by the table and kernels within 1e-12 rad of Math\'s angle', () => {
  // noise-like frames (most quiet bins free), whole hops: the table path; the angle Math gives it, the reference
  let half = 1024, N = 2 * half, w = PI2 / N, A = 511, S = 409, st = lockState(half), seed = 5
  let rnd = () => (seed = (seed * 1664525 + 1013904223) >>> 0) / 2 ** 32 - 0.5
  let frame = () => [Float64Array.from({ length: half + 1 }, rnd), Float64Array.from({ length: half + 1 }, rnd)]
  let worst = 0, free = 0
  for (let t = 0; t < 6; t++) {
    let [re, im] = frame()
    im[0] = im[half] = 0
    let pr = st.xr.slice(), pi = st.xi.slice(), pm = st.xm.slice(), ur = st.ur.slice(), ui = st.ui.slice(), xr = re.slice(), xi = im.slice()
    st.mag.set(re.map((v, k) => Math.hypot(v, im[k])))
    lockAdvance(re, im, st, t === 0, A, S, w, half)
    if (!t) continue
    for (let k = 1; k < half; k++) {
      if (st.owner[k] >= 0 || !st.mag[k] || !pm[k]) continue
      let d = Math.atan2(xi[k] * pr[k] - xr[k] * pi[k], xr[k] * pr[k] + xi[k] * pi[k])
      let dp = wrapPhase(d - k * w * A), a = (k * w + dp / A) * S
      let er = ur[k] * Math.cos(a) - ui[k] * Math.sin(a), ei = ur[k] * Math.sin(a) + ui[k] * Math.cos(a)
      worst = Math.max(worst, Math.hypot(st.ur[k] - er, st.ui[k] - ei)); free++
    }
  }
  ok(free > 1000, `${free} free bins checked`)
  ok(worst < 1e-12, `largest difference ${worst.toExponential(2)}`)
})

test('spectral-pvoc – lockAdvance and lockStep: the WebAssembly kernel gives the JS loop\'s bits', () => {
  // the state lockState makes where WebAssembly does not run: lockAdvance takes the JS loop on it
  let plain = half => { let n = half + 1, f = () => new Float64Array(n); return { mag: f(), owner: new Int32Array(n), xr: f(), xi: f(), xm: f(), ur: f(), ui: f() } }
  let seed = 9, rnd = () => (seed = (seed * 1664525 + 1013904223) >>> 0) / 2 ** 32 - 0.5
  // partials over noise, some bins exactly zero, a silent frame now and then: every path of the step
  let frame = (half, t) => {
    let re = Float64Array.from({ length: half + 1 }, () => t % 7 === 3 ? 0 : 0.01 * rnd()), im = re.map(() => t % 7 === 3 ? 0 : 0.01 * rnd())
    for (let [c, a] of [[half / 8 + 0.3, 1], [half / 3 + 0.7, 0.5]]) for (let k = Math.max(0, Math.floor(c) - 3); k <= Math.min(half, Math.floor(c) + 3); k++) {
      let ph = c * 0.37 * t + 0.2 * k
      re[k] += a / (1 + Math.abs(k - c)) * Math.cos(ph); im[k] += a / (1 + Math.abs(k - c)) * Math.sin(ph)
    }
    for (let k = 0; k <= half; k += 13) re[k] = im[k] = 0
    im[0] = im[half] = 0
    return [re, im]
  }
  let same = (a, b) => a.length === b.length && a.every((v, i) => Object.is(v, b[i]))
  let kernel = 0, bad = 0
  for (let half of [2, 4, 64, 1024]) for (let step of [false, true]) {
    let N = 2 * half, w = PI2 / N, a = plain(half), b = lockState(half)
    if ('k' in b) kernel++
    for (let t = 0; t < 40; t++) {
      let A = Math.max(1, Math.round(N / 4 * (0.7 + 0.3 * Math.sin(t)))), S = Math.max(1, Math.floor(N / 5)), reset = t === 0 || t === 20
      let [re, im] = frame(half, t), re2 = re.slice(), im2 = im.slice()
      a.mag.set(re.map((v, k) => Math.sqrt(v * v + im[k] * im[k])))
      lockAdvance(re, im, a, reset, A, S, w, half)
      if (step) lockStep(re2, im2, b, reset, A, S, w, half)
      else { b.mag.set(re2.map((v, k) => Math.sqrt(v * v + im2[k] * im2[k]))); lockAdvance(re2, im2, b, reset, A, S, w, half) }
      if (!same(re, re2) || !same(im, im2) || !['mag', 'xr', 'xi', 'xm', 'ur', 'ui'].every(k => same(a[k], b[k])) || !reset && !same(a.owner, b.owner)) bad++
    }
  }
  ok(typeof WebAssembly !== 'object' || kernel === 8, 'the states run on the kernel where WebAssembly runs')
  is(bad, 0, 'every frame the same bits: synthesis bins, magnitudes, the state and the lock map')
})

// ── @audio/spectral-zcr ──
function frameZcr (d, N = 2048, HOP = 512) {
	let acc = 0, cnt = 0
	for (let off = 0; off + N <= d.length; off += HOP) { acc += zcr(d.subarray(off, off + N)); cnt++ }
	return cnt ? acc / cnt : 0
}

test('zcr — sine crossings ≈ 2f/sr and scale linearly with frequency', () => {
	// librosa.feature.zero_crossing_rate: a sine crosses zero twice per period → rate = 2f/sr
	almost(frameZcr(sine(440, fs)), 2 * 440 / fs, 2 * 440 / fs * 0.03, '440 Hz ≈ 0.01995')
	almost(frameZcr(sine(4000, fs)), 2 * 4000 / fs, 2 * 4000 / fs * 0.03, '4 kHz ≈ 0.1814, scales linearly')
})

test('zcr — white noise ≈ 0.5; DC-offset constant signal → 0', () => {
	// independent symmetric samples flip sign with p ≈ 0.5 (librosa.feature.zero_crossing_rate docs, noise example)
	almost(frameZcr(whiteNoise(fs)), 0.5, 0.05)
	ok(frameZcr(new Float32Array(fs).fill(0.5)) === 0, 'constant signal never crosses zero')
})

test('zcr — audio.js manifest matches the mono-kernel average on a stereo pair', () => {
	let n = fs
	let l = sine(440, n), r = sine(440, n).map(v => v * 0.3)
	let mono = new Float32Array(n)
	for (let i = 0; i < n; i++) mono[i] = (l[i] + r[i]) / 2
	almost(zcrStat.compute([l, r], { sampleRate: fs }), frameZcr(mono), 1e-9)
})

test('zcr — np.signbit semantics (−0 counts negative) and empty-input guard', () => {
	// librosa zero_crossings uses np.signbit: [0.1, −0, 0.2] → signs [+,−,+] → 2 changes / 3 samples
	almost(zcr(Float32Array.from([0.1, -0, 0.2])), 2 / 3, 1e-9)
	is(zcr(new Float32Array(0)), 0)
})

// ── @audio/spectral-target ──

test('target — speech anchors render exactly (Byrne, Dillon, Tran et al. 1994, JASA 96(4):2108, Table II "Combined" col, p.2116, normalized to 70 dB SPL overall: 500→62.1, 1000→53.7, 2000→48.7, 4000→45.6 dB)', () => {
	let fs = 44100, bins = 442 // n=882 → 50 Hz/bin, so 500/1000/2000/4000 Hz land on exact bins
	let curve = target('speech', { fs, bins })
	let n = 2 * (bins - 1)
	let bin = f => Math.round(f * n / fs)
	let hz = [500, 1000, 2000, 4000]
	let published = [62.1, 53.7, 48.7, 45.6]
	let rendered = hz.map(f => curve[bin(f)])
	let rMean = rendered.reduce((a, b) => a + b, 0) / rendered.length
	let pMean = published.reduce((a, b) => a + b, 0) / published.length
	for (let i = 0; i < hz.length; i++)
		almost(rendered[i] - rMean, published[i] - pMean, 0.05, `${hz[i]} Hz anchor exact (mean-normalization-invariant check)`)
})

test('target — pink is exactly −3.0103 dB/octave (analytic, equal energy per octave)', () => {
	let fs = 44100, bins = 2049
	let curve = target('pink', { fs, bins })
	for (let k of [25, 100, 400]) // doubling the bin index doubles the frequency exactly
		almost(curve[2 * k] - curve[k], -3.0103, 0.01, `bin ${k}→${2 * k}`)
})

test('target — custom anchor table renders an exact straight log-frequency line', () => {
	let fs = 44100, bins = 883 // n=1764 → 25 Hz/bin, so 100/200/400/800/1600 Hz land on exact bins
	let curve = target([[100, 0], [1600, -40]], { fs, bins }) // 4 octaves, −10 dB/octave
	let n = 2 * (bins - 1)
	let bin = f => Math.round(f * n / fs)
	let v100 = curve[bin(100)]
	almost(curve[bin(200)] - v100, -10, 1e-4, '1 octave in')
	almost(curve[bin(400)] - v100, -20, 1e-4, '2 octaves in')
	almost(curve[bin(800)] - v100, -30, 1e-4, '3 octaves in')
	almost(curve[bin(1600)] - v100, -40, 1e-4, '4 octaves in — the far anchor')
})

test('deviation — target against itself is ~0 inside the band, exactly 0 at/outside the band edges', () => {
	let fs = 44100, bins = 2049
	let n = 2 * (bins - 1)
	let td = target('speech', { fs, bins })
	let measured = Float32Array.from(td, db => 10 ** (db / 20)) // dB → linear mag
	let dev = deviation(measured, td, { fs })
	let maxAbs = 0
	for (let k = 0; k < bins; k++) {
		let f = k * fs / n
		if (f > 50 && f < 19845) maxAbs = Math.max(maxAbs, Math.abs(dev[k]))
	}
	ok(maxAbs <= 0.05, `identity deviation ${maxAbs.toFixed(4)} dB inside the band`)
	is(dev[0], 0, 'DC — outside [20 Hz, 0.45·fs], exactly 0')
	is(dev[bins - 1], 0, 'Nyquist — outside [20 Hz, 0.45·fs], exactly 0')
})

test('deviation — flat ("white noise") measured vs pink target slopes correction downward with frequency', () => {
	let fs = 44100, bins = 2049
	let n = 2 * (bins - 1)
	let flat = new Float32Array(bins).fill(1) // flat linear magnitude = flat LTAS
	let td = target('pink', { fs, bins })
	let dev = deviation(flat, td, { fs })
	let bin = f => Math.round(f * n / fs)
	let c500 = dev[bin(500)], c8k = dev[bin(8000)]
	ok(c8k - c500 <= -6, `correction(8kHz)−correction(500Hz) = ${(c8k - c500).toFixed(2)} dB (expect ≤ −6, HF cut)`)
})

test('deviation — clamp bounds the correction everywhere, even under extreme mismatch', () => {
	let fs = 44100, bins = 2049
	let flat = new Float32Array(bins).fill(1)
	let td = target('speech', { fs, bins }) // steep HF falloff vs a flat "measured" → large mismatch
	let clamp = 8
	let dev = deviation(flat, td, { fs, clamp })
	let maxAbs = 0
	for (let k = 0; k < bins; k++) maxAbs = Math.max(maxAbs, Math.abs(dev[k]))
	ok(maxAbs <= clamp + 1e-6, `|correction| ≤ ${clamp} dB everywhere (max ${maxAbs.toFixed(2)})`)
})

// 1.0 levelled both curves by the mean of their bins over [20 Hz, 0.45·fs]: the octave under fs/2 holds half the bins,
// and a preset is extrapolated past its anchors (speech held flat past 10 kHz), so a voice whose spectrum is the target's
// under 10 kHz and falls away above read ~12 dB over the target through the mids, and was given treble to 16 kHz.
test('deviation: a voice band-limited at 10 kHz, the speech target\'s shape under it, reads no correction (the target\'s band)', () => {
	let fs = 48000, bins = 2049, n = 2 * (bins - 1)
	let td = target('speech', { fs, bins })
	let m = Float32Array.from(td, (db, k) => 10 ** ((db - (k * fs / n > 10000 ? 40 : 0)) / 20))
	let dev = deviation(m, td, { fs }), maxAbs = 0
	for (let k = 1; k < bins; k++) maxAbs = Math.max(maxAbs, Math.abs(dev[k]))
	ok(maxAbs <= 0.05, `max |correction| ${maxAbs.toFixed(2)} dB`)
	is(String(td.band), '100,10000', 'the target carries its anchors\' band')
})

test('deviation: levelled per octave: a curve 20 dB down under 1 kHz is met halfway, +10 under and −10 over', () => {
	let fs = 48000, bins = 2049, n = 2 * (bins - 1), bin = f => Math.round(f * n / fs)
	let td = target([[100, 0], [10000, 0]], { fs, bins })     // flat over 100 Hz–10 kHz: 1 decade either side of 1 kHz
	let m = Float32Array.from(td, (db, k) => 10 ** ((db - (k * fs / n < 1000 ? 20 : 0)) / 20))
	let dev = deviation(m, td, { fs, smoothOct: 1 / 3 })
	almost(dev[bin(300)], 10, 0.5, `300 Hz ${dev[bin(300)].toFixed(2)} dB`)
	almost(dev[bin(3000)], -10, 0.5, `3 kHz ${dev[bin(3000)].toFixed(2)} dB`)
})

test('smooth — a single-bin spike is reduced ≥6 dB by 1/3-oct smoothing, total mean preserved', () => {
	let fs = 44100, bins = 2049
	let db = new Float32Array(bins)
	db[500] = 12
	let sm = smooth(db, { fs, oct: 1 / 3 })
	ok(sm[500] <= 6, `spike 12 → ${sm[500].toFixed(2)} dB`)
	let meanBefore = db.reduce((a, b) => a + b, 0) / bins
	let meanAfter = sm.reduce((a, b) => a + b, 0) / bins
	almost(meanAfter, meanBefore, 0.1, `mean preserved: ${meanBefore.toFixed(5)} → ${meanAfter.toFixed(5)}`)
})

test('spectral-target × spectral-ltas — pink noise LTAS deviates ~0 dB from the pink target, 100 Hz–10 kHz', () => {
	let sig = pinkNoise(20 * fs, 11) // seeded, 20 s
	let m = ltas(sig, { frameSize: 4096, hop: 2048 })
	let td = target('pink', { fs, bins: m.length })
	let dev = deviation(m, td, { fs })
	let n = 2 * (m.length - 1)
	let maxAbs = 0
	for (let k = 0; k < m.length; k++) {
		let f = k * fs / n
		if (f >= 100 && f <= 10000) maxAbs = Math.max(maxAbs, Math.abs(dev[k]))
	}
	ok(maxAbs <= 1.5, `pink-noise LTAS vs pink target within ${maxAbs.toFixed(2)} dB, 100 Hz–10 kHz`)
})


// Channels of one sound: partials the channels share, each at its own level and phase in every channel, over noise of
// each channel's own; with `gaps`, a silent frame now and then and some bins exactly zero
const jointFrames = (half, nch, seed = 3, gaps = true) => {
  let rnd = () => (seed = (seed * 1664525 + 1013904223) >>> 0) / 2 ** 32 - 0.5
  let lv = Array.from({ length: nch }, (_, c) => [1 - 0.3 * c, 0.4 + 0.2 * c]), off = Array.from({ length: nch }, (_, c) => [0.9 * c, -1.7 * c])
  return t => {
    let res = [], ims = []
    for (let c = 0; c < nch; c++) {
      let quiet = gaps && t % 7 === 3
      let re = Float64Array.from({ length: half + 1 }, () => quiet ? 0 : 0.01 * rnd()), im = re.map(() => quiet ? 0 : 0.01 * rnd())
      if (!quiet) [[half / 8 + 0.3, 0], [half / 3 + 0.7, 1]].forEach(([f, i]) => {
        for (let k = Math.max(0, Math.floor(f) - 3); k <= Math.min(half, Math.floor(f) + 3); k++) {
          let ph = f * 0.37 * t + 0.2 * k + off[c][i], a = lv[c][i] / (1 + Math.abs(k - f))
          re[k] += a * Math.cos(ph); im[k] += a * Math.sin(ph)
        }
      })
      if (gaps) for (let k = 0; k <= half; k += 13) re[k] = im[k] = 0
      im[0] = im[half] = 0
      res.push(re); ims.push(im)
    }
    return [res, ims]
  }
}

test('spectral-pvoc – lockJoint: every channel keeps its magnitudes, and its level and phase against the others, in every frame', () => {
  let half = 256, N = 2 * half, w = PI2 / N, nch = 3, j = jointState(half, nch), frame = jointFrames(half, nch), worstMag = 0, worstRatio = 0, turned = 0
  for (let t = 0; t < 30; t++) {
    let [res, ims] = frame(t), xr = res.map(a => a.slice()), xi = ims.map(a => a.slice())
    lockJoint(res, ims, j, t === 0, 128, 103, w, half)
    for (let k = 0; k <= half; k++) {
      for (let c = 0; c < nch; c++) worstMag = Math.max(worstMag, Math.abs(Math.hypot(res[c][k], ims[c][k]) - Math.hypot(xr[c][k], xi[c][k])))
      // Y_c·conj(Y_0) against X_c·conj(X_0): the cross spectrum of two channels, their level and phase apart
      for (let c = 1; c < nch; c++) {
        let yr = res[c][k] * res[0][k] + ims[c][k] * ims[0][k], yi = ims[c][k] * res[0][k] - res[c][k] * ims[0][k]
        let ar = xr[c][k] * xr[0][k] + xi[c][k] * xi[0][k], ai = xi[c][k] * xr[0][k] - xr[c][k] * xi[0][k]
        worstRatio = Math.max(worstRatio, Math.hypot(yr - ar, yi - ai))
      }
      if (t && Math.hypot(res[0][k] - xr[0][k], ims[0][k] - xi[0][k]) > 1e-3) turned++
    }
  }
  ok(worstMag < 1e-12, `magnitudes kept: ${worstMag.toExponential(2)}`)
  ok(worstRatio < 1e-12, `the cross spectrum between channels kept: ${worstRatio.toExponential(2)}`)
  ok(turned > 100, `the bins do turn (${turned} of them): the hops differ`)
})

test('spectral-pvoc – lockJoint: one channel is lockStep\'s frame; channels of one signal at different levels take its rotation', () => {
  // no silence in them: a bin that wakes from nothing has no advance to measure, and each step gives it an angle of its own
  let half = 256, N = 2 * half, w = PI2 / N, frame = jointFrames(half, 1, 3, false), one = jointState(half, 1), st = lockState(half), two = jointState(half, 2), worst = 0, worst2 = 0
  for (let t = 0; t < 30; t++) {
    let [[re], [im]] = frame(t), a = re.slice(), b = im.slice(), l = re.slice(), li = im.slice(), r = re.map(v => 0.25 * v), ri = im.map(v => 0.25 * v)
    lockJoint([re], [im], one, t === 0, 128, 103, w, half)
    lockStep(a, b, st, t === 0, 128, 103, w, half)
    lockJoint([l, r], [li, ri], two, t === 0, 128, 103, w, half)
    for (let k = 0; k <= half; k++) {
      worst = Math.max(worst, Math.hypot(re[k] - a[k], im[k] - b[k]))
      worst2 = Math.max(worst2, Math.hypot(l[k] - a[k], li[k] - b[k]), Math.hypot(4 * r[k] - a[k], 4 * ri[k] - b[k]))
    }
  }
  ok(worst < 1e-9, `one channel: within ${worst.toExponential(2)} of lockStep`)
  ok(worst2 < 1e-9, `a signal and a quarter of it: each within ${worst2.toExponential(2)} of lockStep on it`)
})

test('spectral-pvoc – lockJoint: the WebAssembly kernel gives the JS loop\'s bits', () => {
  // the state jointState makes where WebAssembly does not run: lockJoint takes the JS loops on it
  let plain = (half, nch) => { let n = half + 1, f = () => new Float64Array(n), each = () => Array.from({ length: nch }, f); return { st: { mag: f(), owner: new Int32Array(n), xr: f(), xi: f(), xm: f(), ur: f(), ui: f() }, re: f(), im: f(), fr: f().fill(1), fi: f(), xr: each(), xi: each() } }
  let same = (a, b) => a.length === b.length && a.every((v, i) => Object.is(v, b[i]))
  let kernel = 0, bad = 0, frames = 0
  for (let half of [2, 4, 64, 1024]) for (let nch of [1, 2, 3, 5]) {
    let N = 2 * half, w = PI2 / N, a = plain(half, nch), b = jointState(half, nch), frame = jointFrames(half, nch, 11 + nch)
    if (b.st.k?.ref) kernel++
    for (let t = 0; t < 40; t++, frames++) {
      // whole hops (the kernel's), a fractional one now and then (the loops on the kernel's memory), a reset mid-way
      let A = Math.max(1, Math.round(N / 4 * (0.7 + 0.3 * Math.sin(t)))) + (t % 9 === 5 ? 0.5 : 0), S = Math.max(1, Math.floor(N / 5)), reset = t === 0 || t === 20
      let [res, ims] = frame(t), res2 = res.map(x => x.slice()), ims2 = ims.map(x => x.slice())
      lockJoint(res, ims, a, reset, A, S, w, half)
      lockJoint(res2, ims2, b, reset, A, S, w, half)
      let state = ['re', 'im', 'fr', 'fi'].every(k => same(a[k], b[k])) && ['mag', 'xr', 'xi', 'xm', 'ur', 'ui'].every(k => same(a.st[k], b.st[k])) && a.xr.every((x, c) => same(x, b.xr[c]) && same(a.xi[c], b.xi[c]))
      if (!res.every((x, c) => same(x, res2[c]) && same(ims[c], ims2[c])) || !state) bad++
    }
  }
  ok(typeof WebAssembly !== 'object' || kernel === 16, 'the states run on the kernel where WebAssembly runs')
  is(bad, 0, `${frames} frames the same bits: every channel's synthesis bins, the reference and the state`)
})

test('spectral-pvoc — locked scatter tracks source phases across silence and bin moves', () => {
  const ctx = { N: 32, half: 16, hop: 8, freqPerBin: PI2 / 32 }
  const create = () => ({ syn: new Float64Array(17), prev: null })
  const frame = (state, peak, ratio, n, reset = false) => {
    const mag = new Float64Array(17), phase = new Float64Array(17)
    if (peak) { mag[peak] = 1; mag[peak - 1] = mag[peak + 1] = .5 }
    for (let k = 0; k <= 16; k++) phase[k] = wrapPhase(.13 * k + n * k * ctx.freqPerBin * ctx.hop)
    const out = new Float64Array(17), angles = new Float64Array(17)
    scatterLocked(mag, phase, state.prev, reset || !state.prev, peak ? [peak] : [], ratio,
      ctx, state.syn, out, angles, new Int32Array(1), new Float64Array(1), new Float64Array(17))
    state.prev = phase
    return { mag, phase, out, angles }
  }
  const a = create()
  is([...frame(a, 0, 1, 0).out].every(x => x === 0), true, 'silent frame')
  for (const [n, peak] of [[1, 3], [2, 4], [3, 3]]) {
    const f = frame(a, peak, 1, n)
    for (let k = peak - 1; k <= peak + 1; k++) {
      almost(f.out[k], f.mag[k], 1e-12, 'unity magnitude')
      almost(wrapPhase(f.angles[k] - f.phase[k]), 0, 1e-12, 'unity phase after peak moves')
    }
  }
  const b = create(), c = create()
  const first = frame(b, 3, 1, 0)
  frame(c, 8, .8, 0) // Interleaving must not share phase history.
  const moved = frame(b, 3, 1.7, 1), dest = 5
  const expected = first.angles[3] + 3 * Math.PI + 3 * ctx.freqPerBin * 1.7 * ctx.hop - dest * Math.PI
  almost(wrapPhase(moved.angles[dest] - expected), 0, 1e-12, 'destination move preserves center phase')
  frame(b, 0, 1, 2)
  const resumed = frame(b, 6, 1, 3)
  almost(wrapPhase(resumed.angles[6] - resumed.phase[6]), 0, 1e-12, 'new peak after silence initializes from analysis')
  const reset = frame(b, 6, 1, 4, true)
  almost(wrapPhase(reset.angles[6] - reset.phase[6]), 0, 1e-12, 'explicit reset')
})


test('spectral-pvoc — locked scatter: every bin rides its nearest peak (nearestPeak, a tie to the lower), moved by its shift', () => {
  // bins of their own energy, peak sets with ties (a bin halfway between two peaks), one peak, and found ones: the bins
  // a scatter fills are exactly those the nearest peak's shift sends a bin to
  let half = 256, ctx = { N: 512, half, hop: 128, freqPerBin: PI2 / 512 }, seed = 11, bad = 0
  const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647
  for (let t = 0; t < 60; t++) {
    let mag = Float64Array.from({ length: half + 1 }, () => 0.01 + rnd()), phase = Float64Array.from(mag, () => PI2 * rnd() - Math.PI)
    let peaks = t % 3 === 0 ? [10, 20, 40, 60, 61, 200] : t % 3 === 1 ? [128] : findPeaks(mag, half), ratio = 0.6 + t / 50
    let out = new Float64Array(half + 1), dest = new Int32Array(peaks.length)
    scatterLocked(mag, phase, null, true, peaks, ratio, ctx, new Float64Array(half + 1), out, new Float64Array(half + 1), dest, new Float64Array(peaks.length), new Float64Array(half + 1))
    let want = new Set()
    for (let k = 0; k <= half; k++) {
      let pi = nearestPeak(peaks, k), d = dest[pi] + k - peaks[pi]
      if (dest[pi] >= 0 && d >= 0 && d <= half) want.add(d)
    }
    for (let k = 0; k <= half; k++) if ((out[k] > 0) !== want.has(k)) bad++
  }
  is(bad, 0, 'the bins filled, as nearestPeak assigns them')
})

test('spectral-pvoc — unresolved boundary frames preserve energy, including DC-only input', () => {
  for (const half of [0, 1, 16]) {
    const mag = Float64Array.from({ length: half + 1 }, (_, k) => 1 / (k + 1))
    const phase = Float64Array.from(mag, (_, k) => .2 * k)
    const out = new Float64Array(half + 1), angles = out.slice(), syn = out.slice()
    scatterLocked(mag, phase, null, true, [], 1.5,
      { N: Math.max(2, half * 2), half, hop: 1, freqPerBin: Math.PI / Math.max(1, half) },
      syn, out, angles, new Int32Array(0), new Float64Array(0), out.slice())
    ok(out.every((x, k) => x === mag[k]), 'no resolved pitch: preserve magnitude')
    ok(angles.every((x, k) => x === phase[k]), 'preserve phase')
    ok(syn.every(Number.isNaN), 'next resolved peak starts fresh')
  }
})

// frames for the forms a real-time caller uses: partials with their lobes over a floor, a plateau now and then, silence
function pvocFrames (half, n, seed) {
  const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647
  return Array.from({ length: n }, (_, t) => {
    let mag = Float64Array.from({ length: half + 1 }, () => 1e-5 * rnd()), phase = Float64Array.from(mag, () => PI2 * rnd() - Math.PI)
    if (t % 9 === 4) mag.fill(0)
    else for (let p = 0; p < Math.min(24, half >> 2); p++) {
      let c = 1 + rnd() * (half - 2), a = rnd()
      for (let k = Math.max(0, Math.floor(c) - 2); k <= Math.min(half, Math.ceil(c) + 2); k++) mag[k] += a * Math.max(0, 1 - Math.abs(k - c) / 2.5)
    }
    if (t % 7 === 3) mag.fill(0.7, 5, 7)
    return { mag, phase }
  })
}

test('spectral-pvoc – findPeaks(mag, half, out): the peaks of the array form, in out, their count returned; at most half / 2; an out too short is a RangeError', () => {
  let bad = 0, frames = 0
  for (const half of [2, 16, 256, 1024]) for (const { mag } of pvocFrames(half, 30, 5 + half)) {
    let peaks = findPeaks(mag, half), out = new Int32Array(half >> 1).fill(-1), n = findPeaks(mag, half, out)
    if (n !== peaks.length || peaks.some((k, i) => out[i] !== k)) bad++
    frames++
  }
  is(bad, 0, `${frames} frames: the count and the peaks`)
  // every odd bin a peak: half / 2 of them, out filled exactly
  let half = 64, mag = Float64Array.from({ length: half + 1 }, (_, k) => k & 1), out = new Int32Array(half >> 1)
  is(findPeaks(mag, half, out), half >> 1, 'every other bin: half / 2 peaks')
  is(out[out.length - 1], half - 1, 'the last at half − 1')
  let threw = null
  try { findPeaks(mag, half, new Int32Array((half >> 1) - 1)) } catch (e) { threw = e }
  ok(threw instanceof RangeError, 'an out shorter than half / 2: RangeError')
})

test('spectral-pvoc – scatterLocked with findPeaks\' out and its count, scatterGated and scatterLocked with the ratio in a Float64Array: the bits of the array and the number', async () => {
  const { scatterGated } = await import('@audio/spectral-pvoc')
  let bad = 0, frames = 0
  for (const half of [16, 256, 1024]) {
    let N = 2 * half, ctx = { N, half, hop: N >> 2, freqPerBin: PI2 / N }, z = () => new Float64Array(half + 1)
    // two scatterLocked states, one a form, carried frame to frame; a reset now and then
    let a = { syn: z() }, b = { syn: z() }, prev = null, list = pvocFrames(half, 40, 11 + half), cell = new Float64Array(1), out = new Int32Array(half >> 1)
    list.forEach(({ mag, phase }, t) => {
      let ratio = [1, 0.5, 2, 4, 0.71, 1.2599210498948732][t % 6] * (1 + t / 1000), reset = t === 0 || t % 13 === 7
      cell[0] = ratio
      let peaks = findPeaks(mag, half), n = findPeaks(mag, half, out)
      let run = (p, r, st, count) => { let m = z(), q = z(), pk = z(), dest = new Int32Array(half), sp = new Float64Array(half); scatterLocked(mag, phase, reset ? null : prev, reset, p, r, ctx, st.syn, m, q, dest, sp, pk, count); return [m, q, pk, st.syn, dest.subarray(0, peaks.length)] }
      let x = run(peaks, ratio, a), y = run(out, cell, b, n)
      if (!x.every((v, i) => v.length === y[i].length && v.every((w, k) => Object.is(w, y[i][k])))) bad++
      let gated = r => { let m = z(), q = z(), pk = z(); scatterGated(mag, phase, prev, r, ctx, m, q, pk); return [m, q, pk] }
      let g = gated(ratio), h = gated(cell)
      if (!g.every((v, i) => v.every((w, k) => Object.is(w, h[i][k])))) bad++
      prev = phase
      frames++
    })
  }
  is(bad, 0, `${frames} frames, state carried: the same bits`)
})

// The validator's measure (@audio/compile tools/verify.js, `young`): in a child Node with its flags, a case a child, the
// bytes V8's young generation takes in over 2^14 warm calls, read in two forced minor collections, the lesser of two runs.
// The ratio is computed each call, as a frame computes it: a number made once would not show a number allocated a call
const PVOC_ALLOC_CHILD = `
import v8 from 'node:v8'
const { findPeaks, scatterGated, scatterLocked } = await import(${JSON.stringify(new URL('./packages/spectral-pvoc/pvoc.js', import.meta.url).href)})
const minor = { type: 'minor' }, young = s => s.reduce((a, x) => x.spaceName === 'new_space' || x.spaceName === 'new_large_object_space' ? a + x.spaceUsedSize : a, 0)
const count = f => { const p = new v8.GCProfiler(); p.start(); gc(minor); f(); gc(minor); const g = p.stop().statistics; let y = 0; for (let i = 1; i < g.length; i++) y += young(g[i].beforeGC.heapSpaceStatistics) - young(g[i - 1].afterGC.heapSpaceStatistics); return y }
const half = 1024, N = 2048, ctx = { N, half, hop: 512, freqPerBin: 2 * Math.PI / N }, CALLS = 1 << 14
let s = 3; const rnd = () => (s = (s * 16807) % 2147483647) / 2147483647
const frames = Array.from({ length: 8 }, () => {
  let mag = Float64Array.from({ length: half + 1 }, () => 1e-4 * rnd()), phase = Float64Array.from(mag, () => 2 * Math.PI * rnd() - Math.PI)
  for (let p = 0; p < 24; p++) { let c = 2 + rnd() * (half - 4); for (let k = Math.floor(c) - 2; k <= Math.ceil(c) + 2; k++) mag[k] += rnd() * Math.max(0, 1 - Math.abs(k - c) / 2.5) }
  return { mag, phase }
})
const z = () => new Float64Array(half + 1), m = z(), q = z(), pk = z(), syn = z(), dest = new Int32Array(half), sp = new Float64Array(half), cell = new Float64Array(1)
const peaks = frames.map(f => { let o = new Int32Array(half >> 1), n = findPeaks(f.mag, half, o); return [o, n] })
const ratio = i => 1.2599210498948732 + (i & 7) * 1e-3
const per = call => { const run = n => { for (let i = 0; i < n; i++) call(i, frames[i & 7], frames[(i + 7) & 7].phase) }; run(CALLS); run(CALLS); return Math.min(count(() => run(CALLS)), count(() => run(CALLS))) / CALLS }
const cases = {
  'findPeaks into out': () => { let o = new Int32Array(half >> 1); return per((i, f) => findPeaks(f.mag, half, o)) },
  'scatterGated, the ratio in a Float64Array': () => per((i, f, prev) => { cell[0] = ratio(i); m.fill(0); q.fill(0); pk.fill(0); scatterGated(f.mag, f.phase, prev, cell, ctx, m, q, pk) }),
  'scatterLocked, peaks in out with their count, the ratio in a Float64Array': () => per((i, f, prev) => { cell[0] = ratio(i); m.fill(0); q.fill(0); pk.fill(0); scatterLocked(f.mag, f.phase, prev, false, peaks[i & 7][0], cell, ctx, syn, m, q, dest, sp, pk, peaks[i & 7][1]) }),
  'control, scatterGated, the ratio a number': () => per((i, f, prev) => { m.fill(0); q.fill(0); pk.fill(0); scatterGated(f.mag, f.phase, prev, ratio(i), ctx, m, q, pk) }),
}
const name = process.argv[1]
console.log(JSON.stringify(name ? cases[name]() : Object.keys(cases)))
`
test('spectral-pvoc – a real-time frame\'s calls allocate nothing: findPeaks into out, scatterGated and scatterLocked with the ratio in a Float64Array', async () => {
  const { default: v8 } = await import('node:v8')
  if (!v8.GCProfiler) return console.warn(`skipped: Node ${process.version} has no v8.GCProfiler to count with (18.15 and later have it)`)
  const { execFileSync } = await import('node:child_process')
  const child = (...a) => JSON.parse(execFileSync(process.execPath, ['--expose-gc', '--no-concurrent-recompilation', '--no-allocation-site-pretenuring', '--input-type=module', '-e', PVOC_ALLOC_CHILD, ...a], { encoding: 'utf8' }))
  for (const k of child()) {
    const b = child(k)
    // a ratio that is not whole, passed as a number to a function V8 does not inline, is a 16-byte number a call
    if (k.startsWith('control')) ok(b >= 8, `the measure sees a number a call: ${k}, ${b} bytes a call`)
    else ok(b <= 1 / 16, `${k}: ${b} bytes a call (at most 1/16)`)
  }
})

// The family's manifests: audio.js names every member that has a manifest (a package of packages/ the umbrella depends
// on, with an audio.js), and no other; each name a member's manifest exports is the entry's, the same binding
test('audio.js lists every member that has a manifest', async () => {
	const { readdirSync, existsSync, readFileSync } = await import('node:fs')
	const at = p => new URL(p, import.meta.url), { dependencies } = JSON.parse(readFileSync(at('package.json'), 'utf8'))
	const members = readdirSync(at('packages')).filter(d => existsSync(at(`packages/${d}/audio.js`)))
		.map(d => JSON.parse(readFileSync(at(`packages/${d}/package.json`), 'utf8')).name).filter(n => n in dependencies).sort()
	const listed = [...readFileSync(at('audio.js'), 'utf8').matchAll(/^export \* from '(@audio\/[^']+)\/audio'$/gm)].map(m => m[1]).sort()
	is(listed, members, 'one line per member, no other')
	const family = await import(at('audio.js')), missing = []
	for (const m of members) for (const [k, v] of Object.entries(await import(`${m}/audio`))) if (k !== 'default' && family[k] !== v) missing.push(`${m}: ${k}`)
	is(missing, [], 'every name a member\'s manifest exports, the same binding')
})
