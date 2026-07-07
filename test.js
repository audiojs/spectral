import test, { almost, ok, is } from 'tst'
import { fft } from 'fourier-transform'
import { centroid, spread, flatness, rolloff, flux, slope, crest, mfcc, ltas, edit } from './index.js'

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
