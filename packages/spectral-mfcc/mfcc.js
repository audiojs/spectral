// MFCC — mel-frequency cepstral coefficients: FFT → mel filterbank (triangular,
// Davis & Mermelstein 1980) → log → DCT-II. Ported from audio core.

import { fft } from 'fourier-transform'

const hz2mel = f => 2595 * Math.log10(1 + f / 700)
const mel2hz = m => 700 * (10 ** (m / 2595) - 1)

// The filterbank and the DCT depend on the options and the frame's length alone: each built once and kept by its
// parameters, as an FFT keeps its twiddles by its size (the last 32 kept)
const cache = new Map()
function memo (key, make) {
	let v = cache.get(key)
	if (v) return v
	if (cache.size >= 32) cache.delete(cache.keys().next().value)
	cache.set(key, v = make())
	return v
}

// triangular mel filterbank weights over FFT power bins, each band as the bins it covers: its first bin and their
// weights (the bins outside weigh 0)
function melBank (nMel, nBins, fs, n, fmin, fmax) {
	let melMin = hz2mel(fmin), melMax = hz2mel(fmax)
	let centers = new Float64Array(nMel + 2)
	for (let i = 0; i < nMel + 2; i++) centers[i] = mel2hz(melMin + (melMax - melMin) * i / (nMel + 1)) * n / fs
	let bank = []
	for (let m = 0; m < nMel; m++) {
		let lo = centers[m], mid = centers[m + 1], hi = centers[m + 2]
		let k0 = Math.max(0, Math.ceil(lo)), k1 = Math.min(hi, nBins)
		let w = new Float64Array(Math.max(0, Math.ceil(k1) - k0))
		for (let k = k0; k < k1; k++) {
			let v = k < mid ? (k - lo) / (mid - lo) : (hi - k) / (hi - mid)
			w[k - k0] = v < 0 ? 0 : v
		}
		bank.push({ k0, w })
	}
	return bank
}

// DCT-II basis, unnormalized: row k is cos(π k (2m + 1) / 2M) over the M bands
function dctBasis (bins, nMel) {
	let c = new Float64Array(bins * nMel)
	for (let k = 0; k < bins; k++) for (let m = 0; m < nMel; m++) c[k * nMel + m] = Math.cos(Math.PI * k * (2 * m + 1) / (2 * nMel))
	return c
}

/**
 * @param {Float32Array} data — mono PCM frame (power-of-2 length)
 * @param {object} opts — { fs, bins=13, nMel=40, fmin=0, fmax=fs/2 }
 * @returns {Float32Array} cepstral coefficients
 */
export default function mfcc (data, { fs = 44100, bins = 13, nMel = 40, fmin = 0, fmax = fs / 2 } = {}) {
	let [re, im] = fft(data)
	let n = data.length, half = n / 2
	let bank = memo(`bank ${nMel} ${n} ${fs} ${fmin} ${fmax}`, () => melBank(nMel, half + 1, fs, n, fmin, fmax))
	let dct = memo(`dct ${bins} ${nMel}`, () => dctBasis(bins, nMel))
	let logMel = memo(`log ${nMel}`, () => new Float64Array(nMel))

	for (let m = 0; m < nMel; m++) {
		let { k0, w } = bank[m], e = 0
		for (let j = 0; j < w.length; j++) { let k = k0 + j; e += (re[k] * re[k] + im[k] * im[k]) * w[j] }
		logMel[m] = Math.log(e + 1e-10)
	}

	let out = new Float32Array(bins)
	for (let k = 0; k < bins; k++) {
		let sum = 0
		for (let m = 0; m < nMel; m++) sum += logMel[m] * dct[k * nMel + m]
		out[k] = sum
	}
	return out
}
