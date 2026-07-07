// Spectral contrast — per-octave-band difference between peak and valley energy
// percentiles of the magnitude spectrum (Jiang 2002; librosa/essentia parity).

/**
 * @param {Float32Array} mag — magnitude spectrum
 * @param {object} opts — { fs=44100, n=2·(mag.length−1), fmin=200, bands=6, quantile=0.2 }
 * @returns {Float32Array} contrast in dB per band
 */
export default function contrast (mag, { fs = 44100, n = 2 * (mag.length - 1), fmin = 200, bands = 6, quantile = 0.2 } = {}) {
	let out = new Float32Array(bands)
	for (let b = 0; b < bands; b++) {
		let lo = fmin * 2 ** b, hi = fmin * 2 ** (b + 1)
		let k0 = Math.max(1, Math.round(lo * n / fs)), k1 = Math.min(mag.length - 1, Math.round(hi * n / fs))
		if (k1 <= k0 + 2) break
		let seg = Array.from(mag.subarray(k0, k1 + 1)).sort((a, b) => a - b)
		let q = Math.max(1, Math.round(seg.length * quantile))
		let valley = 0, peak = 0
		for (let i = 0; i < q; i++) { valley += seg[i]; peak += seg[seg.length - 1 - i] }
		out[b] = 20 * Math.log10((peak + 1e-12) / (valley + 1e-12))
	}
	return out
}
