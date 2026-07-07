// Spectral editing — apply gain to time×frequency regions (delete / attenuate / boost).
// STFT → per-region bin gains → weighted overlap-add resynthesis (COLA-normalized).
// Audacity spectral-edit / FFmpeg afftfilt class.

import { fft, ifft } from 'fourier-transform'

/**
 * @param {Float32Array} data — mono PCM
 * @param {object} opts — { fs, frameSize=2048, hop=frameSize/4,
 *   regions: [{ t0=0, t1=∞ (seconds), f0=0, f1=fs/2 (Hz), gain=0 (linear) }] }
 * @returns {Float32Array} edited copy
 */
export default function edit (data, { fs = 44100, frameSize = 2048, hop = frameSize / 4, regions = [] } = {}) {
	let half = frameSize / 2
	let win = new Float64Array(frameSize)
	for (let i = 0; i < frameSize; i++) win[i] = 0.5 - 0.5 * Math.cos(2 * Math.PI * i / frameSize)

	let out = new Float64Array(data.length)
	let norm = new Float64Array(data.length)
	let buf = new Float64Array(frameSize)

	for (let pos = 0; pos + frameSize <= data.length + frameSize - hop; pos += hop) {
		for (let i = 0; i < frameSize; i++) buf[i] = (data[pos + i] || 0) * win[i]
		let [re, im] = fft(buf)

		let tMid = (pos + frameSize / 2) / fs
		for (let r of regions) {
			if (tMid < (r.t0 ?? 0) || tMid > (r.t1 ?? Infinity)) continue
			let k0 = Math.max(0, Math.floor((r.f0 ?? 0) * frameSize / fs))
			let k1 = Math.min(half, Math.ceil((r.f1 ?? fs / 2) * frameSize / fs))
			let g = r.gain ?? 0
			for (let k = k0; k <= k1; k++) { re[k] *= g; im[k] *= g }
		}

		let frame = ifft(re, im)
		for (let i = 0; i < frameSize; i++) {
			let j = pos + i
			if (j >= data.length) break
			out[j] += frame[i] * win[i]
			norm[j] += win[i] * win[i]
		}
	}

	let res = new Float32Array(data.length)
	for (let i = 0; i < data.length; i++) res[i] = norm[i] > 1e-9 ? out[i] / norm[i] : data[i]
	return res
}
