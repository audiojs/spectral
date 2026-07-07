// Spectral freeze — sustain one spectral frame: keep its magnitudes, advance phase by
// the bin frequency each hop with per-bin random phase dispersion against buzz
// (Ableton Spectral Time / Paulstretch-freeze class).

import { fft, ifft } from 'fourier-transform'

/**
 * @param {Float32Array} data — mono PCM
 * @param {object} opts — { fs=44100, at=0 (s, frame to freeze), duration=2 (s of output),
 *   frameSize=4096, seed=1 }
 * @returns {Float32Array}
 */
export default function freeze (data, { fs = 44100, at = 0, duration = 2, frameSize = 4096, seed = 1 } = {}) {
	let half = frameSize / 2, hop = frameSize / 4
	let win = new Float64Array(frameSize)
	for (let i = 0; i < frameSize; i++) win[i] = 0.5 - 0.5 * Math.cos(2 * Math.PI * i / frameSize)
	let start = Math.min(Math.max(0, Math.round(at * fs)), Math.max(0, data.length - frameSize))
	let buf = new Float64Array(frameSize)
	for (let i = 0; i < frameSize; i++) buf[i] = (data[start + i] || 0) * win[i]
	let fr = fft(buf)
	let mag = new Float64Array(half + 1)
	for (let k = 0; k <= half; k++) mag[k] = Math.hypot(fr[0][k], fr[1][k])

	let n = Math.round(duration * fs)
	let out = new Float64Array(n + frameSize)
	let norm = new Float64Array(n + frameSize)
	let s = seed >>> 0 || 1
	let rand = () => { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff }
	let phase = new Float64Array(half + 1)
	for (let k = 0; k <= half; k++) phase[k] = rand() * 2 * Math.PI
	let re = new Float64Array(half + 1), im = new Float64Array(half + 1)
	for (let pos = 0; pos < n; pos += hop) {
		for (let k = 0; k <= half; k++) {
			phase[k] += 2 * Math.PI * k * hop / frameSize + (rand() - 0.5) * 0.3
			re[k] = mag[k] * Math.cos(phase[k])
			im[k] = mag[k] * Math.sin(phase[k])
		}
		let y = ifft(re, im)
		for (let i = 0; i < frameSize; i++) {
			out[pos + i] += y[i] * win[i]
			norm[pos + i] += win[i] * win[i]
		}
	}
	let res = new Float32Array(n)
	for (let i = 0; i < n; i++) res[i] = norm[i] > 1e-9 ? out[i] / norm[i] : 0
	return res
}
