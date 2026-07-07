// LTAS — long-term average spectrum via Welch's method: hann-windowed frames,
// power-averaged, returned as RMS magnitude per bin. The adaptive-EQ / match-EQ substrate.

import { fft } from 'fourier-transform'

export default function ltas (data, { frameSize = 4096, hop = frameSize / 2 } = {}) {
	let half = frameSize / 2
	let win = new Float64Array(frameSize)
	for (let i = 0; i < frameSize; i++) win[i] = 0.5 - 0.5 * Math.cos(2 * Math.PI * i / (frameSize - 1))

	let acc = new Float64Array(half + 1)
	let frames = 0
	let buf = new Float64Array(frameSize)
	for (let pos = 0; pos + frameSize <= data.length; pos += hop) {
		for (let i = 0; i < frameSize; i++) buf[i] = data[pos + i] * win[i]
		let [re, im] = fft(buf)
		for (let k = 0; k <= half; k++) acc[k] += re[k] * re[k] + im[k] * im[k]
		frames++
	}
	if (!frames) return new Float32Array(half + 1)

	let out = new Float32Array(half + 1)
	for (let k = 0; k <= half; k++) out[k] = Math.sqrt(acc[k] / frames)
	return out
}
