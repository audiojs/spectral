// Constant-Q transform — log-spaced bins, per-bin Q-scaled windowed DFT (direct
// evaluation; Brown 1991). Offline analysis: chroma_cqt substrate, log-frequency display.

/**
 * @param {Float32Array} data — mono PCM
 * @param {object} opts — { fs=44100, fmin=55, binsPerOctave=12, octaves=6, at=0 (s, analysis center) }
 * @returns {{ freqs: Float32Array, mag: Float32Array }}
 */
export default function cqt (data, { fs = 44100, fmin = 55, binsPerOctave = 12, octaves = 6, at = 0 } = {}) {
	let nBins = binsPerOctave * octaves
	let Q = 1 / (2 ** (1 / binsPerOctave) - 1)
	let freqs = new Float32Array(nBins)
	let mag = new Float32Array(nBins)
	let center = Math.round(at * fs)
	for (let b = 0; b < nBins; b++) {
		let f = fmin * 2 ** (b / binsPerOctave)
		freqs[b] = f
		let N = Math.min(data.length, Math.round(Q * fs / f))
		let start = Math.max(0, Math.min(data.length - N, center - (N >> 1)))
		let reS = 0, imS = 0
		for (let i = 0; i < N; i++) {
			let w = 0.5 - 0.5 * Math.cos(2 * Math.PI * i / N)
			let ph = 2 * Math.PI * f * i / fs
			let x = data[start + i] * w
			reS += x * Math.cos(ph)
			imS += x * Math.sin(ph)
		}
		mag[b] = 2 * Math.hypot(reS, imS) / N
	}
	return { freqs, mag }
}
