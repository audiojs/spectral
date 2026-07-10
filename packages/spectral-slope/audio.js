// stat manifest — spectral slope (linear-regression tilt of the magnitude spectrum).
// The kernel is frame-domain (magnitude spectrum in) — the stat form frames the
// mono-folded signal (hann, 2048/512), runs the kernel per frame, and averages.

import kernelFn from './slope.js'
import rfft from 'fourier-transform'

const N = 2048, HOP = 512

export const slope = {
	stat: 'slope',
	compute: (channels, { sampleRate, ...opts }) => {
		const n = channels[0]?.length || 0
		if (!n) return 0
		// mono fold
		const mono = new Float32Array(n)
		for (const ch of channels) for (let i = 0; i < n; i++) mono[i] += ch[i] / channels.length
		const win = new Float64Array(N)
		for (let i = 0; i < N; i++) win[i] = 0.5 - 0.5 * Math.cos(2 * Math.PI * i / (N - 1))
		const frame = new Float64Array(N)
		let acc = 0, cnt = 0
		for (let off = 0; off + N <= n; off += HOP) {
			for (let i = 0; i < N; i++) frame[i] = mono[off + i] * win[i]
			const mag = rfft(frame)
			const v = kernelFn(mag, { fs: sampleRate, n: N, ...opts })
			if (Number.isFinite(v)) { acc += v; cnt++ }
		}
		return cnt ? acc / cnt : 0
	},
}
