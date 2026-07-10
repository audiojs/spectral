// stat manifest — spectral contrast (per-band peak/valley dB spread; librosa class).
// Frame-domain kernel returning one value per band — the stat form frames the
// mono-folded signal and averages each band across frames → Float32Array[bands+1].

import kernelFn from './contrast.js'
import rfft from 'fourier-transform'

const N = 2048, HOP = 512

export const contrast = {
	stat: 'contrast',
	compute: (channels, { sampleRate, ...opts }) => {
		const n = channels[0]?.length || 0
		if (!n) return new Float32Array(0)
		const mono = new Float32Array(n)
		for (const ch of channels) for (let i = 0; i < n; i++) mono[i] += ch[i] / channels.length
		const win = new Float64Array(N)
		for (let i = 0; i < N; i++) win[i] = 0.5 - 0.5 * Math.cos(2 * Math.PI * i / (N - 1))
		const frame = new Float64Array(N)
		let acc = null, cnt = 0
		for (let off = 0; off + N <= n; off += HOP) {
			for (let i = 0; i < N; i++) frame[i] = mono[off + i] * win[i]
			const v = kernelFn(rfft(frame), { fs: sampleRate, n: N, ...opts })
			if (!acc) acc = new Float64Array(v.length)
			for (let b = 0; b < v.length; b++) acc[b] += v[b]
			cnt++
		}
		if (!acc) return new Float32Array(0)
		return Float32Array.from(acc, x => x / cnt)
	},
}
