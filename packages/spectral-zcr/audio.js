// stat manifest — zero-crossing rate (per-frame sign-change fraction; librosa zero_crossing_rate class).
// The kernel is time-domain (raw samples in, no window/FFT — librosa does not window zcr) — the
// stat form frames the mono-folded signal (2048/512), runs the kernel per frame, and averages.

import kernelFn from './zcr.js'

const N = 2048, HOP = 512

export const zcr = {
	stat: 'zcr',
	compute: (channels, { sampleRate, ...opts }) => {
		const n = channels[0]?.length || 0
		if (!n) return 0
		// mono fold
		const mono = new Float32Array(n)
		for (const ch of channels) for (let i = 0; i < n; i++) mono[i] += ch[i] / channels.length
		const frame = new Float32Array(N)
		const kOpts = { fs: sampleRate, ...opts }
		let acc = 0, cnt = 0
		for (let off = 0; off + N <= n; off += HOP) {
			for (let i = 0; i < N; i++) frame[i] = mono[off + i]
			const v = kernelFn(frame, kOpts)
			if (Number.isFinite(v)) { acc += v; cnt++ }
		}
		return cnt ? acc / cnt : 0
	},
}
