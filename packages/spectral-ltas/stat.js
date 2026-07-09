// stat manifest — long-term average spectrum (mean magnitude spectrum of the take).

import ltasFn from './ltas.js'

export const ltas = {
	stat: 'ltas',
	compute: (channels, { sampleRate, ...opts }) => {
		const n = channels[0]?.length || 0
		if (!n) return new Float32Array(0)
		const mono = new Float32Array(n)
		for (const ch of channels) for (let i = 0; i < n; i++) mono[i] += ch[i] / channels.length
		return ltasFn(mono, opts)
	},
}
