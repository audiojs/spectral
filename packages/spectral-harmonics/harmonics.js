// Harmonic descriptors from f0 — partial amplitudes at k·f0 give inharmonicity,
// tristimulus (T1/T2/T3) and odd/even energy ratio (essentia parity).

import { fft } from 'fourier-transform'

/**
 * @param {Float32Array} data — mono PCM frame (power of 2)
 * @param {object} opts — { fs=44100, f0 (required, Hz), nHarmonics=10 }
 * @returns {{ amps: Float32Array, inharmonicity: number, tristimulus: [n,n,n], oddEven: number }}
 */
export default function harmonics (data, { fs = 44100, f0, nHarmonics = 10 } = {}) {
	if (!f0) throw new RangeError('harmonics: opts.f0 required')
	let n = data.length, half = n / 2
	let buf = new Float64Array(n)
	for (let i = 0; i < n; i++) buf[i] = data[i] * (0.5 - 0.5 * Math.cos(2 * Math.PI * i / n))
	let [re, im] = fft(buf)
	let mag = k => Math.hypot(re[k], im[k])

	let amps = new Float32Array(nHarmonics)
	let inh = 0, total = 0
	for (let h = 1; h <= nHarmonics; h++) {
		let kc = h * f0 * n / fs
		if (kc >= half - 2) break
		let b = Math.round(kc), best = b, m = 0
		for (let j = Math.max(1, b - 2); j <= Math.min(half, b + 2); j++) if (mag(j) > m) { m = mag(j); best = j }
		amps[h - 1] = m
		inh += Math.abs(best - kc) / kc * m
		total += m
	}
	let e = amps.map(a => a * a)
	let sum = e.reduce((a, b) => a + b, 0) + 1e-20
	let odd = 0, even = 0
	for (let h = 2; h <= nHarmonics; h++) {
		if (h % 2) odd += e[h - 1] || 0
		else even += e[h - 1] || 0
	}
	return {
		amps,
		inharmonicity: total > 0 ? inh / total : 0,
		tristimulus: [e[0] / sum, (e[1] + e[2] + e[3]) / sum, e.slice(4).reduce((a, b) => a + b, 0) / sum],
		oddEven: odd / (even + 1e-20),
	}
}
