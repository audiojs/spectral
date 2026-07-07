// Spectral slope — least-squares linear regression of magnitude over frequency. Peeters 2004 §6.6.
export default function slope (mag, { fs = 44100, n = 2 * (mag.length - 1) } = {}) {
	let N = mag.length
	if (N < 2) return 0
	let sf = 0, sm = 0, sfm = 0, sff = 0
	for (let k = 0; k < N; k++) {
		let f = k * fs / n
		sf += f; sm += mag[k]; sfm += f * mag[k]; sff += f * f
	}
	let den = N * sff - sf * sf
	return den !== 0 ? (N * sfm - sf * sm) / den : 0
}
