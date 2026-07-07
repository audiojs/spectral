// Spectral centroid — magnitude-weighted mean frequency (Hz). Peeters 2004 §6.1.
export default function centroid (mag, { fs = 44100, n = 2 * (mag.length - 1) } = {}) {
	let num = 0, den = 0
	for (let k = 0; k < mag.length; k++) { num += (k * fs / n) * mag[k]; den += mag[k] }
	return den > 0 ? num / den : 0
}
