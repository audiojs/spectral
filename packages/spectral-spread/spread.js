// Spectral spread — magnitude-weighted standard deviation around the centroid (Hz). Peeters 2004 §6.1.
export default function spread (mag, { fs = 44100, n = 2 * (mag.length - 1) } = {}) {
	let num = 0, den = 0
	for (let k = 0; k < mag.length; k++) { num += (k * fs / n) * mag[k]; den += mag[k] }
	if (den <= 0) return 0
	let c = num / den, v = 0
	for (let k = 0; k < mag.length; k++) { let d = k * fs / n - c; v += d * d * mag[k] }
	return Math.sqrt(v / den)
}
