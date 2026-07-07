// Spectral crest — peak-to-mean ratio of the POWER spectrum (≥1; high = peaky/tonal).
// FFmpeg aspectralstats / Peeters 2004 §6.4.
export default function crest (mag) {
	let n = mag.length
	if (!n) return 0
	let max = 0, sum = 0
	for (let k = 0; k < n; k++) {
		let p = mag[k] * mag[k]
		if (p > max) max = p
		sum += p
	}
	return sum > 0 ? max / (sum / n) : 0
}
