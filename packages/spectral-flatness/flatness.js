// Spectral flatness — geometric/arithmetic mean ratio of the POWER spectrum, 0 (tonal) … 1 (noise).
// Peeters 2004 §6.3; power spectrum per the audio-core fix.
export default function flatness (mag) {
	let n = mag.length
	if (!n) return 0
	let logSum = 0, sum = 0
	for (let k = 0; k < n; k++) {
		let p = mag[k] * mag[k] + 1e-12
		logSum += Math.log(p); sum += p
	}
	return Math.exp(logSum / n) / (sum / n)
}
