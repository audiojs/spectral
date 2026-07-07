// Spectral flux — frame-to-frame magnitude change. Euclidean by default (FFmpeg aspectralstats);
// { rectified: true } sums positive differences only (onset-detection flavour).
export default function flux (mag, prev, { rectified = false } = {}) {
	if (!prev) return 0
	let v = 0
	for (let k = 0; k < mag.length; k++) {
		let d = mag[k] - (prev[k] || 0)
		if (rectified) { if (d > 0) v += d }
		else v += d * d
	}
	return rectified ? v : Math.sqrt(v)
}
