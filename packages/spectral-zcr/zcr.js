// Zero-crossing rate — fraction of consecutive samples whose sign bit differs
// (librosa zero_crossing_rate: np.signbit diff, so −0 counts negative; mean over frame length).
export default function zcr (x) {
	let n = x.length, c = 0
	if (!n) return 0
	let prev = x[0] < 0 || (x[0] === 0 && 1 / x[0] < 0)
	for (let i = 1; i < n; i++) {
		const v = x[i], cur = v < 0 || (v === 0 && 1 / v < 0)
		if (cur !== prev) c++
		prev = cur
	}
	return c / n
}
