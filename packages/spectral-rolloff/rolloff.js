// Spectral rolloff — frequency below which `p` of total spectral energy lies (default 85%).
export default function rolloff (mag, { fs = 44100, n = 2 * (mag.length - 1), p = 0.85 } = {}) {
	let total = 0
	for (let k = 0; k < mag.length; k++) total += mag[k] * mag[k]
	if (total <= 0) return 0
	let acc = 0
	for (let k = 0; k < mag.length; k++) {
		acc += mag[k] * mag[k]
		if (acc >= p * total) return k * fs / n
	}
	return (mag.length - 1) * fs / n
}
