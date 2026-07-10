// Target-curve library — per-content-type long-term-average-spectrum (LTAS) targets, plus the
// deviation math that turns (measured LTAS, target) into an EQ correction curve. Curve data +
// pure math, no DSP kernel of its own — companion of @audio/spectral-ltas, which returns linear
// RMS magnitude per FFT bin (length half+1, bin k ↔ k·fs/frameSize); every curve here renders
// onto that exact bin grid. See README for full citations, the voice-music convention, and edge
// policies.

// ── cited anchor tables (Hz, dB — relative; target() mean-normalizes each to 0 dB over its own
// defined band, so the reference level chosen below is arbitrary and cancels out) ──

// Byrne, Dillon, Tran et al. 1994, "An international comparison of long-term average speech
// spectra," JASA 96(4):2108–2120 — Table II, p.2116, "Combined" column (male+female average
// across 17 samples/12 languages, normalized to 70 dB SPL overall level). 1/3-octave, the
// paper's own 100 Hz–10 kHz range (it also tabulates 63/80 Hz and 12.5/16 kHz, omitted here).
const SPEECH = [
	[100, 54.4], [125, 57.7], [160, 56.8], [200, 60.2], [250, 60.3],
	[315, 59.0], [400, 62.1], [500, 62.1], [630, 60.5], [800, 56.8],
	[1000, 53.7], [1250, 53.0], [1600, 52.0], [2000, 48.7], [2500, 48.1],
	[3150, 46.8], [4000, 45.6], [5000, 44.5], [6300, 44.3], [8000, 43.4],
	[10000, 43.4],
]

// Pestana, Ma, Reiss, Barbosa & Black 2013, "Spectral Characteristics of Popular Commercial
// Recordings 1950–2010," AES Convention 135, Paper 8960, p.5: mean LTAS of 772 number-one
// UK/US hits shows "a linearly decaying distribution of around 5 dB per octave between 100 and
// 4000 Hz" (quoted phrase verified via Schedin, "Target Spectrums For Mastering," Luleå
// University of Technology bachelor's thesis, 2021, p.7, which cites the same page — the
// primary AES paper was not directly reachable while authoring this module). The paper
// additionally describes a low cut "around 60 Hz" and a slope that gets "gradually steeper"
// above 4 kHz, but publishes no numeric rate for either segment, so neither is modeled with an
// invented number: below 100 Hz this module's standard low-anchor policy holds the 100 Hz level
// flat (reads as the paper's own "low-shelf plateau" language), and above 4 kHz it continues the
// one cited −5 dB/oct rate (standard high-anchor policy) rather than guess a steeper figure.
const MUSIC_DB_PER_OCT = -5
const MUSIC = [[100, 0], [4000, MUSIC_DB_PER_OCT * Math.log2(4000 / 100)]]

// Analytic pink (equal energy per octave = −10·log10(2) ≈ −3.0103 dB/octave). Definitional, no
// citation needed — the neutral tonal-balance reference and the exact test anchor. Two anchors
// spanning the audible range fully determine the line (log-linear interpolation of a log-linear
// function is exact at every intermediate bin, not an approximation).
const PINK_DB_PER_OCT = -10 * Math.log10(2)
const PINK = [[20, PINK_DB_PER_OCT * Math.log2(20 / 1000)], [20000, PINK_DB_PER_OCT * Math.log2(20000 / 1000)]]

const PRESETS = { speech: SPEECH, music: MUSIC, pink: PINK }

// Defined-band ceiling as a fraction of fs — keeps normalization/deviation out of the
// ultrasonic fringe near Nyquist regardless of anchor table or measured input.
const NYQUIST_MARGIN = 0.45

// audiojs-defined convention (not literature): voice-over-music content leans on whichever of
// speech/music is louder at each frequency — protects intelligibility in the speech formant
// bands while still tracking music's wider-band energy — pulled down 6 dB so the blend doesn't
// read hotter than either source target. See README.
const VOICE_MUSIC_HEADROOM_DB = 6

/**
 * Render an anchor table onto a bin grid: log-frequency linear interpolation between control
 * points. Below the lowest anchor, holds its value (bin 0 / DC included — never −Infinity).
 * Above the highest anchor, continues the last segment's slope.
 */
function renderAnchors (anchors, fs, bins) {
	let n = 2 * (bins - 1)
	let out = new Float32Array(bins)
	let last = anchors.length - 1
	let i = 0
	for (let k = 0; k < bins; k++) {
		let f = k * fs / n
		if (f <= anchors[0][0]) { out[k] = anchors[0][1]; continue }
		if (f >= anchors[last][0]) {
			if (last === 0) out[k] = anchors[0][1]
			else {
				let [f0, d0] = anchors[last - 1]
				let [f1, d1] = anchors[last]
				let slope = (d1 - d0) / Math.log2(f1 / f0)
				out[k] = d1 + slope * Math.log2(f / f1)
			}
			continue
		}
		while (i < last - 1 && anchors[i + 1][0] < f) i++
		let [f0, d0] = anchors[i]
		let [f1, d1] = anchors[i + 1]
		let t = Math.log2(f / f0) / Math.log2(f1 / f0)
		out[k] = d0 + t * (d1 - d0)
	}
	return out
}

// Arithmetic mean of `db[k]` over bins whose frequency falls in [lo, hi] (inclusive).
function bandMean (db, fs, n, lo, hi) {
	let k0 = Math.max(0, Math.ceil(lo * n / fs))
	let k1 = Math.min(db.length - 1, Math.floor(hi * n / fs))
	let sum = 0, cnt = 0
	for (let k = k0; k <= k1; k++) { sum += db[k]; cnt++ }
	return cnt ? sum / cnt : 0
}

function normalizeToBand (db, fs, n, lo, hi) {
	let m = bandMean(db, fs, n, lo, hi)
	let out = new Float32Array(db.length)
	for (let k = 0; k < db.length; k++) out[k] = db[k] - m
	return out
}

/**
 * @param {'speech'|'music'|'pink'|'voice-music'|Array<[number,number]>} name — preset name or a
 *   custom `[[freqHz, dB], ...]` anchor table (ascending frequency)
 * @param {object} opts — { fs=44100, bins=2049 }
 * @returns {Float32Array} dB per bin, mean-normalized to 0 dB over the target's defined band
 *   (= [lowest anchor, min(highest anchor, 0.45·fs)])
 */
export default function target (name, opts = {}) {
	if (name === 'voice-music') return renderVoiceMusic(opts)

	let anchors = Array.isArray(name) ? name : PRESETS[name]
	if (!anchors) throw new Error(`@audio/spectral-target: unknown target "${name}"`)

	let { fs = 44100, bins = 2049 } = opts
	let n = 2 * (bins - 1)
	let raw = renderAnchors(anchors, fs, bins)
	let lo = anchors[0][0]
	let hi = Math.min(anchors[anchors.length - 1][0], NYQUIST_MARGIN * fs)
	return normalizeToBand(raw, fs, n, lo, hi)
}

function renderVoiceMusic ({ fs = 44100, bins = 2049 } = {}) {
	let n = 2 * (bins - 1)
	let s = renderAnchors(SPEECH, fs, bins)
	let m = renderAnchors(MUSIC, fs, bins)
	let raw = new Float32Array(bins)
	for (let k = 0; k < bins; k++) raw[k] = Math.max(s[k], m[k]) - VOICE_MUSIC_HEADROOM_DB
	let lo = Math.min(SPEECH[0][0], MUSIC[0][0])
	let hi = Math.min(Math.max(SPEECH[SPEECH.length - 1][0], MUSIC[MUSIC.length - 1][0]), NYQUIST_MARGIN * fs)
	return normalizeToBand(raw, fs, n, lo, hi)
}

/**
 * Fractional-octave (log-frequency) smoothing of a per-bin dB curve: each bin becomes the mean
 * of all bins within ±oct/2 octaves of it. Exported because match-EQ reuses it directly.
 * @param {Float32Array} db
 * @param {object} opts — { fs=44100, oct=1/3 }
 */
export function smooth (db, { fs = 44100, oct = 1 / 3 } = {}) {
	let bins = db.length
	let n = 2 * (bins - 1)
	let half = oct / 2
	let out = new Float32Array(bins)
	out[0] = db[0]
	for (let k = 1; k < bins; k++) {
		let f = k * fs / n
		let k0 = Math.max(1, Math.round(f / 2 ** half * n / fs))
		let k1 = Math.min(bins - 1, Math.round(f * 2 ** half * n / fs))
		let sum = 0
		for (let j = k0; j <= k1; j++) sum += db[j]
		out[k] = sum / (k1 - k0 + 1)
	}
	return out
}

// Fixed-width (in octaves) raised-cosine fade at each edge of [lo, hi]; 0 at and beyond the
// edge, 1 once TAPER_OCT octaves inside — so an EQ correction curve never jumps discontinuously
// to 0 at the band boundary, and is exactly 0 at and outside it (no ultrasonic/DC boosts).
const TAPER_OCT = 0.5

function rcos (x) {
	return 0.5 - 0.5 * Math.cos(Math.PI * Math.min(1, Math.max(0, x)))
}

function taper (f, lo, hi) {
	if (f <= lo || f >= hi) return 0
	let g = 1
	let loEdge = lo * 2 ** TAPER_OCT
	if (f < loEdge) g = Math.min(g, rcos(Math.log2(f / lo) / TAPER_OCT))
	let hiEdge = hi / 2 ** TAPER_OCT
	if (f > hiEdge) g = Math.min(g, rcos(Math.log2(hi / f) / TAPER_OCT))
	return g
}

/**
 * (measured LTAS, target) → EQ correction curve. Both sides are independently band-limited and
 * mean-normalized over [20 Hz, 0.45·fs] (so absolute level never leaks into tone), subtracted,
 * fractional-octave smoothed, clamped, then tapered to exactly 0 at and beyond the band edges.
 * @param {Float32Array} measuredLtas — linear RMS magnitude per bin, e.g. from @audio/spectral-ltas
 * @param {Float32Array} targetDb — dB per bin, e.g. from target()
 * @param {object} opts — { fs=44100, smoothOct=1/3, clamp=12 }
 * @returns {Float32Array} correction dB per bin (target − measured), finite everywhere
 */
export function deviation (measuredLtas, targetDb, { fs = 44100, smoothOct = 1 / 3, clamp = 12 } = {}) {
	let bins = targetDb.length
	let n = 2 * (bins - 1)
	let lo = 20, hi = NYQUIST_MARGIN * fs

	let measuredDb = new Float32Array(bins)
	for (let k = 0; k < bins; k++) measuredDb[k] = 20 * Math.log10(Math.max(measuredLtas[k] || 0, 1e-12))

	let t = normalizeToBand(targetDb, fs, n, lo, hi)
	let m = normalizeToBand(measuredDb, fs, n, lo, hi)
	let raw = new Float32Array(bins)
	for (let k = 0; k < bins; k++) raw[k] = t[k] - m[k]

	let sm = smooth(raw, { fs, oct: smoothOct })

	let out = new Float32Array(bins)
	for (let k = 0; k < bins; k++) {
		let f = k * fs / n
		let g = taper(f, lo, hi)
		if (g > 0) out[k] = Math.max(-clamp, Math.min(clamp, sm[k])) * g
	}
	return out
}
