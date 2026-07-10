# @audio/spectral

> Spectral features — centroid, spread, flatness, rolloff, flux, slope, crest, zcr. All planned.

Also planned: `spectral-mfcc` (librosa parity), `spectral-ltas` (adaptive/match-EQ substrate), `spectral-edit` (time×frequency region editing — Audacity/afftfilt class).

Reference implementations exist in the `audio` package stats (`centroid`, `flatness`, `crest` shipped there — extract, don't rewrite). Parity: FFmpeg `aspectralstats`, MIREX spectral features. LTAS / spectral-stats for the Auto-Chain analyser build on these atoms.

## `spectral-target`

Target-curve library — per-content-type LTAS targets (`'speech'`, `'music'`, `'pink'`,
`'voice-music'`, or a custom anchor table) plus `deviation()`, the math that turns a measured
`spectral-ltas` curve and a target into an EQ correction curve. The adaptive-EQ data+math layer
behind Mix Analyser's spectral-balance report, Speech/Music Enhancer's adaptive EQ, Auto-Chain's
Stage-2 EQ, and Match-by-Reference.

- `'speech'` — Byrne, Dillon, Tran et al. 1994, *JASA* 96(4):2108–2120, Table II "Combined"
  column (12-language universal LTASS, 100 Hz–10 kHz).
- `'music'` — Pestana, Ma, Reiss, Barbosa & Black 2013, AES Convention 135 Paper 8960, p.5:
  "a linearly decaying distribution of around 5 dB per octave between 100 and 4000 Hz."
- `'pink'` — analytic −3.0103 dB/octave, no citation needed.
- `'voice-music'` — **audiojs-defined convention**, not literature: bin-wise `max(speech, music)
  − 6 dB`.

Edge policies: below the lowest anchor holds flat, above the highest continues the last
segment's slope, interpolation is log-frequency linear. Every curve is mean-normalized to 0 dB
over its own defined band; `deviation()` band-limits and tapers to exactly 0 outside
`[20 Hz, 0.45·fs]`. Full citations, exact tables, and the manifest-less rationale are in
`packages/spectral-target/README.md`.
