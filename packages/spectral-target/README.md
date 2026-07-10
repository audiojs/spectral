# @audio/spectral-target

> Target-curve library — per-content-type LTAS targets + the deviation math that turns
> (measured LTAS, target) into an EQ correction curve.

Curve data + pure math, no DSP kernel of its own. Companion of [`@audio/spectral-ltas`](https://github.com/audiojs/spectral/tree/main/packages/spectral-ltas)
(returns linear RMS magnitude per FFT bin, length `half+1`, bin `k ↔ k·fs/frameSize`) — every
curve here renders onto that exact bin grid.

```js
import target, { deviation, smooth } from '@audio/spectral-target'

let curve = target('speech', { fs: 44100, bins: 2049 })      // Float32Array dB per bin
let custom = target([[100, 0], [1000, -10], [8000, -25]], { fs, bins })

let correction = deviation(measuredLtas, curve, { fs })      // target − measured, smoothed, clamped
let smoothed = smooth(someDbCurve, { fs, oct: 1 / 3 })
```

| export | role |
|---|---|
| `target(name \| anchors, { fs=44100, bins=2049 })` | render `'speech'` \| `'music'` \| `'pink'` \| `'voice-music'` or a custom `[[freqHz, dB], …]` table to a dB-per-bin curve |
| `deviation(measuredLtas, targetDb, { fs=44100, smoothOct=1/3, clamp=12 })` | `target − measured` → EQ correction curve, dB per bin |
| `smooth(db, { fs=44100, oct=1/3 })` | fractional-octave smoothing of any per-bin dB curve |

## Sources

**`'speech'`** — Byrne, Dillon, Tran et al. 1994, ["An international comparison of long-term
average speech spectra,"](https://doi.org/10.1121/1.410152) *JASA* 96(4):2108–2120. Table II,
p.2116, **"Combined"** column (male+female average across 17 samples / 12 languages, normalized
to 70 dB SPL overall) — the paper's own "universal LTASS." 1/3-octave, 100 Hz–10 kHz (the paper
also tabulates 63/80 Hz and 12.5/16 kHz; this module uses its stated 100 Hz–10 kHz range):

| Hz | 100 | 125 | 160 | 200 | 250 | 315 | 400 | 500 | 630 | 800 | 1000 | 1250 | 1600 | 2000 | 2500 | 3150 | 4000 | 5000 | 6300 | 8000 | 10000 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| dB | 54.4 | 57.7 | 56.8 | 60.2 | 60.3 | 59.0 | 62.1 | 62.1 | 60.5 | 56.8 | 53.7 | 53.0 | 52.0 | 48.7 | 48.1 | 46.8 | 45.6 | 44.5 | 44.3 | 43.4 | 43.4 |

**`'music'`** — Pestana, Ma, Reiss, Barbosa & Black 2013, ["Spectral Characteristics of Popular
Commercial Recordings 1950–2010,"](https://aes.org/publications/elibrary-page/?id=17010) AES
Convention 135, Paper 8960, p.5: mean LTAS of 772 UK/US number-one hits shows *"a linearly
decaying distribution of around 5 dB per octave between 100 and 4000 Hz."* (Quoted phrase
verified via Schedin, ["Target Spectrums For Mastering,"](https://www.diva-portal.org/smash/get/diva2:1557981/FULLTEXT01.pdf)
Luleå University of Technology bachelor's thesis, 2021, p.7, which cites the same page — the
primary AES paper wasn't directly reachable while sourcing this module; flagging the path
honestly rather than implying a direct read.) The paper also describes a low cut "around 60 Hz"
and a slope that gets "gradually steeper" above 4 kHz, but publishes no numeric rate for either,
so neither is modeled with an invented number — see Edge policies below for what happens instead.

**`'pink'`** — analytic, −10·log10(2) ≈ −3.0103 dB/octave (equal energy per octave). Definitional,
no citation needed: the neutral tonal-balance reference and the exact test anchor.

**`'voice-music'`** — **audiojs-defined convention, not literature.** Bin-wise `max(speech, music)
− 6 dB`: leans on whichever source curve is louder at each frequency (protects the speech
formant bands from being swamped by a music-shaped target, while still tracking music's wider
low/high extension), pulled down 6 dB so the blend doesn't read hotter than either source. Pragmatic choice for voice-over-music content; revisit if a cited source surfaces.

## Edge policies

- **Below the lowest anchor**: hold its value (flat). Applies to bin 0 / DC — never returns
  `-Infinity`.
- **Above the highest anchor**: continue the last segment's slope, not flat.
- **Interpolation**: log-frequency linear between anchors (exact for `'pink'`, since a
  log-linear function sampled and re-interpolated log-linearly reproduces itself exactly at
  every bin — not an approximation).
- **Normalization**: every curve `target()` returns is mean-normalized to 0 dB over its own
  *defined band* = `[lowest anchor, min(highest anchor, 0.45·fs)]`. For `'music'` that band is
  `[100, 4000]` Hz — the paper's own measured range — so `deviation()` guidance built from a
  *rendered* `'music'` curve is naturally quietest right at that band's edges; combine with
  `'pink'` or a custom table if you need sourced full-band guidance.
- **`deviation()`'s own band** is fixed at `[20 Hz, 0.45·fs]` regardless of which target was
  rendered (it has no anchor metadata to read — `targetDb` is a plain array by design). Outside
  that band, correction tapers to exactly 0 over a 0.5-octave raised-cosine ramp at each edge
  (no ultrasonic/DC boosts, and no discontinuous jump in the correction curve).

## Manifest

No `audio.js` — this atom is not a processor and not a per-signal stat, it's curve data + pure
math (same category as `@audio/spectral-pvoc`). A `balance` stat (deviation vs a named target)
was considered but doesn't drop out cleanly: computing it from raw `channels` would mean either
depending on `@audio/spectral-ltas` (no `@audio/spectral-*` package depends on a sibling today)
or re-implementing Welch averaging here (duplicating it). Per the atom contract, a package
without a manifest is still an atom — consumers that already hold a measured LTAS call
`deviation()` directly.

## Consumers

Mix Analyser's spectral-balance report, Speech/Music Enhancer's adaptive EQ, Auto-Chain's
Stage-2 EQ, and Match-by-Reference (where the "target" is a reference track's own measured LTAS
via `@audio/spectral-ltas`, passed straight into `deviation()` in place of a preset). MIT.
