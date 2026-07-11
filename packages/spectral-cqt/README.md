# @audio/spectral-cqt

> Constant-Q transform — log-frequency spectrum, one Q-scaled window per bin, direct DFT evaluation (Brown 1991).

`npm install @audio/spectral-cqt`

```js
import cqt from '@audio/spectral-cqt'

let { freqs, mag } = cqt(data, { fs: 44100, fmin: 55, binsPerOctave: 12, octaves: 6, at: 0 })
// freqs: Float32Array(72) — bin center frequencies, log-spaced (12/octave over 6 octaves from 55 Hz)
// mag:   Float32Array(72) — magnitude at each bin, analyzed at t = `at` seconds
```

Raw-signal in: each bin gets its own window, long at low frequencies (better frequency
resolution, worse time resolution) and short at high frequencies (Q = `1 / (2^(1/binsPerOctave) − 1)`
constant across bins) — the point of a constant-Q analysis over a fixed-window FFT. One call
returns one analysis frame centered at `at`; call again per frame for a CQT-gram.

Options: - `fs` — sample rate (default 44100, Hz) · `fmin` — lowest bin center (default 55, Hz — A1) · `binsPerOctave` — bins per octave (default 12 — semitone resolution) · `octaves` — octave span (default 6) · `at` — analysis center, seconds (default 0)

Use when: log-frequency display (spectrogram with musical pitch spacing), chroma/CQT feature
substrates — not a drop-in FFT replacement (no fixed hop/frame grid, no `ifft`).

Part of [@audio/spectral](https://github.com/audiojs/spectral).
