# @audio/spectral-freeze

> Spectral freeze — sustain one spectral frame indefinitely, phase advanced per bin with random dispersion against buzz (Ableton Spectral Time / Paulstretch-freeze class).

`npm install @audio/spectral-freeze`

```js
import freeze from '@audio/spectral-freeze'

let out = freeze(data, { fs: 44100, at: 1.5, duration: 4, frameSize: 4096 })
// Float32Array(round(duration * fs)) — new audio, the frame at `at` sustained for `duration` seconds
```

Takes the magnitude spectrum of one analysis frame (windowed at `at`), then resynthesizes it hop
by hop: each bin's phase advances by its own bin frequency plus a small random offset (seeded by
`seed`), so the sustained tone doesn't beat/buzz from perfectly periodic phase increments the way
a naive frame-repeat would.

Options: - `fs` — sample rate (default 44100, Hz) · `at` — frame to freeze, seconds (default 0) · `duration` — output length, seconds (default 2) · `frameSize` — STFT frame (default 4096) · `seed` — phase-dispersion RNG seed (default 1, deterministic)

Use when: pad/drone textures from a single frozen instant, spectral "hold" effects — not a
general time-stretch (use `@audio/spectral-pvoc` / `@audio/stretch` for that).

Part of [@audio/spectral](https://github.com/audiojs/spectral).
