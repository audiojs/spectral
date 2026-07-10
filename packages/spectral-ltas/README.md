# @audio/spectral-ltas

> LTAS — long-term average spectrum via Welch's method: hann-windowed frames, power-averaged, returned as RMS magnitude per bin. Adaptive-EQ / match-EQ substrate.

`npm install @audio/spectral-ltas`

```js
import ltas from '@audio/spectral-ltas'

let curve = ltas(data, { frameSize: 4096, hop: 2048 })
// Float32Array, length frameSize/2 + 1 — RMS magnitude per FFT bin, bin k ↔ k·fs/frameSize
```

Options: - `frameSize` — Welch analysis window, samples (default 4096) · `hop` — frame hop, samples (default `frameSize / 2`, i.e. 50% overlap)

Also exported as an `audio.js` stat manifest (`./audio` — mono-folds multichannel input, then runs the kernel).

Part of [@audio/spectral](https://github.com/audiojs/spectral).
