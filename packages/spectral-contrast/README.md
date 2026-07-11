# @audio/spectral-contrast

> Spectral contrast — per-octave-band peak/valley energy difference (Jiang 2002; librosa/essentia parity).

`npm install @audio/spectral-contrast`

```js
import { fft } from 'fourier-transform'
import contrast from '@audio/spectral-contrast'

let win = frame.map((x, i) => x * (0.5 - 0.5 * Math.cos(2 * Math.PI * i / (frame.length - 1))))
let [re, im] = fft(win)
let mag = re.map((r, k) => Math.hypot(r, im[k]))

let bands = contrast(mag, { fs: 44100 })   // Float32Array(6), dB per octave band
```

Each band spans `[fmin·2ᵇ, fmin·2ᵇ⁺¹)` (default `fmin = 200` Hz); within a band, contrast is
`20·log10(peak / valley)` where peak/valley are the mean of the top/bottom `quantile` of bins by
magnitude — high contrast marks tonal energy standing above a noisy floor in that band.

Options: - `fs` — sample rate (default 44100, Hz) · `n` — FFT size the spectrum was taken at (default `2 * (mag.length - 1)`, inferred from `mag.length`) · `fmin` — lowest band edge (default 200, Hz) · `bands` — number of octave bands (default 6) · `quantile` — peak/valley fraction of each band's bins (default 0.2)

Also exported as an `audio.js` stat manifest (`./audio` — mono-fold, 2048/512 hann frames, averaged over the take).

Part of [@audio/spectral](https://github.com/audiojs/spectral).
