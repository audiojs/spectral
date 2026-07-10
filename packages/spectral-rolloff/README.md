# @audio/spectral-rolloff

> Spectral rolloff — frequency below which `p` of total spectral energy lies (default 85%).

`npm install @audio/spectral-rolloff`

```js
import { fft } from 'fourier-transform'
import rolloff from '@audio/spectral-rolloff'

let win = frame.map((x, i) => x * (0.5 - 0.5 * Math.cos(2 * Math.PI * i / (frame.length - 1))))
let [re, im] = fft(win)
let mag = re.map((r, k) => Math.hypot(r, im[k]))

let hz = rolloff(mag, { fs: 44100, p: 0.85 })
```

Options: - `fs` — sample rate (default 44100, Hz) · `n` — FFT size the spectrum was taken at (default `2 * (mag.length - 1)`, inferred from `mag.length`) · `p` — energy fraction (default 0.85)

Also exported as an `audio.js` stat manifest (`./audio` — mono-fold, 2048/512 hann frames, averaged over the take).

Part of [@audio/spectral](https://github.com/audiojs/spectral).
