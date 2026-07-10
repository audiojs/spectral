# @audio/spectral-centroid

> Spectral centroid — magnitude-weighted mean frequency (Hz). Peeters 2004 §6.1.

`npm install @audio/spectral-centroid`

```js
import { fft } from 'fourier-transform'
import centroid from '@audio/spectral-centroid'

let win = frame.map((x, i) => x * (0.5 - 0.5 * Math.cos(2 * Math.PI * i / (frame.length - 1))))
let [re, im] = fft(win)                          // frame.length a power of 2
let mag = re.map((r, k) => Math.hypot(r, im[k])) // magnitude spectrum, length frame.length/2+1

let hz = centroid(mag, { fs: 44100 })
```

Options: - `fs` — sample rate (default 44100, Hz) · `n` — FFT size the spectrum was taken at (default `2 * (mag.length - 1)`, inferred from `mag.length`)

Part of [@audio/spectral](https://github.com/audiojs/spectral).
