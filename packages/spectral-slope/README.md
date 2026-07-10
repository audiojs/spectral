# @audio/spectral-slope

> Spectral slope — least-squares linear regression of magnitude over frequency. Peeters 2004 §6.6.

`npm install @audio/spectral-slope`

```js
import { fft } from 'fourier-transform'
import slope from '@audio/spectral-slope'

let win = frame.map((x, i) => x * (0.5 - 0.5 * Math.cos(2 * Math.PI * i / (frame.length - 1))))
let [re, im] = fft(win)
let mag = re.map((r, k) => Math.hypot(r, im[k]))

let m = slope(mag, { fs: 44100 })   // magnitude/Hz regression coefficient; negative = downward tilt
```

Options: - `fs` — sample rate (default 44100, Hz) · `n` — FFT size the spectrum was taken at (default `2 * (mag.length - 1)`, inferred from `mag.length`)

Also exported as an `audio.js` stat manifest (`./audio` — mono-fold, 2048/512 hann frames, averaged over the take).

Part of [@audio/spectral](https://github.com/audiojs/spectral).
