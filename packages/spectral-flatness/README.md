# @audio/spectral-flatness

> Spectral flatness — geometric/arithmetic mean ratio of the power spectrum, 0 (tonal) … 1 (noise). Peeters 2004 §6.3.

`npm install @audio/spectral-flatness`

```js
import { fft } from 'fourier-transform'
import flatness from '@audio/spectral-flatness'

let win = frame.map((x, i) => x * (0.5 - 0.5 * Math.cos(2 * Math.PI * i / (frame.length - 1))))
let [re, im] = fft(win)
let mag = re.map((r, k) => Math.hypot(r, im[k]))

let f = flatness(mag)   // 0 = tonal, 1 = noise-like
```

No options — `flatness(mag)` takes only the magnitude spectrum.

Part of [@audio/spectral](https://github.com/audiojs/spectral).
