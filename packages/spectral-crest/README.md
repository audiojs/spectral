# @audio/spectral-crest

> Spectral crest — peak-to-mean ratio of the power spectrum (≥1; high = peaky/tonal). FFmpeg `aspectralstats` / Peeters 2004 §6.4.

`npm install @audio/spectral-crest`

```js
import { fft } from 'fourier-transform'
import crest from '@audio/spectral-crest'

let win = frame.map((x, i) => x * (0.5 - 0.5 * Math.cos(2 * Math.PI * i / (frame.length - 1))))
let [re, im] = fft(win)
let mag = re.map((r, k) => Math.hypot(r, im[k]))

let c = crest(mag)   // 1 = flat spectrum, higher = one dominant bin
```

No options — `crest(mag)` takes only the magnitude spectrum.

Part of [@audio/spectral](https://github.com/audiojs/spectral).
