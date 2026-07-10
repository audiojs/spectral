# @audio/spectral-flux

> Spectral flux — frame-to-frame magnitude change. Euclidean by default (FFmpeg `aspectralstats`); `{ rectified: true }` sums positive differences only (onset-detection flavour).

`npm install @audio/spectral-flux`

```js
import { fft } from 'fourier-transform'
import flux from '@audio/spectral-flux'

let prev = null
for (const frame of frames) {              // consecutive, equal-length frames
	let win = frame.map((x, i) => x * (0.5 - 0.5 * Math.cos(2 * Math.PI * i / (frame.length - 1))))
	let [re, im] = fft(win)
	let mag = re.map((r, k) => Math.hypot(r, im[k]))
	let v = flux(mag, prev)                  // 0 on the first frame (no prev)
	prev = mag
}
```

Options: - `rectified` — sum only positive frame-to-frame increases instead of Euclidean distance (default `false`)

Also exported as an `audio.js` stat manifest (`./audio` — mono-fold, 2048/512 hann frames, averaged over the take).

Part of [@audio/spectral](https://github.com/audiojs/spectral).
