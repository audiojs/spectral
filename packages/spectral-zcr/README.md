# @audio/spectral-zcr

> Zero-crossing rate — fraction of consecutive samples whose sign differs (librosa `zero_crossing_rate` class).

`npm install @audio/spectral-zcr`

```js
import zcr from '@audio/spectral-zcr'

let r = zcr(frame)   // 0…1, sign-change fraction over the frame
```

Time-domain, not spectral — no FFT, takes raw PCM directly. `−0` counts as negative (matches
librosa's `np.signbit`), so a hard-zero sample still registers a crossing against a preceding
positive sample. High ZCR marks noisy/fricative/percussive content; low ZCR marks tonal/voiced
content.

No options — `zcr(data)` takes only the signal.

Also exported as an `audio.js` stat manifest (`./audio` — mono-fold, 2048/512 frames, unwindowed, averaged over the take).

Part of [@audio/spectral](https://github.com/audiojs/spectral).
