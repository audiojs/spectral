# @audio/spectral-harmonics

> Harmonic descriptors from f0 — partial amplitudes at k·f0 give inharmonicity, tristimulus (T1/T2/T3) and odd/even energy ratio (essentia parity).

`npm install @audio/spectral-harmonics`

```js
import harmonics from '@audio/spectral-harmonics'

let h = harmonics(frame, { fs: 44100, f0: 220 })
// { amps: Float32Array(10), inharmonicity: number, tristimulus: [n, n, n], oddEven: number }
```

Options: - `fs` — sample rate (default 44100, Hz) · `f0` — fundamental frequency, Hz (required, throws `RangeError` if omitted) · `nHarmonics` — partials tracked (default 10)

Part of [@audio/spectral](https://github.com/audiojs/spectral).
