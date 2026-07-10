# @audio/spectral-edit

> Spectral editing — gain on time×frequency regions via COLA-normalized STFT resynthesis (Audacity spectral edit; FFmpeg `afftfilt` class).

`npm install @audio/spectral-edit`

```js
import edit from '@audio/spectral-edit'

let out = edit(data, {
	fs: 44100,
	regions: [{ t0: 0.5, t1: 1.2, f0: 2000, f1: 8000, gain: 0 }],   // delete 2-8 kHz between 0.5-1.2 s
})
```

Options: - `fs` — sample rate (default 44100, Hz) · `frameSize` — STFT frame, samples (default 2048) · `hop` — STFT hop, samples (default `frameSize / 4`) · `regions` — array of `{ t0=0, t1=∞ (seconds), f0=0, f1=fs/2 (Hz), gain=0 (linear) }`

Part of [@audio/spectral](https://github.com/audiojs/spectral).
