# @audio/spectral-pvoc

> Phase-vocoder engine — the spectral scatter/lock machinery behind pitch shifting and time stretching.

Frame-level primitives over STFT magnitude/phase arrays. Pure math, no dependencies; pair with [`fourier-transform/stft`](https://github.com/scijs/fourier-transform) (ana/syn hop) or [`@audio/stft`](https://github.com/audiojs/stft) (AMS).

```js
import { findPeaks, scatterLocked, lockPhase, makeFrameRatio } from '@audio/spectral-pvoc'
```

| export | role |
|---|---|
| `scatterGated(mag, phase, prevPhase, ratio, ctx, …scratch)` | Bernsee/SMB peak-gated bin scatter — pitch shift by ratio, RMS-preserving collisions |
| `scatterLocked(mag, phase, prevPhase, reset, peaks, ratio, ctx, …scratch, count?)` | Laroche-Dolson rigid-ROI scatter — phase-coherent shift; `count`: the peaks in `findPeaks`' `out` |
| `lockPhase(phase, propPhase, mag, half)` | lock non-peak bins to their peak's rotation (time-stretch coherence) |
| `lockMap(mag, half, owner)` | which peak each bin locks to (−1: free): the regions `lockPhase` uses |
| `lockState(half)` / `lockAdvance(re, im, st, reset, anaHop, synHop, freqPerBin, half)` | the phase-locked vocoder frame on complex bins: locked bins rotate by one complex multiply, only free bins take trig; with whole hops, a WebAssembly SIMD kernel on the state's own memory, the JS loop's bits |
| `lockStep(re, im, st, reset, anaHop, synHop, freqPerBin, half)` | `lockAdvance` with the frame's magnitudes taken from its bins first |
| `jointState(half, channels)` / `lockJoint(res, ims, j, reset, anaHop, synHop, freqPerBin, half)` | that frame for several channels at once, one rotation for all: the step runs on a reference (the channels' root power, its phase advanced by the angle of Σ X·conj(X′) over them) and its rotation turns every channel's bin, so the level and phase between channels stay the input's; the kernel's passes with whole hops, the JS loops' bits |
| `findPeaks(mag, half, out?)` / `nearestPeak(peaks, k)` | local-maxima peak indices + nearest lookup; given `out` (an `Int32Array` of `half / 2`), the peaks go into it and their count comes back |
| `makeFrameRatio(ratio)` | scalar or time-varying ratio → per-frame resolver |
| `wrapPhase(p)` / `WIN_GAIN` / `PI2` | phase wrap, Hann OLA gain constant |

A frame of a real-time loop allocates nothing: peaks into `findPeaks`' `out`, and the ratio of `scatterGated` and `scatterLocked` in a `Float64Array` holding it at [0], as hosts hand a live value. A ratio that is not whole, passed as a number, is a number V8 allocates at each call: the two are too long to inline.

Used by [`@audio/shift`](https://github.com/audiojs/shift)'s pvoc/pvoc-lock/sms/formant/transient/hpss and [`@audio/stretch`](https://github.com/audiojs/stretch)'s pvoc-lock/transient. MIT.
