/** Phase-vocoder engine — frame-level primitives over STFT magnitude/phase arrays. Pure math, no defaults. */

/** Per-frame geometry the scatter/lock primitives read from. */
export interface PvocCtx {
  /** Nyquist bin index (frameSize / 2) */
  half: number
  /** hop size, samples */
  hop: number
  /** Hz per FFT bin (fs / N) */
  freqPerBin: number
  /** analysis window / FFT size */
  N: number
}

export const PI2: number
/** rms(w)/mean(w) for the engine's periodic Hann — energy-preserving OLA gain for scatterGated. */
export const WIN_GAIN: number

/** Wrap a phase (radians) into (−π, π]. */
export function wrapPhase(p: number): number

/** Local magnitude maxima of `mag[0..half]` above a fraction of the frame's peak. Returns bin indices, ascending. */
export function findPeaks(mag: Float32Array | Float64Array, half: number): number[]
/** The same peaks into `out` (at least `half / 2` long: a peak's right neighbour is never one; shorter is a RangeError), their count returned: nothing allocated. */
export function findPeaks(mag: Float32Array | Float64Array, half: number, out: Int32Array): number

/** Binary-search the peak (from `findPeaks`) nearest bin `k`. Returns −1 if `peaks` is empty. */
export function nearestPeak(peaks: number[], k: number): number

/**
 * Bernsee/SMB peak-gated bin scatter — pitch shift by `ratio`, RMS-preserving collisions.
 * `newMag`/`newFreq`/`peakMag` are caller-owned scratch sized `half+1`, zero-filled before the call.
 * Writes into `newMag`/`newFreq` in place.
 * `ratio` may be a typed array holding it at [0], as hosts hand a live value (a Float64Array keeps every bit): a ratio that
 * is not whole, passed as a number each frame, is a number V8 allocates each frame (the function is too long to inline).
 */
export function scatterGated(
  mag: Float32Array | Float64Array,
  phase: Float32Array | Float64Array,
  prevPhase: Float32Array | Float64Array | null,
  ratio: number | Float64Array,
  ctx: PvocCtx,
  newMag: Float64Array,
  newFreq: Float64Array,
  peakMag: Float64Array,
): void

/**
 * Laroche-Dolson rigid-ROI peak-locked bin scatter — phase-coherent shift by `ratio`.
 * `syn` is the caller-owned running frame-center phase accumulator indexed by source bin (NaN means untracked), mutated in place.
 * `newMag`/`newPhase`/`peakMag` are caller-owned scratch sized `half+1`; `peakDest`/`peakSynPhase` sized `peaks.length`. All zero-filled before the call.
 * `count`: the peaks there are, the first ones of `peaks` (all of them when omitted), for peaks in `findPeaks`' `out`.
 * `ratio` may be a typed array holding it at [0], as `scatterGated`'s.
 */
export function scatterLocked(
  mag: Float32Array | Float64Array,
  phase: Float32Array | Float64Array,
  prevPhase: Float32Array | Float64Array | null,
  reset: boolean,
  peaks: ArrayLike<number>,
  ratio: number | Float64Array,
  ctx: PvocCtx,
  syn: Float64Array,
  newMag: Float64Array,
  newPhase: Float64Array,
  peakDest: Int32Array,
  peakSynPhase: Float64Array,
  peakMag: Float64Array,
  count?: number,
): void

/** Scalar or time-varying ratio → per-frame resolver. */
export function makeFrameRatio(ratio: number | ((t: number) => number)): {
  /** ratio at t = 0 */
  scalar: number
  /** ratio at a given frame position */
  at(frameStart: number, sampleRate: number): number
}

/** Lock non-peak bin phases to their nearest peak's rotation, in place on `propPhase` (time-stretch coherence). */
export function lockPhase(
  phase: Float32Array | Float64Array,
  propPhase: Float32Array | Float64Array,
  mag: Float32Array | Float64Array,
  half: number,
): void

/**
 * Which peak each bin locks to, into `owner` (sized ≥ `half+1`): the peak whose rotation bin k takes, or −1 for a
 * free bin (a peak itself, or under 3% of its peak's magnitude). The regions `lockPhase` uses; returns `owner`.
 */
export function lockMap(mag: Float32Array | Float64Array, half: number, owner: Int32Array): Int32Array

/** Complex-domain locked-vocoder state over `half+1` bins; the caller fills `mag` with each frame's magnitudes. */
export interface LockState {
  mag: Float64Array
  owner: Int32Array
  xr: Float64Array
  xi: Float64Array
  xm: Float64Array
  ur: Float64Array
  ui: Float64Array
}

/** State for `lockAdvance` over `half+1` bins. */
export function lockState(half: number): LockState

/**
 * One phase-locked vocoder frame on complex bins, in place: `re`/`im` hold the analysis bins, then the synthesis bins
 * (`st.mag` = this frame's magnitudes). `reset` restarts synthesis phases at the analysis phases (first frame,
 * transient); otherwise free bins advance at their instantaneous frequency and locked bins take their peak's rotation.
 */
export function lockAdvance(
  re: Float64Array, im: Float64Array, st: LockState, reset: boolean,
  anaHop: number, synHop: number, freqPerBin: number, half: number,
): void

/** `lockAdvance` with this frame's magnitudes, √(re² + im²), taken from its bins into `st.mag` first. */
export function lockStep(
  re: Float64Array, im: Float64Array, st: LockState, reset: boolean,
  anaHop: number, synHop: number, freqPerBin: number, half: number,
): void

/** State for `lockJoint` over `half+1` bins of `channels` channels. */
export interface JointState {
  st: LockState
  re: Float64Array
  im: Float64Array
  fr: Float64Array
  fi: Float64Array
  xr: Float64Array[]
  xi: Float64Array[]
}

export function jointState(half: number, channels: number): JointState

/**
 * One phase-locked vocoder frame of every channel, in place, one rotation for all: `res[c]`/`ims[c]` hold channel c's
 * analysis bins, then its synthesis bins. The step runs on a reference (the channels' root power, its phase advanced by
 * the angle of Σ X·conj(X′) over the channels), so the level and phase between channels stay the input's in every frame.
 */
export function lockJoint(
  res: Float64Array[], ims: Float64Array[], j: JointState, reset: boolean,
  anaHop: number, synHop: number, freqPerBin: number, half: number,
): void
