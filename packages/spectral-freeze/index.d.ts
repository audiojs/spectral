/** Spectral freeze — sustain one spectral frame with per-bin phase dispersion (Ableton Spectral Time / Paulstretch-freeze class). */
export interface FreezeOptions {
  /** sample rate, default 44100 */
  fs?: number
  /** frame to freeze, seconds, default 0 */
  at?: number
  /** output length, seconds, default 2 */
  duration?: number
  /** STFT frame, samples, default 4096 */
  frameSize?: number
  /** phase-dispersion RNG seed, default 1 */
  seed?: number
}

/** Signal-domain effect: pass raw PCM, returns new audio of length round(duration * fs). */
export default function freeze(data: Float32Array | Float64Array, options?: FreezeOptions): Float32Array
