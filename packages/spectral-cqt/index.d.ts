/** Constant-Q transform — log-frequency spectrum, one Q-scaled window per bin (Brown 1991). */
export interface CqtOptions {
  /** sample rate, default 44100 */
  fs?: number
  /** lowest bin center, default 55 Hz (A1) */
  fmin?: number
  /** bins per octave, default 12 (semitone resolution) */
  binsPerOctave?: number
  /** octave span, default 6 */
  octaves?: number
  /** analysis center, seconds, default 0 */
  at?: number
}

export interface CqtResult {
  /** bin center frequencies, Hz, log-spaced, length binsPerOctave * octaves */
  freqs: Float32Array
  /** magnitude at each bin, same length as freqs */
  mag: Float32Array
}

/** Signal-domain: pass raw PCM, computes its own per-bin windowed DFT. One analysis frame per call. */
export default function cqt(data: Float32Array | Float64Array, options?: CqtOptions): CqtResult
