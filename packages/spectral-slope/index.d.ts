/** Spectral slope — least-squares linear regression of magnitude over frequency. Peeters 2004 §6.6. */
export interface SlopeOptions {
  /** sample rate, default 44100 */
  fs?: number
  /** FFT size the spectrum was taken at, default 2 * (mag.length - 1) */
  n?: number
}

/** Magnitude-domain: pass an already-computed FFT magnitude spectrum. Returns the regression coefficient, magnitude/Hz (negative = downward tilt). */
export default function slope(mag: Float32Array | Float64Array, options?: SlopeOptions): number
