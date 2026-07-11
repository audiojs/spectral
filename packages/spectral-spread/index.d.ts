/** Spectral spread — magnitude-weighted standard deviation around the centroid (Hz). Peeters 2004 §6.1. */
export interface SpreadOptions {
  /** sample rate, default 44100 */
  fs?: number
  /** FFT size the spectrum was taken at, default 2 * (mag.length - 1) */
  n?: number
}

/** Magnitude-domain: pass an already-computed FFT magnitude spectrum. Returns the spread, Hz. */
export default function spread(mag: Float32Array | Float64Array, options?: SpreadOptions): number
