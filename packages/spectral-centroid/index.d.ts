/** Spectral centroid — magnitude-weighted mean frequency (Hz). Peeters 2004 §6.1. */
export interface CentroidOptions {
  /** sample rate, default 44100 */
  fs?: number
  /** FFT size the spectrum was taken at, default 2 * (mag.length - 1) */
  n?: number
}

/** Magnitude-domain: pass an already-computed FFT magnitude spectrum. Returns the centroid, Hz. */
export default function centroid(mag: Float32Array | Float64Array, options?: CentroidOptions): number
