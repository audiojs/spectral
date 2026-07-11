/** Spectral rolloff — frequency below which `p` of total spectral energy lies. */
export interface RolloffOptions {
  /** sample rate, default 44100 */
  fs?: number
  /** FFT size the spectrum was taken at, default 2 * (mag.length - 1) */
  n?: number
  /** energy fraction, default 0.85 */
  p?: number
}

/** Magnitude-domain: pass an already-computed FFT magnitude spectrum. Returns the rolloff frequency, Hz. */
export default function rolloff(mag: Float32Array | Float64Array, options?: RolloffOptions): number
