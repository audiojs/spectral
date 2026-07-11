/** Spectral contrast — per-octave-band peak/valley energy difference. Jiang 2002, librosa parity. */
export interface ContrastOptions {
  /** sample rate, default 44100 */
  fs?: number
  /** FFT size the spectrum was taken at, default 2 * (mag.length - 1) */
  n?: number
  /** lowest band edge, default 200 Hz */
  fmin?: number
  /** number of octave bands, default 6 */
  bands?: number
  /** peak/valley fraction of each band's bins, default 0.2 */
  quantile?: number
}

/** Magnitude-domain: pass an already-computed FFT magnitude spectrum. Returns dB per band. */
export default function contrast(mag: Float32Array | Float64Array, options?: ContrastOptions): Float32Array
