/** Spectral editing — gain on a time×frequency region (delete / attenuate / boost). */
export interface SpectralEditRegion {
  /** region start, seconds, default 0 */
  t0?: number
  /** region end, seconds, default Infinity */
  t1?: number
  /** region low edge, Hz, default 0 */
  f0?: number
  /** region high edge, Hz, default fs/2 */
  f1?: number
  /** linear gain applied to the region, default 0 (delete) */
  gain?: number
}

export interface EditOptions {
  /** sample rate, default 44100 */
  fs?: number
  /** STFT frame, samples, default 2048 */
  frameSize?: number
  /** STFT hop, samples, default frameSize / 4 */
  hop?: number
  /** gain regions, applied where they overlap a frame's center time and bin range */
  regions?: SpectralEditRegion[]
}

/** Signal-domain effect: pass raw PCM, returns a new edited buffer of the same length. */
export default function edit(data: Float32Array | Float64Array, options?: EditOptions): Float32Array
