/** LTAS — long-term average spectrum via Welch's method: hann-windowed frames, power-averaged, returned as RMS magnitude per bin. */
export interface LtasOptions {
  /** Welch analysis window, samples, default 4096 */
  frameSize?: number
  /** frame hop, samples, default frameSize / 2 (50% overlap) */
  hop?: number
}

/** Signal-domain: pass raw PCM, computes its own framed FFT. Returns Float32Array, length frameSize/2 + 1. */
export default function ltas(data: Float32Array | Float64Array, options?: LtasOptions): Float32Array
