/** Harmonic descriptors from f0 — partial amplitudes at k·f0 give inharmonicity, tristimulus, odd/even ratio (essentia parity). */
export interface HarmonicsOptions {
  /** sample rate, default 44100 */
  fs?: number
  /** fundamental frequency, Hz — required, throws RangeError if omitted */
  f0: number
  /** partials tracked, default 10 */
  nHarmonics?: number
}

export interface HarmonicsResult {
  /** per-partial peak amplitude, length nHarmonics */
  amps: Float32Array
  /** energy-weighted mean deviation of each partial from an exact integer multiple of f0 */
  inharmonicity: number
  /** [T1, T2, T3] — fundamental / partials 2-4 / partials 5+ energy fractions */
  tristimulus: [number, number, number]
  /** odd-harmonic energy / even-harmonic energy */
  oddEven: number
}

/** Signal-domain: pass a raw PCM frame (power-of-2 length), computes its own FFT. */
export default function harmonics(data: Float32Array | Float64Array, options: HarmonicsOptions): HarmonicsResult
