/** Spectral flatness — geometric/arithmetic mean ratio of the power spectrum, 0 (tonal) … 1 (noise). Peeters 2004 §6.3. */
/** Magnitude-domain: pass an already-computed FFT magnitude spectrum. No options. */
export default function flatness(mag: Float32Array | Float64Array): number
