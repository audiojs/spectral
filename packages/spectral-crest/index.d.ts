/** Spectral crest — peak-to-mean ratio of the power spectrum (≥1; high = peaky/tonal). Peeters 2004 §6.4. */
/** Magnitude-domain: pass an already-computed FFT magnitude spectrum. No options. */
export default function crest(mag: Float32Array | Float64Array): number
