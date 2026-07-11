/** Spectral flux — frame-to-frame magnitude change. Euclidean by default (FFmpeg aspectralstats). */
export interface FluxOptions {
  /** sum only positive frame-to-frame increases instead of Euclidean distance, default false */
  rectified?: boolean
}

/** Magnitude-domain: current frame's magnitude spectrum plus the previous frame's (or null/undefined on the first frame, returns 0). */
export default function flux(mag: Float32Array | Float64Array, prev?: Float32Array | Float64Array | null, options?: FluxOptions): number
