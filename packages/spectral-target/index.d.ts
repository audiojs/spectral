/** Target-curve library — per-content-type LTAS targets + the deviation math for adaptive/match EQ. Curve data + pure math, no DSP kernel. */
export type TargetPreset = 'speech' | 'music' | 'pink' | 'voice-music'
/** Custom anchor table, ascending frequency: [[freqHz, dB], ...] */
export type TargetAnchors = Array<[number, number]>

export interface TargetOptions {
  /** sample rate, default 44100 */
  fs?: number
  /** output length (bins), default 2049 */
  bins?: number
}

/**
 * Render a preset or custom anchor table to a dB-per-bin curve, mean-normalized to 0 dB over
 * its own defined band = [lowest anchor, min(highest anchor, 0.45·fs)].
 */
export default function target(name: TargetPreset | TargetAnchors, options?: TargetOptions): Float32Array

export interface SmoothOptions {
  /** sample rate, default 44100 */
  fs?: number
  /** smoothing width, octaves, default 1/3 */
  oct?: number
}

/** Fractional-octave (log-frequency) smoothing of a per-bin dB curve. */
export function smooth(db: Float32Array, options?: SmoothOptions): Float32Array

export interface DeviationOptions {
  /** sample rate, default 44100 */
  fs?: number
  /** smoothing width, octaves, default 1/3 */
  smoothOct?: number
  /** correction clamp, dB, default 12 */
  clamp?: number
}

/**
 * (measured LTAS, target) → EQ correction curve, dB per bin (target − measured, smoothed,
 * clamped, tapered to exactly 0 outside [20 Hz, 0.45·fs]).
 */
export function deviation(measuredLtas: Float32Array, targetDb: Float32Array, options?: DeviationOptions): Float32Array
