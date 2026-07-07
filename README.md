# @audio/spectral

> Spectral features — centroid, spread, flatness, rolloff, flux, slope, crest. All planned.

Also planned: `spectral-mfcc` (librosa parity), `spectral-ltas` (adaptive/match-EQ substrate), `spectral-edit` (time×frequency region editing — Audacity/afftfilt class).

Reference implementations exist in the `audio` package stats (`centroid`, `flatness`, `crest` shipped there — extract, don't rewrite). Parity: FFmpeg `aspectralstats`, MIREX spectral features. LTAS / spectral-stats for the Auto-Chain analyser build on these atoms.
