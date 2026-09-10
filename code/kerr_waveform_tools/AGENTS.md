# Kerr Waveform Tools Guidance

- Keep this package independent of `sminus2_point_particle`. Shared APIs should expose units, conventions, background parameters, and provenance when they affect interpretation.
- Add an abstraction when an implemented object has a real cross-task use or removes meaningful duplication. Do not prebuild QNM-filter, radial-cache, or SPA interfaces before their tasks establish the required behavior.
- Solver-specific grids, boundary rows, sources, evolution state, and checkpoints normally remain in the solver package.
- Preserve compatibility re-exports only while they help existing reviewed artifacts; avoid maintaining two implementations of the same physical formula.
