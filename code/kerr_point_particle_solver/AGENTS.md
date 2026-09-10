# Kerr Point-Particle Solver Guidance

- Treat the equations in `findings/kerr_point_particle_evolution/` as the scientific reference; do not maintain a second editable derivation here.
- Keep fixed-$m$ field/source, spatial closure, RK evolution, checkpointing, and run orchestration in this package. Prefer `kerr_waveform_tools` for shared Kerr, worldline, waveform, and strain objects.
- Repository production runs should normally use the direct-to-data staging workflow. Explicit output paths remain useful for tests and external reuse.
- Changes that may alter the waveform should usually be developed and reviewed in workspace first. Schema and I/O changes should preserve useful historical readers when practical.
