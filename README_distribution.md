# Distribution Record

Prepared on 2026-09-10 for owner review. This record covers organization and reproducibility checks, not new scientific validation.

## File Organization

The following existing trees were moved here with their relative hierarchy intact. Their parent-repository entry paths now contain relative symbolic links:

- `code/kerr_point_particle_solver/`
- `code/kerr_waveform_tools/`
- `code/kerr_point_particle_verification/`
- `code/kerr_point_particle_solver_migration.md`
- `findings/kerr_point_particle_evolution/`
- `data/kerr_point_particle_verification/`
- `data/kerr_point_particle_evolution/sxs0305_medium/`
- `workspace/plan/W02_point_particle_evolution/`

The shared `REFERENCE.md`, its two `resource_details/` documents, and the data-directory README were copied unchanged to retain their original context. They remain historical snapshots in this distribution. No Python or WL interface edits were needed. This distribution contains actual files, with no symlinks pointing outside it.

Generated Python/pytest caches were separated into `/tmp/particle-publication-caches/`; they are not publication material. The original files were checked against a pre-move SHA256 inventory: all 153 non-cache files matched, including 40 Python files and 26 WL files (25 symbolic artifacts plus the flux provider). Only new README files were written.

## Executed Checks

- Existing solver and shared-tool tests: **63 passed**, executed from the solver package directory with both source directories on `PYTHONPATH`.
- Reference-data manifest: **16 files passed** the existing SHA256 validator.
- Existing spin-batch dependency manifest: all **32 hashes matched**, and the generated batch specification equalled its saved plan after relocation.
- No heavy Mathematica checks or full production waveforms were recomputed for this organization task.

Environment: Python 3.13.1, NumPy 2.2.5, SciPy 1.16.0, Matplotlib 3.10.0, pytest 8.4.0. The first test invocation from the distribution root produced 61 passes and two failures because subprocess tests expect `scripts/` relative to the solver package. Running from that package resolved both failures without source changes. Follow the top-level README commands; the byte-preserved older solver README has a root-level test invocation with this working-directory limitation.

## Archival Link Limitations

The current theory, code, symbolic-result, and reference-data navigation is included. A file-target scan identified these five unresolved references in preserved historical/supporting documents:

| Document | Historical target |
| --- | --- |
| W02 main plan | `workspace/W02_point_particle_evolution/README.md` (already absent in the parent repository) |
| W02 main plan | `findings/W01_sxs_waveform_inputs.md` (separate W01 study) |
| Solver migration record | `findings/kerr_point_particle_evolution/migration_plan.md` (already absent in the parent repository) |
| Local-resource index | `workspace/plan/W08_end_to_end_error/task3_plan.md` (later study) |
| Paper index | `workspace/W08_lu_spin_scan/interpretation.md` (later study) |

Original absolute local paths and historical Markdown anchors also remain as provenance; external resources and all section anchors were not validated. The independent cross-code waveform referenced by an absolute local path is not bundled. These archival references are not used by the Python example or the bundled symbolic-summary interfaces. They were not silently rewritten because this task permits content edits only in README files and script path interfaces.

## Publication and Review

The directory can be used as a standalone GitHub repository root. No Git repository was initialized, staged, committed, or published. No license or author/citation metadata was invented. Review the retained local-path provenance and choose licensing/citation metadata before public release.

Parent PROJECT.md, TASK.md, work_log.md, and AGENTS/FutureSKILL rules were inspected but not edited under this task's README-only content boundary. The parent README supplies the new navigation entry while the compatibility links preserve existing paths. The original T006--T007 identifiers describe scientific provenance; no new owner task identifier was assigned here. Please commit the organization after review.
