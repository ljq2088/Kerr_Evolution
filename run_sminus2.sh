#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
export OPENBLAS_NUM_THREADS=1
exec .venv/bin/python code/kerr_point_particle_solver/scripts/run_case.py \
  code/kerr_point_particle_solver/configs/smoke/evolution.json \
  code/kerr_point_particle_solver/configs/smoke/output.json \
  code/kerr_point_particle_solver/tests/fixtures/chi0p8_l6_flux.json "$@"
