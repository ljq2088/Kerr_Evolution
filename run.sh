#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
export OPENBLAS_NUM_THREADS=1
exec .venv/bin/python src/kerr_scalar.py "$@"
