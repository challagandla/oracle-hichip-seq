#!/usr/bin/env bash
# Run the repository tests in the same activation-free environment as the workflow.
set -Eeuo pipefail

ENV_NAME="oracle-hichip-runner"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

die() { printf '[ERROR] %s\n' "$1" >&2; exit 1; }

CONDA_BIN=""
for candidate in \
    "${CONDA_EXE:-}" \
    "$(type -P conda 2>/dev/null || true)" \
    "$HOME/miniforge3/bin/conda" \
    "$HOME/miniconda3/bin/conda"; do
    if [[ -n "$candidate" && -x "$candidate" ]]; then
        CONDA_BIN="$candidate"
        break
    fi
done

[[ -n "$CONDA_BIN" ]] || die "Conda is missing; run: bash setup.sh"
"$CONDA_BIN" run --name "$ENV_NAME" python -c \
    'import cooler, h5py, numpy, pandas, pytest, yaml' >/dev/null 2>&1 || \
    die "Runner $ENV_NAME is missing test dependencies; run: bash setup.sh"

exec env -u R_LIBS -u R_LIBS_USER -u PYTHONPATH \
    R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null \
    "$CONDA_BIN" run --no-capture-output --name "$ENV_NAME" \
    pytest -q tests
