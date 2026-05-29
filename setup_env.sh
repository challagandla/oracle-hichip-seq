#!/usr/bin/env bash
# ORACLE HiChIP — one-command environment setup
# Run from code/hichip/ directory.
set -euo pipefail

ENV_NAME="${ENV_NAME:-oracle-hichip}"

echo "==> Building conda environment '${ENV_NAME}' with mamba ..."
if command -v mamba >/dev/null 2>&1; then
    SOLVER=mamba
elif command -v micromamba >/dev/null 2>&1; then
    SOLVER=micromamba
else
    echo "Installing micromamba (mamba/conda not found)..."
    curl -L micro.mamba.pm/install.sh | bash
    SOLVER=micromamba
fi

$SOLVER env create -y -f environment.yml -n "${ENV_NAME}" || \
    $SOLVER env update -y -f environment.yml -n "${ENV_NAME}"

# Lock the env for reproducibility
echo "==> Locking environment to environment.lock.yml ..."
$SOLVER env export -n "${ENV_NAME}" --no-builds > environment.lock.yml

# Activate hint
cat <<'EOF'

==> Done.

Activate with:
    mamba activate oracle-hichip
or
    micromamba activate oracle-hichip

Verify install:
    snakemake --version
    cooler --version
    pairtools --version
    macs2 --version
    fithichip --help | head -5

Run pipeline:
    snakemake -n --configfile config/config.yaml   # dry run / DAG
    snakemake --cores 32 --configfile config/config.yaml

EOF
