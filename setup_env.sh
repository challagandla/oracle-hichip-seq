#!/usr/bin/env bash
# ORACLE HiChIP — one-command environment setup
set -euo pipefail

ENV_NAME="${ENV_NAME:-oracle-hichip}"

echo "==> Building conda environment '${ENV_NAME}' with mamba/micromamba ..."
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

echo "==> Locking environment to environment.lock.yml ..."
$SOLVER env export -n "${ENV_NAME}" --no-builds > environment.lock.yml

cat <<'EOF'

==> Done.

Activate with:
    mamba activate oracle-hichip
or:
    micromamba activate oracle-hichip

Verify install:
    snakemake --version
    cooler --version
    pairtools --version
    macs2 --version
    FitHiChIP_HiCPro.sh -h | head -5 || fithichip --help | head -5

Run from repository root:
    snakemake -s workflow/Snakefile -n --configfile config/config.yaml
    snakemake -s workflow/Snakefile --cores 32 --configfile config/config.yaml --use-conda

EOF
