#!/usr/bin/env bash
#SBATCH --job-name=ffd_entropy
#SBATCH --output=logs/ffd_entropy_%A_%a.out
#SBATCH --error=logs/ffd_entropy_%A_%a.err
#SBATCH --partition=m1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --array=0-5

set -euo pipefail

# Edit these two paths for your cluster checkout/scratch layout.
REPO_DIR="${REPO_DIR:-$PWD}"
OUT_DIR="${OUT_DIR:-$REPO_DIR/tutorials/ffd_production/data}"
BIN="$REPO_DIR/tutorials/ffd_production/build/bin/ffd_production"

mkdir -p "$OUT_DIR" "$REPO_DIR/logs"

cd "$REPO_DIR"

# Load modules here if your cluster uses environment modules, for example:
# module purge
# module load astro/gcc/10.3.0 astro/blas/3.11.0 astro/lapack/3.11

make -C tutorials/ffd_production print-config
make ffd-production CXX="${CXX:-g++}"

# Avoid oversubscription. Match these to --cpus-per-task unless benchmarking says otherwise.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export VECLIB_MAXIMUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

CHI_VALUES=(64 128 192 256 128 192)
DT_VALUES=(0.1 0.1 0.1 0.1 0.05 0.05)
STEP_VALUES=(50 50 50 50 100 100)
ORDER_VALUES=(2 2 2 2 2 2)

IDX="${SLURM_ARRAY_TASK_ID:-0}"
CHI="${CHI_VALUES[$IDX]}"
DT="${DT_VALUES[$IDX]}"
STEPS="${STEP_VALUES[$IDX]}"
ORDER="${ORDER_VALUES[$IDX]}"

TAG="slurm_${SLURM_JOB_ID:-local}_${IDX}"

RUN_CMD=()
if command -v /usr/bin/time >/dev/null 2>&1; then
  RUN_CMD=(/usr/bin/time -v)
elif command -v gtime >/dev/null 2>&1; then
  RUN_CMD=(gtime -v)
fi

"${RUN_CMD[@]}" "$BIN" \
  --mode entropy \
  --N 40 \
  --dt "$DT" \
  --steps "$STEPS" \
  --chi "$CHI" \
  --cutoff 1e-10 \
  --order "$ORDER" \
  --ells 5,10,15,20 \
  --entropy-measure-every 5 \
  --progress-every 5 \
  --out-dir "$OUT_DIR" \
  --tag "$TAG"
