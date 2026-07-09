#!/usr/bin/env bash
set -euo pipefail

BIN="${1:-./tutorials/ffd_production/build/bin/ffd_production}"
OUT_DIR="${2:-tutorials/ffd_production/data}"

echo "== scheduler tools =="
command -v sbatch || true
command -v srun || true
command -v sinfo || true
command -v qsub || true
command -v qstat || true

echo
echo "== host =="
hostname
date

echo
echo "== CPU and memory =="
lscpu || true
nproc --all || true
free -h || true
numactl --hardware || true
ulimit -a || true

echo
echo "== compiler and libraries =="
which c++ || true
c++ --version || true
which g++ || true
g++ --version || true
ldd "$BIN" || true

echo
echo "== threading environment =="
echo "OMP_NUM_THREADS=${OMP_NUM_THREADS:-}"
echo "OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-}"
echo "MKL_NUM_THREADS=${MKL_NUM_THREADS:-}"
echo "VECLIB_MAXIMUM_THREADS=${VECLIB_MAXIMUM_THREADS:-}"

echo
echo "== filesystem =="
df -h . || true
df -h "${TMPDIR:-/tmp}" || true
mkdir -p "$OUT_DIR"

RUN_CMD=()
if command -v /usr/bin/time >/dev/null 2>&1; then
  RUN_CMD=(/usr/bin/time -v)
elif command -v gtime >/dev/null 2>&1; then
  RUN_CMD=(gtime -v)
fi

echo
echo "== smoke benchmark =="
"${RUN_CMD[@]}" "$BIN" \
  --mode entropy \
  --N 40 \
  --dt 0.1 \
  --steps 5 \
  --chi 32 \
  --order 2 \
  --ell-count 5 \
  --entropy-measure-every 5 \
  --out-dir "$OUT_DIR" \
  --tag "bench_smoke_$(date +%Y%m%d_%H%M%S)"

echo
echo "== thread scaling benchmark =="
for threads in 1 2 4 8 16; do
  export OMP_NUM_THREADS="$threads"
  export OPENBLAS_NUM_THREADS="$threads"
  export MKL_NUM_THREADS="$threads"
  export VECLIB_MAXIMUM_THREADS="$threads"
  echo "-- threads=$threads"
  "${RUN_CMD[@]}" "$BIN" \
    --mode entropy \
    --N 40 \
    --dt 0.1 \
    --steps 10 \
    --chi 64 \
    --order 2 \
    --ell-count 5 \
    --entropy-measure-every 5 \
    --out-dir "$OUT_DIR" \
    --tag "bench_threads${threads}_$(date +%Y%m%d_%H%M%S)"
done
