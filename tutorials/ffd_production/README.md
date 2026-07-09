# FFD Production TEBD Runner

This executable is the cluster-oriented version of the manual FFD TEBD tutorial.
The tutorial code remains in `tutorials/ffd_manual`; this directory is for
parameter scans, benchmarking, and long-running entropy/observable production
jobs.

## Build

From the repository root:

```bash
make ffd-production
```

On macOS the Makefile links against Accelerate. On Linux it links against
LAPACK/BLAS:

```bash
make ffd-production CXX=g++ LDFLAGS="-llapack -lblas"
```

If your cluster provides MKL, override the link flags according to the module
documentation, for example:

```bash
make ffd-production CXX=g++ LDFLAGS="$MKL_LINK_FLAGS"
```

## Smoke Test

```bash
./tutorials/ffd_production/build/bin/ffd_production \
  --mode entropy \
  --N 40 \
  --dt 0.1 \
  --steps 5 \
  --chi 32 \
  --order 2 \
  --ell-count 5 \
  --entropy-measure-every 5 \
  --tag smoke
```

Output goes to:

```text
tutorials/ffd_production/data/
```

## Entropy Scan Example

```bash
./tutorials/ffd_production/build/bin/ffd_production \
  --mode entropy \
  --N 40 \
  --dt 0.1 \
  --steps 50 \
  --chi 128 \
  --order 2 \
  --ells 5,10,15,20 \
  --entropy-measure-every 5 \
  --progress-every 5 \
  --tag scanA
```

## Observable Scan Example

```bash
./tutorials/ffd_production/build/bin/ffd_production \
  --mode hm \
  --N 40 \
  --dt 0.1 \
  --steps 50 \
  --chi 128 \
  --order 2 \
  --probes 18,19,20,21,22 \
  --hm-measure-every 5 \
  --progress-every 5 \
  --tag hm_scanA
```

## Cluster Benchmarking

On an interactive compute node:

```bash
bash tutorials/ffd_production/cluster/bench_node.sh
```

For SLURM job arrays, edit the partition/account/time/memory directives in:

```text
tutorials/ffd_production/cluster/slurm_entropy_scan.sbatch
```

Then submit:

```bash
sbatch tutorials/ffd_production/cluster/slurm_entropy_scan.sbatch
```

## Threading Notes

Start with one independent TEBD run per job-array task. For each task, use a
small number of BLAS/LAPACK threads:

```bash
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
```

Do not combine many job-array tasks per node with many BLAS threads until
`bench_node.sh` shows that this is actually faster on the cluster.
