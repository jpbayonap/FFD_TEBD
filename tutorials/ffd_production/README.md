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

On macOS the Makefile links against Accelerate. On Linux it first looks for a
system OpenBLAS shared library at common cluster paths such as
`/usr/lib64/libopenblas.so.0`. If that is present, a plain build should work:

```bash
make -C tutorials/ffd_production print-config
make ffd-production CXX=g++
```

If the cluster provides normal linker names, you can override the link flags:

```bash
make ffd-production CXX=g++ LDFLAGS="-llapack -lblas"
```

If only a versioned OpenBLAS library is visible, use the full path:

```bash
make ffd-production CXX=g++ LDFLAGS="/usr/lib64/libopenblas.so.0"
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
mkdir -p logs
sbatch tutorials/ffd_production/cluster/slurm_entropy_scan.sbatch
```

On the observed cluster, the matrix partitions are the best starting point for
CPU TEBD. Use `m1` for first production scans: it has 56 CPUs and about 257 GB
per node, which is enough for several moderate independent jobs without
occupying the largest memory nodes.

## Threading Notes

Start with one independent TEBD run per job-array task. On the observed matrix
partition, a short `N=40`, `chi=96`, second-order benchmark gave nearly flat
runtime for 1, 2, 4, and 8 OpenBLAS threads, with 2 threads slightly fastest.
Use 2 threads as the current default:

```bash
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2
```

Do not combine many job-array tasks per node with many BLAS threads until
`bench_node.sh` shows that this is actually faster on the cluster.

## Useful Accounting

```bash
sacct -j JOBID --format=JobID,JobName,Partition,State,Elapsed,MaxRSS,AllocCPUS,ReqMem
```
