# FFD_GGE

C++ tensor-network project built around the ITensor C++ library.

## Project strategy

We will work in three layers:

1. Build and verify ITensor C++ locally.
2. Reproduce standard MPO calculations in two passes:
   - library-assisted constructions first;
   - manual constructions second.
3. Turn the folder into a clean standalone GitHub project.

The first concrete track will be a standard Ising-chain MPO example. After that, we can add a Hubbard track using the same pattern.

## Layout

- `docs/` setup notes and tutorial plans
- `tutorials/ising_builtin/` library-driven reproduction
- `tutorials/ising_manual/` manual MPO reconstruction
- `vendor/` external dependencies kept outside the main source tree

## First setup commands

Run these yourself from the project root:

```bash
cd /Users/juan/Desktop/Git/FFD_GGE
mkdir -p vendor
git clone --branch v3.2.0 --depth 1 https://github.com/ITensor/ITensor.git vendor/ITensor
cd vendor/ITensor
cp options.mk.sample options.mk
```

Then edit `vendor/ITensor/options.mk` on macOS:

- comment the default `g++` `CCCOM=...` line
- uncomment the `clang++` line
- keep `PLATFORM=macos`
- keep `BLAS_LAPACK_LIBFLAGS=-framework Accelerate`

Then build:

```bash
cd /Users/juan/Desktop/Git/FFD_GGE/vendor/ITensor
make -j1
```

`-j1` is intentional because the official install notes warn that parallel builds can be unreliable for ITensor.

Once that finishes, the project helper targets in the root `Makefile` will work.

## Next milestone

After ITensor builds, we will:

1. compile a tiny verification program;
2. implement the built-in Ising example;
3. reproduce the same MPO by hand.
