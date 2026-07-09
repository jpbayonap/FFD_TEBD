# ITensor C++ setup on this machine

This project assumes:

- macOS
- Apple Clang
- the official ITensor C++ repository cloned into `vendor/ITensor`

## 1. Clone the library

```bash
cd /Users/juan/Desktop/Git/FFD_GGE
mkdir -p vendor
git clone --branch v3.2.0 --depth 1 https://github.com/ITensor/ITensor.git vendor/ITensor
```

## 2. Create the user options file

```bash
cd /Users/juan/Desktop/Git/FFD_GGE/vendor/ITensor
cp options.mk.sample options.mk
```

## 3. Edit `options.mk`

For this macOS setup, use the Clang configuration and the Accelerate framework:

- disable the default `g++` `CCCOM=...` line
- enable the `clang++` `CCCOM=...` line
- keep `PLATFORM=macos`
- keep `BLAS_LAPACK_LIBFLAGS=-framework Accelerate`

The official `options.mk.sample` in the ITensor `v3` branch uses:

- default GCC line: `CCCOM=g++ -m64 -std=c++17 -fconcepts -fPIC`
- macOS Clang alternative: `#CCCOM=clang++ -std=c++17 -fPIC -Wno-gcc-compat`
- macOS BLAS/LAPACK: `PLATFORM=macos`
- Accelerate link flags: `BLAS_LAPACK_LIBFLAGS=-framework Accelerate`

## 4. Build ITensor

```bash
cd /Users/juan/Desktop/Git/FFD_GGE/vendor/ITensor
make -j1
```

## 5. Return to the project root

```bash
cd /Users/juan/Desktop/Git/FFD_GGE
make check-itensor
```

If step 4 succeeds, the next command to try is:

```bash
make ising-builtin
```
## 6. Run a executable

./tutorials/ising_builtin/build/bin/ising_builtin
