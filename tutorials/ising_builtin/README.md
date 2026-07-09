# Ising built-in track

This track follows the standard ITensor workflow:

1. define a `SpinHalf` site set;
2. build the transverse-field Ising Hamiltonian with `AutoMPO`;
3. convert it with `toMPO`;
4. run DMRG and inspect local observables.

The local driver is [main.cc](/Users/juan/Desktop/Git/FFD_GGE/tutorials/ising_builtin/src/main.cc).

## What to study in the code

- why the Hamiltonian uses `-4*J*Sz_j*Sz_{j+1}` and `-2*h*Sx_j`
- how `AutoMPO` translates operator strings into an MPO
- how `psi.position(j)` simplifies local expectation values

## Suggested experiments

1. Change `h` and observe the crossover from ordered to field-polarized behavior.
2. Change `N` and watch how `m_z(center)` and `m_x(center)` vary.
3. Compare `energy` against `inner(psi,H,psi)` and inspect `maxLinkDim(psi)`.
