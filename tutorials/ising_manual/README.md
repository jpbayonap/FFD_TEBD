# Ising manual track

Purpose: reproduce the same Ising-chain Hamiltonian without relying on the high-level MPO constructor.

Planned work:

1. define the local physical basis and operator tensors explicitly;
2. build site tensors and bond/link structure by hand;
3. assemble the MPO;
4. compare observables and matrix elements against the built-in version.
