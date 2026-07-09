#include "itensor/all.h"
#include "itensor/util/print_macro.h"

// Run example : make ising-builtin
// then ./tutorials/ising_builtin/build/bin/ising_builtin

using namespace itensor;

struct IsingParams
{
    int N = 24;
    Real J = 1.0;
    Real h = 2.0;
    int nsweeps = 6;
};

MPO makeTransverseFieldIsing(SiteSet const &sites, Real J, Real h)
{
    auto ampo = AutoMPO(sites);
    for (int j = 1; j < length(sites); ++j)
    {
        // SpinHalf uses S operators, not Pauli matrices:
        // sigma^z sigma^z = 4 Sz Sz and sigma^x = 2 Sx.
        ampo += -4.0 * J, "Sz", j, "Sz", j + 1;
    }
    for (int j = 1; j <= length(sites); ++j)
    {
        ampo += -2.0 * h, "Sx", j;
    }
    return toMPO(ampo);
}

MPS makeInitialState(SiteSet const &sites)
{
    auto state = InitState(sites);
    for (int j = 1; j <= length(sites); ++j)
    {
        state.set(j, (j % 2 == 1) ? "Up" : "Dn");
    }
    return randomMPS(state);
    //  strict Neel state with MPS(state);
}

Real measurePauliComponent(MPS &psi, SiteSet const &sites, std::string const &opname, int site)
{
    psi.position(site);
    auto ket = psi(site);
    auto bra = dag(prime(ket, "Site"));
    return 2.0 * elt(bra * op(sites, opname, site) * ket);
}

Real meanPauliComponent(MPS &psi, SiteSet const &sites, std::string const &opname)
{
    Real total = 0.0;
    for (int j = 1; j <= length(sites); ++j)
    {
        total += measurePauliComponent(psi, sites, opname, j);
    }
    return total / length(sites);
}

int main()
{
    auto p = IsingParams{};

    auto sites = SpinHalf(p.N, {"ConserveQNs=", false});
    auto H = makeTransverseFieldIsing(sites, p.J, p.h);
    auto psi0 = makeInitialState(sites);

    auto sweeps = Sweeps(p.nsweeps);
    sweeps.maxdim() = 10, 20, 40, 80, 120, 160;
    sweeps.cutoff() = 1E-10;
    sweeps.noise() = 1E-6, 1E-7, 1E-8, 0.0;

    auto [energy, psi] = dmrg(H, psi0, sweeps, {"Quiet=", true});

    auto center = p.N / 2;
    auto mz_center = measurePauliComponent(psi, sites, "Sz", center);
    auto mx_center = measurePauliComponent(psi, sites, "Sx", center);
    auto mz_mean = meanPauliComponent(psi, sites, "Sz");
    auto mx_mean = meanPauliComponent(psi, sites, "Sx");

    printfln("Transverse-field Ising chain");
    printfln("N = %d, J = %.6f, h = %.6f", p.N, p.J, p.h);
    printfln("Ground-state energy = %.12f", energy);
    printfln("Energy check         = %.12f", inner(psi, H, psi));
    printfln("maxLinkDim(psi)      = %d", maxLinkDim(psi));
    printfln("m_z(center)          = %.12f", mz_center);
    printfln("m_x(center)          = %.12f", mx_center);
    printfln("mean m_z             = %.12f", mz_mean);
    printfln("mean m_x             = %.12f", mx_mean);

    return 0;
}
