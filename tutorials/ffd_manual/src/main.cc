#include <iostream>
#include <fstream>
#include <iomanip>
#ifndef ACCELERATE_NEW_LAPACK
#define ACCELERATE_NEW_LAPACK
#endif
#include <Accelerate/Accelerate.h>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

using Complex = std::complex<double>;
// --------COMPLEX LINEAR ALGEBRA struct/functions---------

struct ComplexMatrix_col
{
    int rows = 0;
    int cols = 0;
    std::vector<Complex> data;

    ComplexMatrix_col() = default;

    ComplexMatrix_col(int r, int c)
        : rows(r), cols(c), data(r * c, Complex{0.0, 0.0})
    {
    }

    Complex &operator()(int i, int j)
    {
        if (i < 0 || i >= rows || j < 0 || j >= cols)
        {
            throw std::out_of_range("ComplexMatrix index out of range");
        }
        return data[j * rows + i];
    }

    Complex operator()(int i, int j) const
    {
        if (i < 0 || i >= rows || j < 0 || j >= cols)
        {
            throw std::out_of_range("ComplexMatrix index out of range");
        }
        return data[j * rows + i];
    }
};

struct ComplexSiteTensor
{
    int Dl = 1;
    int d = 2;
    int Dr = 1;
    std::vector<Complex> data;

    ComplexSiteTensor() = default;

    ComplexSiteTensor(int Dl_, int d_, int Dr_)
        : Dl(Dl_), d(d_), Dr(Dr_), data(Dl_ * d_ * Dr_, 0.0)
    {
    }
    // Dim= Dr Dl d = (Dr-1)dD_l +(d-1)D_L +(D_l-1)
    Complex &operator()(int l, int s, int r)
    {
        if (l < 0 || l >= Dl || s < 0 || s >= d || r < 0 || r >= Dr)
        {
            throw std::out_of_range("SiteTensor index out of range");
        }
        return data[(r * d + s) * Dl + l];
    }

    Complex operator()(int l, int s, int r) const
    {
        if (l < 0 || l >= Dl || s < 0 || s >= d || r < 0 || r >= Dr)
        {
            throw std::out_of_range("SiteTensor index out of range");
        }
        return data[(r * d + s) * Dl + l];
    }
};

struct CMPS
{
    int N = 0;
    // Spin 1/2 states
    int d = 2;
    // Ortogonality center
    int center = 0;
    std::vector<ComplexSiteTensor> A;
};

CMPS makeComplexProductState(const std::vector<int> &local_states, int d = 2)
{
    if (local_states.empty())
    {
        throw std::runtime_error("makeComplexProductState expects at least one site state");
    }

    CMPS psi;
    psi.N = static_cast<int>(local_states.size());
    psi.d = d;
    psi.center = 0;

    for (int s : local_states)
    {
        if (s < 0 || s >= d)
        {
            throw std::runtime_error("makeComplexProductState found invalid local state");
        }

        ComplexSiteTensor T(1, d, 1);
        T(0, s, 0) = Complex{1.0, 0.0};
        psi.A.push_back(T);
    }

    return psi;
}

// bond dimensions keeper
std::vector<int> ComplexBondDimensions(const CMPS &psi)
{
    if (static_cast<int>(psi.A.size()) != psi.N)
    {
        throw std::runtime_error("ComplexBondDimensions: inconsistent CMPS size");
    }

    std::vector<int> dims;
    for (int j = 0; j + 1 < psi.N; ++j)
    {
        dims.push_back(psi.A[j].Dr);
    }
    return dims;
}

CMPS makeTiltedProductState(int N, double theta)
{
    if (N <= 0)
    {
        throw std::runtime_error("makeTiltedProducState: uncompatible system size");
    }

    CMPS psi;
    psi.N = N;
    psi.d = 2;
    psi.center = 0;

    double amp0 = std::sin(theta); // |0 > coefficient
    double amp1 = std::cos(theta); // |1 > coefficient

    for (int j = 0; j < N; ++j)
    {
        ComplexSiteTensor T(1, 2, 1);
        T(0, 0, 0) = Complex{amp0, 0.0};
        T(0, 1, 0) = Complex{amp1, 0.0};
        // Update MPS
        psi.A.push_back(T);
    }

    return psi;
}

std::vector<double> makeStaggeredSqrtCouplings(
    int N,
    double alpha,
    double beta,
    double gamma)
{
    if (N <= 0)
    {
        throw std::runtime_error("makeStaggeredSqrtCouplings: invalid N");
    }

    std::vector<double> b(N);

    for (int m = 0; m < N; ++m)
    {
        int r = m % 3;
        b[m] = (r == 0)   ? std::sqrt(alpha)
               : (r == 1) ? std::sqrt(beta)
                          : std::sqrt(gamma);
    }
    return b;
}

// print helper
void printIntVector(const std::vector<int> &v, const std::string &name)
{
    std::cout << name << ":";
    for (int x : v)
    {
        std::cout << " " << x;
    }
    std::cout << "\n";
}

void printDoubleVector(const std::vector<double> &v, const std::string &name)
{
    std::cout << name << ":";
    for (double x : v)
    {
        std::cout << " " << x;
    }
    std::cout << "\n";
}

ComplexMatrix_col cmatmul(const ComplexMatrix_col &A, const ComplexMatrix_col &B)
{
    if (A.cols != B.rows)
    {
        throw std::runtime_error("cmatmul dimension mismatch");
    }

    ComplexMatrix_col C(A.rows, B.cols);
    Complex alpha{1.0, 0.0};
    Complex beta{0.0, 0.0};

    cblas_zgemm(
        CblasColMajor,
        CblasNoTrans,
        CblasNoTrans,
        A.rows, // m
        B.cols, // n
        A.cols, // k
        &alpha,
        A.data.data(), // data pointer A
        A.rows,        // dla
        B.data.data(), // data pointer B
        B.rows,        // ldb
        &beta,
        C.data.data(), // data pointer C
        C.rows         // ldf
    );

    // O(mnk) reference implementation
    //  for (int j = 0; j < B.cols; ++j)
    //  {
    //    for (int k = 0; k < A.cols; ++k)
    //    {
    //      Complex Bkj = B(k, j);
    //      for (int i = 0; i < A.rows; ++i)
    //      {
    //        C(i, j) += A(i, k) * Bkj;
    //      }
    //    }
    //  }

    return C;
}

ComplexMatrix_col complexKron(const ComplexMatrix_col &A, const ComplexMatrix_col &B)
{
    int r = A.rows * B.rows;
    int c = A.cols * B.cols;
    ComplexMatrix_col ComplexKron_AB(r, c);
    for (int iA = 0; iA < A.rows; ++iA)
    {
        for (int iB = 0; iB < B.rows; ++iB)
        {
            for (int jA = 0; jA < A.cols; ++jA)
            {
                for (int jB = 0; jB < B.cols; ++jB)
                {
                    int row = iA * B.rows + iB;
                    int col = jA * B.cols + jB;
                    ComplexKron_AB(row, col) = A(iA, jA) * B(iB, jB);
                }
            }
        }
    }
    return ComplexKron_AB;
}

ComplexMatrix_col add(const ComplexMatrix_col &A, const ComplexMatrix_col &B)
{

    if (A.rows != B.rows || A.cols != B.cols)
    {
        throw std::runtime_error(" dimension mismatch");
    }

    ComplexMatrix_col C(A.rows, A.cols);
    for (int i = 0; i < A.rows; ++i)
    {
        for (int j = 0; j < A.cols; ++j)
        {
            C(i, j) = A(i, j) + B(i, j);
        }
    }
    return C;
}

ComplexMatrix_col scale(double alpha, const ComplexMatrix_col &A)
{
    ComplexMatrix_col B(A.rows, A.cols);
    for (int i = 0; i < A.rows; ++i)
    {
        for (int j = 0; j < A.cols; ++j)
        {
            B(i, j) = alpha * A(i, j);
        }
    }
    return B;
}
// ------- END COMPLEX LINEAR ALGEBRA struct/functions ---------

// ----SIMPLE HELPERS----
void printComplexMatrix_col(const ComplexMatrix_col &A, const std::string &name)
{
    std::cout << name << "(" << A.rows << "x" << A.cols << ")\n";
    for (int i = 0; i < A.rows; ++i)
    {
        for (int j = 0; j < A.cols; ++j)
        {
            Complex z = A(i, j);
            std::cout << std::setw(10) << z.real()
                      << (z.imag() < 0.0 ? " - " : " + ")
                      << std::setw(10) << std::abs(z.imag()) << "i ";
        }
        std::cout << "\n";
    }
}

ComplexMatrix_col makeComplexIdentity(int n)
{
    if (n < 0)
    {
        throw std::runtime_error("makeComplexIdentity expects n >= 0");
    }
    ComplexMatrix_col I(n, n);
    for (int i = 0; i < n; ++i)
    {
        I(i, i) = Complex{1.0, 0.0};
    }
    return I;
}

ComplexMatrix_col conjTranspose(const ComplexMatrix_col &A)
{
    ComplexMatrix_col AH(A.cols, A.rows);
    for (int i = 0; i < A.rows; ++i)
    {
        for (int j = 0; j < A.cols; ++j)
        {
            AH(j, i) = std::conj(A(i, j));
        }
    }
    return AH;
}

double maxAbsDiffComplex(const ComplexMatrix_col &A, const ComplexMatrix_col &B)
{
    if (A.rows != B.rows || A.cols != B.cols)
    {
        throw std::runtime_error("maxAbsDiffComplex dimension mismatch");
    }

    double ans = 0.0;
    for (int i = 0; i < A.rows; ++i)
    {
        for (int j = 0; j < A.cols; ++j)
        {
            ans = std::max(ans, std::abs(A(i, j) - B(i, j)));
        }
    }
    return ans;
}

struct ComplexSVDResult
{
    ComplexMatrix_col U;
    std::vector<double> S;
    ComplexMatrix_col VH;
};

struct ComplexSplitStateResult
{
    ComplexMatrix_col left;
    ComplexMatrix_col right;
    std::vector<double> schmidt;
    double discarded_weight = 0.0;
    double discarded_fraction = 0.0;
};

struct ComplexThreeSiteSplitResult
{
    ComplexSiteTensor A0;
    ComplexSiteTensor A1;
    ComplexSiteTensor A2;

    std::vector<double> schmidt01;
    std::vector<double> schmidt12;

    double discarded_weight01 = 0.0;
    double discarded_weight12 = 0.0;

    double discarded_fraction01 = 0.0;
    double discarded_fraction12 = 0.0;
};

void printComplexSplitReport(
    const ComplexSplitStateResult &split,
    const std::string &name)
{
    std::cout << name
              << " chi =" << split.schmidt.size()
              << " discarded_weight=" << split.discarded_weight
              << " discarded_fraction=" << split.discarded_fraction
              << "\n";

    printDoubleVector(split.schmidt, name + " singular values");
}

ComplexMatrix_col
makeComplexDiagonalFromReal(const std::vector<double> &s)
{
    int k = static_cast<int>(s.size());
    ComplexMatrix_col D(k, k);
    for (int i = 0; i < k; ++i)
    {
        D(i, i) = Complex{s[i], 0.0};
    }
    return D;
}

ComplexSVDResult complexSvd(const ComplexMatrix_col &A)
{
    using LapackInt = int;

    LapackInt m = A.rows;
    LapackInt n = A.cols;
    LapackInt k = std::min(m, n);

    ComplexMatrix_col Acopy = A;
    ComplexMatrix_col U(m, k);
    ComplexMatrix_col VH(k, n);
    std::vector<double> S(k);

    char jobu = 'S';
    char jobvt = 'S';

    LapackInt lda = m;
    LapackInt ldu = m;
    LapackInt ldvh = k;
    LapackInt info = 0;

    Complex work_query{0.0, 0.0};
    LapackInt lwork = -1;
    std::vector<double> rwork(std::max(1, 5 * k));

    zgesvd_(
        &jobu,
        &jobvt,
        &m,
        &n,
        Acopy.data.data(),
        &lda,
        S.data(),
        U.data.data(),
        &ldu,
        VH.data.data(),
        &ldvh,
        &work_query,
        &lwork,
        rwork.data(),
        &info);

    if (info != 0)
    {
        throw std::runtime_error("complexSvd workspace query failed");
    }

    lwork = static_cast<LapackInt>(work_query.real());
    if (lwork < 1)
    {
        lwork = 1;
    }

    std::vector<Complex> work(lwork);

    zgesvd_(
        &jobu,
        &jobvt,
        &m,
        &n,
        Acopy.data.data(),
        &lda,
        S.data(),
        U.data.data(),
        &ldu,
        VH.data.data(),
        &ldvh,
        work.data(),
        &lwork,
        rwork.data(),
        &info);

    if (info != 0)
    {
        throw std::runtime_error("complexSvd failed");
    }

    return {U, S, VH};
}

double schmidtNorm2(const std::vector<double> &s)
{
    double norm2 = 0.0;
    for (double lambda : s)
    {
        norm2 += lambda * lambda;
    }
    return norm2;
}

// Keep by Cutoff bond dimension

int chooseKeepByDiscardedFraction(
    const std::vector<double> &s,
    int chi_max,
    double cutoff)
{
    if (chi_max <= 0)
    {
        throw ::std::runtime_error("chooseKeepByDiscardedFraction: chi_max >0 is expected");
    }

    int kfull = static_cast<int>(s.size());
    int hard_cap = std::min(chi_max, kfull);

    double total_weight = schmidtNorm2(s);
    if (total_weight <= 0.0)
    {
        return 1;
    }

    for (int keep = 1; keep <= hard_cap; ++keep)
    {
        double discarded = 0.0;
        for (int a = keep; a < kfull; ++a)
        {
            discarded += s[a] * s[a];
        }

        double fraction = discarded / total_weight;
        if (fraction <= cutoff)
        {
            return keep;
        }
    }

    return hard_cap;
}

ComplexSplitStateResult ComplexSplitState(
    const ComplexMatrix_col &Theta,
    int chi,
    bool left_to_right,
    double cutoff = -1.0)
{
    if (chi <= 0)
    {
        throw std::runtime_error("ComplexSplitState expects chi > 0");
    }

    auto full = complexSvd(Theta);
    int kfull = static_cast<int>(full.S.size());
    //     cutoff < 0       old fixed behavior
    // cutoff >= 0      dynamic bond dimension behavior
    int keep = (cutoff < 0.0) ? std::min(chi, kfull) : chooseKeepByDiscardedFraction(full.S, chi, cutoff);

    ComplexMatrix_col U_trunc(full.U.rows, keep);
    ComplexMatrix_col VH_trunc(keep, full.VH.cols);
    std::vector<double> S_trunc(keep);

    for (int a = 0; a < keep; ++a)
    {
        S_trunc[a] = full.S[a];
    }
    for (int i = 0; i < full.U.rows; ++i)
    {
        for (int a = 0; a < keep; ++a)
        {
            U_trunc(i, a) = full.U(i, a);
        }
    }
    for (int a = 0; a < keep; ++a)
    {
        for (int j = 0; j < full.VH.cols; ++j)
        {
            VH_trunc(a, j) = full.VH(a, j);
        }
    }

    double discarded_weight = 0.0;
    for (int a = keep; a < kfull; ++a)
    {
        discarded_weight += full.S[a] * full.S[a];
    }

    double total_weight = schmidtNorm2(full.S);
    double discarded_fraction = (total_weight > 0.0) ? discarded_weight / total_weight : 0.0;

    ComplexMatrix_col S = makeComplexDiagonalFromReal(S_trunc);
    ComplexMatrix_col left_trunc;
    ComplexMatrix_col right_trunc;

    if (left_to_right)
    {
        left_trunc = U_trunc;
        right_trunc = cmatmul(S, VH_trunc);
    }
    else
    {
        left_trunc = cmatmul(U_trunc, S);
        right_trunc = VH_trunc;
    }

    return {left_trunc, right_trunc, S_trunc, discarded_weight, discarded_fraction};
}

ComplexMatrix_col makePauliX()
{
    ComplexMatrix_col X(2, 2);
    X(0, 1) = 1.0;
    X(1, 0) = 1.0;
    return X;
}

ComplexMatrix_col makePauliY()
{
    ComplexMatrix_col Y(2, 2);
    Y(0, 1) = Complex{0.0, -1.0};
    Y(1, 0) = Complex{0.0, 1.0};
    return Y;
}

ComplexMatrix_col makePauliZ()
{
    ComplexMatrix_col Z(2, 2);
    Z(0, 0) = 1.0;
    Z(1, 1) = -1.0;
    return Z;
}

ComplexMatrix_col makeFFD3SiteHamiltonian(double c_zzx, double c_xzz)
{

    auto X = makePauliX();
    auto Z = makePauliZ();

    ComplexMatrix_col ZZX = complexKron(complexKron(Z, Z), X);
    ComplexMatrix_col XZZ = complexKron(complexKron(X, Z), Z);

    return add(scale(c_zzx, ZZX), scale(c_xzz, XZZ));
}

ComplexMatrix_col makeRealTimePauliStringGate(
    const ComplexMatrix_col &P,
    double coeff,
    double dt)
{
    if (P.rows != P.cols)
    {
        throw std::runtime_error("makeRealTimePauliStringGate expects square matrices");
    }

    auto I = makeComplexIdentity(P.rows);
    ComplexMatrix_col G(P.rows, P.cols);

    Complex minus_i{0.0, -1.0};
    double theta = coeff * dt;
    for (int row = 0; row < P.rows; ++row)
    {
        for (int col = 0; col < P.cols; ++col)
        {
            G(row, col) = std::cos(theta) * I(row, col) + minus_i * std::sin(theta) * P(row, col);
        }
    }
    return G;
}

// Fendley Hamiltonian convention
// m = 0: h_0 = X_0
// m = 1: h_1 = Z_0 X_1
// m >= 2: h_m = Z_{m-2} Z_{m-1} X_m

ComplexMatrix_col makeFFDOneSiteTermOperator()
{
    return makePauliX();
}

ComplexMatrix_col makeFFDOneSiteTermGate(double b, double dt)
{
    auto X = makePauliX();
    return makeRealTimePauliStringGate(X, b, dt);
}

ComplexMatrix_col makeFFDTwoSiteTermOperator()
{
    auto Z = makePauliZ();
    auto X = makePauliX();
    auto ZX = complexKron(Z, X);
    return ZX;
}

ComplexMatrix_col makeFFDTwoSiteTermGate(double b, double dt)
{
    auto Z = makePauliZ();
    auto X = makePauliX();
    auto ZX = complexKron(Z, X);
    return makeRealTimePauliStringGate(ZX, b, dt);
}

ComplexMatrix_col makeFFDThreeSiteTermOperator()
{
    auto Z = makePauliZ();
    auto X = makePauliX();
    auto ZZX = complexKron(Z, complexKron(Z, X));
    return ZZX;
}

ComplexMatrix_col makeFFDThreeSiteTermGate(double b, double dt)
{
    auto Z = makePauliZ();
    auto X = makePauliX();
    auto ZZX = complexKron(Z, complexKron(Z, X));
    return makeRealTimePauliStringGate(ZZX, b, dt);
}

ComplexMatrix_col makeFFD3SiteGate(double c_zzx, double c_xzz, double dt)
{
    auto X = makePauliX();
    auto Z = makePauliZ();

    auto ZZX = complexKron(complexKron(Z, Z), X);
    auto XZZ = complexKron(complexKron(X, Z), Z);

    auto G_zzx = makeRealTimePauliStringGate(ZZX, c_zzx, dt);
    auto G_xzz = makeRealTimePauliStringGate(XZZ, c_xzz, dt);

    return cmatmul(G_zzx, G_xzz);
}

// ----------------- END SIMPLE HELPERS ---------------

// ----------------- TEBD HELPERS -----------------------
ComplexMatrix_col ComplexBuildThreeSiteTheta(
    const ComplexSiteTensor &A0,
    const ComplexSiteTensor &A1,
    const ComplexSiteTensor &A2)
{
    if (A0.Dr != A1.Dl || A1.Dr != A2.Dl)
    {
        throw std::runtime_error("ComplexBuildThreeSiteTheta bond dimension mismatch");
    }

    if (A0.d != A1.d || A1.d != A2.d)
    {
        throw std::runtime_error("ComplexBuildThreeSiteTheta physical dimension mismatch");
    }

    int Dl = A0.Dl;
    int D1 = A0.Dr;
    int D2 = A1.Dr;
    int Dr = A2.Dr;
    int d = A0.d;
    // Group in (Dls1, s2s3Dr) matrix
    ComplexMatrix_col Theta(Dl * d, d * d * Dr);

    for (int l = 0; l < Dl; ++l)
    {
        for (int r = 0; r < Dr; ++r)
        {
            for (int s0 = 0; s0 < d; ++s0)
            {
                for (int s1 = 0; s1 < d; ++s1)
                {
                    for (int s2 = 0; s2 < d; ++s2)
                    {
                        Complex val = 0.0;
                        // double sum over m1/2 in [0,D1/2)
                        for (int m1 = 0; m1 < D1; ++m1)
                        {
                            for (int m2 = 0; m2 < D2; ++m2)
                            {
                                val += A0(l, s0, m1) * A1(m1, s1, m2) * A2(m2, s2, r);
                            }
                        }
                        int row = l * d + s0;
                        // Same indexing as Complex site tensor
                        int col = (s1 * d + s2) * Dr + r;
                        Theta(row, col) = val;
                    }
                }
            }
        }
    }

    return Theta;
}

void ComplexWriteBondFromSplit(CMPS &psi, int j, const ComplexSplitStateResult &split)
{
    int Dl = psi.A[j].Dl;
    int d = psi.A[j].d;
    int Dr = psi.A[j + 1].Dr;
    int chi = split.schmidt.size(); // BOND DIMENSION!

    // updated tensors
    ComplexSiteTensor Aj_new(Dl, d, chi);
    ComplexSiteTensor Aj1_new(chi, d, Dr);

    // Shape checks
    if (split.left.rows != Dl * d)
    {
        throw std::runtime_error("Split mismatch for the left tensor row entries");
    }
    if (split.left.cols != chi)
    {
        throw std::runtime_error("Split mismatch for the left tensor col entries ");
    }
    if (split.right.rows != chi)
    {
        throw std::runtime_error("Split mismatch for the right tensor row entries");
    }
    if (split.right.cols != d * Dr)
    {
        throw std::runtime_error("Split mismatch for the right tensor col entries ");
    }

    // left tensor update Aj
    for (int l = 0; l < Dl; ++l)
    {
        for (int sj = 0; sj < d; ++sj)
        {
            for (int a = 0; a < chi; ++a)
            {
                int row = l * d + sj;
                Aj_new(l, sj, a) = split.left(row, a);
            }
        }
    }
    // right tensor update Aj1
    for (int a = 0; a < chi; ++a)
    {
        for (int sj1 = 0; sj1 < d; ++sj1)
        {
            for (int r = 0; r < Dr; ++r)
            {
                int col = sj1 * Dr + r;
                Aj1_new(a, sj1, r) = split.right(a, col);
            }
        }
    }

    // Update state psi
    psi.A[j] = Aj_new;
    psi.A[j + 1] = Aj1_new;
}

ComplexMatrix_col ComplexBuildTwoSiteTheta(
    const ComplexSiteTensor &A0,
    const ComplexSiteTensor &A1)
{
    // left side (l,s1)
    // right side (s2,r)

    ComplexMatrix_col Theta(A0.Dl * A0.d, A1.d * A1.Dr);

    // Check
    if (A0.Dr != A1.Dl)
    {
        throw std::runtime_error("buildBondTheta bond mismatch");
    }

    int S1 = A0.d;
    int S2 = A1.d;
    // sum index for Theta(l,s1,s2,r) = sum_m A0(l,s1,m) * A1(m,s2,r)
    int Dm = A0.Dr;

    for (int l = 0; l < A0.Dl; ++l)
    {
        for (int r = 0; r < A1.Dr; ++r)
        {
            for (int s1 = 0; s1 < S1; ++s1)
            {
                for (int s2 = 0; s2 < S2; ++s2)
                {
                    Complex val = 0.0;
                    for (int m = 0; m < Dm; ++m)
                    {
                        val += A0(l, s1, m) * A1(m, s2, r);
                    }
                    int row = l * S1 + s1;
                    int col = s2 * A1.Dr + r;
                    Theta(row, col) = val;
                }
            }
        }
    }
    return Theta;
}

ComplexMatrix_col ComplexApplyTwoSiteGateToTheta(int Dl, int Dr, int d, const ComplexMatrix_col &Theta, const ComplexMatrix_col &G2)
{
    if (Theta.rows != Dl * d || Theta.cols != d * Dr)
    {
        throw std::runtime_error("applyGateToTheta shape mismatch");
    }
    if (G2.rows != d * d || G2.cols != d * d)
    {
        throw std::runtime_error("applyGateToTheta expects a (d*d)x(d*d) gate");
    }

    ComplexMatrix_col Theta_new(Dl * d, d * Dr);

    for (int l = 0; l < Dl; ++l)
    {
        for (int r = 0; r < Dr; ++r)
        {
            //  vector to apply G2
            ComplexMatrix_col theta_lr(d * d, 1);

            for (int s1 = 0; s1 < d; ++s1)
            {
                for (int s2 = 0; s2 < d; ++s2)
                {
                    int row = l * d + s1;
                    int col = s2 * Dr + r;
                    int flat = d * s1 + s2;
                    theta_lr(flat, 0) = Theta(row, col);
                }
            }
            // Apply gate G2 |theta_lr >

            ComplexMatrix_col theta_lr_new = cmatmul(G2, theta_lr);

            for (int s1 = 0; s1 < d; ++s1)
            {
                for (int s2 = 0; s2 < d; ++s2)
                {
                    int flat = d * s1 + s2;
                    int row = l * d + s1;
                    int col = s2 * Dr + r;
                    Theta_new(row, col) = theta_lr_new(flat, 0);
                }
            }
        }
    }

    return Theta_new;
}

// Forward declaration: implemented below after the canonical site helpers.
void prepareCanonicalBlock(CMPS &psi, int first_site, int block_size);
void prepareCanonicalBlockFast(CMPS &psi, int first_site, int block_size);

ComplexSplitStateResult ComplexUpdateTwoSite(
    CMPS &psi,
    int j,
    const ComplexMatrix_col &G2,
    int chi,
    double cutoff = -1.0)
{
    if (j < 0 || j + 1 >= psi.N)
    {
        throw std::runtime_error("ComplexUpdateTwoSite bond index out of range. ");
    }
    // Debugg: Slow
    // prepareCanonicalBlock(psi, j, 2);

    // Fast implementation
    prepareCanonicalBlockFast(psi, j, 2);

    auto Theta = ComplexBuildTwoSiteTheta(psi.A[j], psi.A[j + 1]);

    int Dl = psi.A[j].Dl;
    int d = psi.A[j].d;
    int Dr = psi.A[j + 1].Dr;

    // SVD with chi bond dimension
    auto Theta_gate = ComplexApplyTwoSiteGateToTheta(Dl, Dr, d, Theta, G2);
    auto split = ComplexSplitState(Theta_gate, chi, true, cutoff);

    ComplexWriteBondFromSplit(psi, j, split);
    // Update orthogonality center
    psi.center = j + 1;

    return split;
}

ComplexMatrix_col ComplexApplyThreeSiteGateToTheta(
    int Dl,
    int Dr,
    int d,
    const ComplexMatrix_col &Theta,
    const ComplexMatrix_col &G3)
{
    if (Theta.rows != Dl * d || Theta.cols != d * d * Dr)
    {
        throw std::runtime_error("ComplexApplyThreeSiteGateToTheta shape mismatch");
    }
    if (G3.rows != d * d * d || G3.cols != d * d * d)
    {
        throw std::runtime_error("ComplexApplyThreeSiteGateToTheta expects 8x8 Gate");
    }

    ComplexMatrix_col Theta_new(Dl * d, d * d * Dr);

    for (int l = 0; l < Dl; ++l)
    {
        for (int r = 0; r < Dr; ++r)
        {
            // Vector used to apply gate
            ComplexMatrix_col theta_lr(d * d * d, 1);

            for (int s0 = 0; s0 < d; ++s0)
            {
                for (int s1 = 0; s1 < d; ++s1)
                {
                    for (int s2 = 0; s2 < d; ++s2)
                    {
                        int flat = (s0 * d + s1) * d + s2;
                        int row = l * d + s0;
                        int col = (s1 * d + s2) * Dr + r;
                        // Update vector
                        theta_lr(flat, 0) = Theta(row, col);
                    }
                }
            }

            // Apply Gate
            auto theta_lr_new = cmatmul(G3, theta_lr);
            // update Theta
            for (int s0 = 0; s0 < d; ++s0)
            {
                for (int s1 = 0; s1 < d; ++s1)
                {
                    for (int s2 = 0; s2 < d; ++s2)
                    {
                        int flat = (s0 * d + s1) * d + s2;
                        int row = l * d + s0;
                        int col = (s1 * d + s2) * Dr + r;
                        Theta_new(row, col) = theta_lr_new(flat, 0);
                    }
                }
            }
        }
    }
    return Theta_new;
}

void ComplexApplyOneSiteGate(
    CMPS &psi,
    int j,
    const ComplexMatrix_col &G1)
{
    if (j < 0 || j >= psi.N)
    {
        throw std::runtime_error("ComplexApplyOneSiteGate: invalid site index");
    }

    auto &A = psi.A[j];
    int Dl = A.Dl;
    int d = A.d;
    int Dr = A.Dr;

    if (G1.rows != d || G1.cols != d)
    {
        throw std::runtime_error("ComplexApplyOneSiteGate: gate shape mismatch");
    }

    ComplexSiteTensor A_new(Dl, d, Dr);

    for (int l = 0; l < Dl; ++l)
    {
        for (int r = 0; r < Dr; ++r)
        {
            for (int s_new = 0; s_new < d; ++s_new)
            {
                Complex val = 0.0;
                for (int s_old = 0; s_old < d; ++s_old)
                {
                    val += G1(s_new, s_old) * A(l, s_old, r);
                }
                // update
                A_new(l, s_new, r) = val;
            }
        }
    }
    psi.A[j] = A_new;
}

// ------------- CANONICALIZATION HELPERS ------------

void leftCanonicalizeSite(CMPS &psi, int j)
{
    if (j < 0 || j + 1 >= psi.N)
    {
        throw std::runtime_error("leftCanonicalizeSite: index out of range");
    }
    // Left canonical if : A_j^\dagger A_j = I
    // M(row, r) = A_j(l, s, r)
    // row = l * d + s

    auto A = psi.A[j];
    auto Ap = psi.A[j + 1];
    if (A.Dr != Ap.Dl)
    {
        throw std::runtime_error("leftCanonicalizeSite: bond mismatch");
    }

    int d = A.d;
    int Dl = A.Dl;
    int Dr = A.Dr;

    int R = Dl * d;
    ComplexMatrix_col M(R, Dr);

    for (int l = 0; l < Dl; ++l)
    {
        for (int s = 0; s < d; ++s)
        {
            for (int r = 0; r < Dr; ++r)
            {
                int row = l * d + s;
                M(row, r) = A(l, s, r);
            }
        }
    }

    auto svdM = complexSvd(M);
    int Dnew = static_cast<int>(svdM.S.size());
    auto U = svdM.U;
    auto B = cmatmul(makeComplexDiagonalFromReal(svdM.S), svdM.VH);
    auto A_j_new = ComplexSiteTensor(Dl, d, Dnew);
    auto A_jp1_new = ComplexSiteTensor(Dnew, Ap.d, Ap.Dr);

    for (int l = 0; l < Dl; ++l)
    {
        for (int s = 0; s < d; ++s)
        {
            for (int a = 0; a < Dnew; ++a)
            {
                // Set the current site to U
                int row = l * d + s;
                A_j_new(l, s, a) = U(row, a);
            }
        }
    }
    for (int a = 0; a < Dnew; ++a)
    {
        for (int sp = 0; sp < d; ++sp)
        {
            for (int rp = 0; rp < Ap.Dr; ++rp)
            {
                Complex val = 0.0;
                //  absorb SVH into the next tensor
                for (int m = 0; m < Dr; ++m)
                {
                    val += B(a, m) * Ap(m, sp, rp);
                }
                A_jp1_new(a, sp, rp) += val;
            }
        }
    }
    // update
    psi.A[j] = A_j_new;
    psi.A[j + 1] = A_jp1_new;
}

void rightCanonicalizeSite(CMPS &psi, int j)
{
    if (j <= 0 || j >= psi.N)
    {
        throw std::runtime_error("rightCanonicalizeSite: index out of range");
    }
    // Left canonical if : A_j^\dagger A_j = I
    // M(l, col) = A_j(l, s, r)
    // col = s * Dr + r

    auto A = psi.A[j];
    auto Am = psi.A[j - 1];
    if (A.Dl != Am.Dr)
    {
        throw std::runtime_error("rightCanonicalizeSite: bond mismatch");
    }

    int d = A.d;
    int Dl = A.Dl;
    int Dr = A.Dr;

    int C = d * Dr;
    ComplexMatrix_col M(Dl, C);

    for (int l = 0; l < Dl; ++l)
    {

        for (int s = 0; s < d; ++s)
        {
            for (int r = 0; r < Dr; ++r)
            {
                int col = s * Dr + r;
                M(l, col) = A(l, s, r);
            }
        }
    }

    auto svdM = complexSvd(M);
    int Dnew = static_cast<int>(svdM.S.size());
    auto VH = svdM.VH;
    auto B = cmatmul(svdM.U, makeComplexDiagonalFromReal(svdM.S));
    auto A_j_new = ComplexSiteTensor(Dnew, d, Dr);
    auto A_jm1_new = ComplexSiteTensor(Am.Dl, Am.d, Dnew);

    for (int s = 0; s < d; ++s)
    {
        for (int r = 0; r < Dr; ++r)
        {
            for (int a = 0; a < Dnew; ++a)
            {
                // Set the current site to VH
                int col = s * Dr + r;
                A_j_new(a, s, r) = VH(a, col);
            }
        }
    }

    for (int a = 0; a < Dnew; ++a)
    {
        for (int sm = 0; sm < d; ++sm)
        {
            for (int lm = 0; lm < Am.Dl; ++lm)
            {
                Complex val = 0.0;
                //  absorb SVH into the next tensor
                for (int m = 0; m < Dl; ++m)
                {
                    val += Am(lm, sm, m) * B(m, a);
                }
                A_jm1_new(lm, sm, a) += val;
            }
        }
    }
    // update
    psi.A[j] = A_j_new;
    psi.A[j - 1] = A_jm1_new;
}

void prepareCanonicalBlock(CMPS &psi, int first_site, int block_size)
{
    if (block_size <= 0)
    {
        throw std::runtime_error("prepareCanonicalBlock: block_size must be positive");
    }

    int last_site = first_site + block_size - 1;

    if (first_site < 0 || last_site >= psi.N)
    {
        throw std::runtime_error("prepareCanonicalBlock: block out of range");
    }

    // Sites to the left of the block become left canonical
    for (int j = 0; j < first_site; ++j)
    {
        leftCanonicalizeSite(psi, j);
    }
    // Sites to the right of the block become right canonical
    for (int j = psi.N - 1; j > last_site; --j)
    {
        rightCanonicalizeSite(psi, j);
    }
}

void moveOrthogonalityCenter(CMPS &psi, int target)
{
    if (target < 0 || target >= psi.N)
    {
        throw std::runtime_error("moveOrthogonalityCenter: target out of range");
    }

    if (psi.center < 0 || psi.center >= psi.N)
    {
        throw std::runtime_error("moveOrthogonalityCenter: MPS center out of range");
    }

    while (psi.center < target)
    {
        leftCanonicalizeSite(psi, psi.center);
        psi.center += 1;
    }

    while (psi.center > target)
    {
        rightCanonicalizeSite(psi, psi.center);
        psi.center -= 1;
    }
}

void prepareCanonicalBlockFast(CMPS &psi, int first_site, int block_size)
{
    if (block_size <= 0)
    {
        throw std::runtime_error("prepareCanonicalBlockFast: block_size should be positive");
    }

    int last_site = first_site + block_size - 1;

    if (first_site < 0 || last_site >= psi.N)
    {
        throw std::runtime_error("prepareCanonicalBlockFast: block out of range");
    }

    moveOrthogonalityCenter(psi, first_site);
}

ComplexThreeSiteSplitResult ComplexSplitStateThreeSiteThetaLeftToRight(
    const ComplexMatrix_col &Theta,
    int Dl,
    int d,
    int Dr,
    int chi,
    double cutoff = -1.0)
{
    if (Theta.rows != Dl * d || Theta.cols != d * d * Dr)
    {
        throw std::runtime_error("ComplexApplyThreeSiteThetaLeftToRight shape mismatch");
    }
    // First split
    auto split01 = ComplexSplitState(Theta, chi, true, cutoff);
    int chi01 = static_cast<int>(split01.schmidt.size());
    // First tensor: split01.left has shape (Dl*d) x chi01.
    // (A0_new) (Dl, d, chi01)
    ComplexSiteTensor A0_new(Dl, d, chi01);

    for (int l = 0; l < Dl; ++l)
    {
        for (int s0 = 0; s0 < d; ++s0)
        {
            for (int a = 0; a < chi01; ++a)
            {
                int row = l * d + s0;
                A0_new(l, s0, a) = split01.left(row, a);
            }
        }
    }

    // split01.right shape = chi01 x (d * d * Dr)
    // reshape 01.right to  (a, s1) | (s2, r)
    ComplexMatrix_col Theta12(chi01 * d, d * Dr);

    for (int a = 0; a < chi01; ++a)
    {
        for (int s1 = 0; s1 < d; ++s1)
        {
            for (int s2 = 0; s2 < d; ++s2)
            {
                for (int r = 0; r < Dr; ++r)
                {
                    int old_col = (s1 * d + s2) * Dr + r;
                    int new_row = a * d + s1;
                    int new_col = s2 * Dr + r;
                    Theta12(new_row, new_col) = split01.right(a, old_col);
                }
            }
        }
    }

    auto split12 = ComplexSplitState(Theta12, chi, true, cutoff);
    int chi12 = static_cast<int>(split12.schmidt.size());

    // Second tensor: split12.left has shape (chi01 *d) x chi12.
    // (A1_new) (chi01, d, chi12)
    ComplexSiteTensor A1_new(chi01, d, chi12);

    for (int a = 0; a < chi01; ++a)
    {
        for (int s1 = 0; s1 < d; ++s1)
        {
            for (int b = 0; b < chi12; ++b)
            {
                int row = a * d + s1;
                A1_new(a, s1, b) = split12.left(row, b);
            }
        }
    }
    // Third tensor: split12.right has shape chi12 x (d*Dr).

    // (A2_new) (chi12, d, Dr)
    ComplexSiteTensor A2_new(chi12, d, Dr);

    for (int b = 0; b < chi12; ++b)
    {
        for (int s2 = 0; s2 < d; ++s2)
        {
            for (int r = 0; r < Dr; ++r)
            {
                int col = s2 * Dr + r;
                A2_new(b, s2, r) = split12.right(b, col);
            }
        }
    }

    return {
        A0_new,
        A1_new,
        A2_new,
        split01.schmidt,
        split12.schmidt,
        split01.discarded_weight,
        split12.discarded_weight,
        split01.discarded_fraction,
        split12.discarded_fraction};
}

ComplexThreeSiteSplitResult ComplexUpdateThreeSite(
    CMPS &psi,
    int j,
    const ComplexMatrix_col &G3,
    int chi,
    double cutoff = -1.0)
{
    // 1. check j is valid:

    if (j < 0 || j + 2 >= psi.N)
    {
        throw std::runtime_error("ComplexUpdateThreeSite: invalid three-site index");
    }
    if (static_cast<int>(psi.A.size()) != psi.N)
    {
        throw std::runtime_error("ComplexUpdateThreeSite: inconsistent CMPS size");
    }

    // 2. Build theta
    // Slow canonicalization
    // prepareCanonicalBlock(psi, j, 3);
    // Fast orthogonality center
    prepareCanonicalBlockFast(psi, j, 3);

    auto Theta = ComplexBuildThreeSiteTheta(psi.A[j], psi.A[j + 1], psi.A[j + 2]);

    // 3. Read dimensions Dl,  d, Dr
    int Dl = psi.A[j].Dl;
    int d = psi.A[j].d;
    int Dr = psi.A[j + 2].Dr;

    // 4. Apply G4

    auto Theta_gate = ComplexApplyThreeSiteGateToTheta(Dl, Dr, d, Theta, G3);

    // 5. apply left to right

    auto split = ComplexSplitStateThreeSiteThetaLeftToRight(Theta_gate, Dl, d, Dr, chi, cutoff);

    // 6. write back

    psi.A[j] = split.A0;
    psi.A[j + 1] = split.A1;
    psi.A[j + 2] = split.A2;

    // shit orthogonality center
    psi.center = j + 2;

    return split;
}

void ComplexApplyThreeSiteLayer(
    CMPS &psi,
    const ComplexMatrix_col &G3,
    int chi,
    int start)
{
    if (start < 0 || start > 2)
    {
        throw std::runtime_error("ComplexApplyThreeSiteLayer: start must be 0, 1, or 2");
    }

    // update
    for (int j = start; j + 2 < psi.N; j += 3)
    {
        ComplexUpdateThreeSite(psi, j, G3, chi);
    }
}

// Report truncation error per time step
struct FFDStepReport
{
    double max_discarded_fraction = 0.0;
    double total_discarded_weight = 0.0;
    int max_bond_dim = 1;
};
// Accumulation helpers
void accumulateSplitReport(
    FFDStepReport &report,
    const ComplexSplitStateResult &split)
{
    report.max_discarded_fraction = std::max(report.max_discarded_fraction, split.discarded_fraction);

    report.total_discarded_weight += split.discarded_weight;
}

void accumulateThreeSiteSplitReport(
    FFDStepReport &report,
    const ComplexThreeSiteSplitResult &split)
{
    report.max_discarded_fraction = std::max(report.max_discarded_fraction, std::max(split.discarded_fraction01, split.discarded_fraction12));

    report.total_discarded_weight += split.discarded_weight01 + split.discarded_weight12;
}

void accumulateStepReport(
    FFDStepReport &report,
    const FFDStepReport &other)
{
    report.max_discarded_fraction =
        std::max(report.max_discarded_fraction, other.max_discarded_fraction);

    report.total_discarded_weight += other.total_discarded_weight;

    report.max_bond_dim = std::max(report.max_bond_dim, other.max_bond_dim);
}

// helper for max bond dimension
void updateReportMaxBondDim(FFDStepReport &report,
                            const CMPS &psi)
{
    auto dims = ComplexBondDimensions(psi);
    for (int D : dims)
    {
        report.max_bond_dim = std::max(report.max_bond_dim, D);
    }
}

FFDStepReport ComplexApplyFFDThreeSiteBulkLayer(
    CMPS &psi,
    const std::vector<double> &b,
    double dt,
    int chi,
    int residue,
    double cutoff = -1.0)
{
    if (static_cast<int>(b.size()) != psi.N)
    {
        throw std::runtime_error("ComplexApplyFFDThreeSiteBulkLayer coefficient size mismatch ");
    }

    if (residue < 0 || residue > 2)
    {
        throw std::runtime_error("ComplexApplyFFDThreeSiteBulkLayer: residue must be a mod3");
    }

    FFDStepReport report;

    for (int m = residue; m < psi.N; m += 3)
    {
        if (m < 2)
        {
            continue;
        }

        int left = m - 2;
        auto G3 = makeFFDThreeSiteTermGate(b[m], dt);
        auto split = ComplexUpdateThreeSite(psi, left, G3, chi, cutoff);
        // update report
        accumulateThreeSiteSplitReport(report, split);
    }
    // update max bond dimension
    updateReportMaxBondDim(report, psi);

    return report;
}

FFDStepReport ComplexApplyFFDStep(
    CMPS &psi,
    const std::vector<double> &b,
    double dt,
    int chi,
    double cutoff = -1.0)
{
    if (static_cast<int>(b.size()) != psi.N)
    {
        throw std::runtime_error("ComplexApplyFFDStep  b size mismatch");
    }
    FFDStepReport report;

    // h_0, h_1, then bulk residues 2, 0, 1

    if (psi.N >= 1)
    {
        auto G0 = makeFFDOneSiteTermGate(b[0], dt);
        ComplexApplyOneSiteGate(psi, 0, G0);
    }

    if (psi.N >= 2)
    {
        auto G1 = makeFFDTwoSiteTermGate(b[1], dt);
        auto split01 = ComplexUpdateTwoSite(psi, 0, G1, chi, cutoff);
        accumulateSplitReport(report, split01);
    }

    // Apply bulk layers
    auto r2 = ComplexApplyFFDThreeSiteBulkLayer(psi, b, dt, chi, 2, cutoff);
    auto r0 = ComplexApplyFFDThreeSiteBulkLayer(psi, b, dt, chi, 0, cutoff);
    auto r1 = ComplexApplyFFDThreeSiteBulkLayer(psi, b, dt, chi, 1, cutoff);

    // update discarded fraction
    report.max_discarded_fraction = std::max(report.max_discarded_fraction, std::max(r2.max_discarded_fraction, std::max(r0.max_discarded_fraction, r1.max_discarded_fraction)));
    // update discarded weight
    report.total_discarded_weight += r2.total_discarded_weight + r0.total_discarded_weight + r1.total_discarded_weight;

    updateReportMaxBondDim(report, psi);
    return report;
}

FFDStepReport ComplexApplyFFDStepSecondOrder(
    CMPS &psi,
    const std::vector<double> &b,
    double dt,
    int chi,
    double cutoff = -1.0)
{

    if (static_cast<int>(b.size()) != psi.N)
    {
        throw std::runtime_error("ComplexApplyFFDStepSecondOrder: coupling size mismatch ");
    }

    FFDStepReport report;
    double half_dt = 0.5 * dt;

    // h0 half step

    if (psi.N >= 1)
    {
        auto G0 = makeFFDOneSiteTermGate(b[0], half_dt);
        ComplexApplyOneSiteGate(psi, 0, G0);
    }
    // h1 half step
    if (psi.N >= 2)
    {
        auto G1 = makeFFDTwoSiteTermGate(b[1], half_dt);
        auto split = ComplexUpdateTwoSite(psi, 0, G1, chi, cutoff);
        accumulateSplitReport(report, split);
    }

    // forward bulk half steps
    auto r2a = ComplexApplyFFDThreeSiteBulkLayer(psi, b, half_dt, chi, 2, cutoff);
    accumulateStepReport(report, r2a);

    auto r0a = ComplexApplyFFDThreeSiteBulkLayer(psi, b, half_dt, chi, 0, cutoff);
    accumulateStepReport(report, r0a);

    // center bulk full step
    auto r1 = ComplexApplyFFDThreeSiteBulkLayer(psi, b, dt, chi, 1, cutoff);
    accumulateStepReport(report, r1);

    // reverse bulk half steps
    auto r0b = ComplexApplyFFDThreeSiteBulkLayer(psi, b, half_dt, chi, 0, cutoff);
    accumulateStepReport(report, r0b);

    auto r2b = ComplexApplyFFDThreeSiteBulkLayer(psi, b, half_dt, chi, 2, cutoff);
    accumulateStepReport(report, r2b);

    // h1 half step
    if (psi.N >= 2)
    {
        auto G1 = makeFFDTwoSiteTermGate(b[1], half_dt);
        auto split = ComplexUpdateTwoSite(psi, 0, G1, chi, cutoff);
        accumulateSplitReport(report, split);
    }

    // h0 half step
    if (psi.N >= 1)
    {
        auto G0 = makeFFDOneSiteTermGate(b[0], half_dt);
        ComplexApplyOneSiteGate(psi, 0, G0);
    }
    updateReportMaxBondDim(report, psi);

    return report;
}

// <AB>= Tr(A^dagger B)
Complex complexThetaInnerProduct(
    const ComplexMatrix_col &A,
    const ComplexMatrix_col &B)
{
    if (A.rows != B.rows || A.cols != B.cols)
    {
        throw std::runtime_error("complexThetaInnerProduct: shape mismatch");
    }

    Complex val = 0.0;

    for (int col = 0; col < A.cols; ++col)
    {
        for (int row = 0; row < B.rows; ++row)
        {
            val += std::conj(A(row, col)) * B(row, col);
        }
    }

    return val;
}

// One-Site Expectation

Complex expectationOneSiteOperator(
    const CMPS &psi_in,
    int j,
    const ComplexMatrix_col &O1)
{
    if (j < 0 || j >= psi_in.N)
    {
        throw std::runtime_error("expectationOneStieOperator: unvalid site index");
    }

    CMPS psi = psi_in;

    // prepareCanonicalBlock(psi, j, 1);
    prepareCanonicalBlockFast(psi, j, 1);

    // Make a copy of the tensor to contract
    const auto &A = psi.A[j];

    Complex numerator = 0.0;
    Complex denominator = 0.0;

    for (int l = 0; l < A.Dl; ++l)
    {
        for (int r = 0; r < A.Dr; ++r)
        {
            for (int s = 0; s < A.d; ++s)
            {
                denominator += std::conj(A(l, s, r)) * A(l, s, r);

                for (int t = 0; t < A.d; ++t)
                {
                    numerator += std::conj(A(l, s, r)) * O1(s, t) * A(l, t, r);
                }
            }
        }
    }
    return numerator / denominator;
}

// Two-site expectation

Complex expectationTwoSiteOperator(
    const CMPS &psi_in,
    int j,
    const ComplexMatrix_col &O2)
{
    CMPS psi = psi_in;

    // prepareCanonicalBlock(psi, j, 2);
    prepareCanonicalBlockFast(psi, j, 2);

    auto Theta = ComplexBuildTwoSiteTheta(psi.A[j], psi.A[j + 1]);

    int Dl = psi.A[j].Dl;
    int d = psi.A[j].d;
    int Dr = psi.A[j + 1].Dr;

    auto OTheta = ComplexApplyTwoSiteGateToTheta(Dl, Dr, d, Theta, O2);

    return complexThetaInnerProduct(Theta, OTheta) / complexThetaInnerProduct(Theta, Theta);
}

// Three-site Expectation

Complex expectationThreeSiteOperator(
    const CMPS &psi_in,
    int j,
    const ComplexMatrix_col &O3)
{
    CMPS psi = psi_in;
    // prepareCanonicalBlock(psi, j, 3);
    prepareCanonicalBlockFast(psi, j, 3);

    auto Theta = ComplexBuildThreeSiteTheta(psi.A[j], psi.A[j + 1], psi.A[j + 2]);

    int Dl = psi.A[j].Dl;
    int d = psi.A[j].d;
    int Dr = psi.A[j + 2].Dr;

    auto OTheta = ComplexApplyThreeSiteGateToTheta(Dl, Dr, d, Theta, O3);

    return complexThetaInnerProduct(Theta, OTheta) / complexThetaInnerProduct(Theta, Theta);
}

// bond = N/2-1  // half-chain entropy

std::vector<double> schmidtValuesAcrossBond(const CMPS &psi_in, int bond)
{
    if (bond < 0 || bond >= psi_in.N - 1)
    {
        throw std::runtime_error("schmidtValuesAcrossBond: unvalid bond site");
    }

    CMPS psi = psi_in;

    // Center site = bond.
    // Reshape as (left block + physical site) |(right block)
    // |psi> = sum_{(l,s),r} M_{(l,s),r} |(l,s)_L> |r_R>

    moveOrthogonalityCenter(psi, bond);

    const auto &A = psi.A[bond];

    ComplexMatrix_col M(A.Dl * A.d, A.Dr);

    for (int l = 0; l < A.Dl; ++l)
    {
        for (int s = 0; s < A.d; ++s)
        {
            for (int r = 0; r < A.Dr; ++r)
            {
                int row = l * A.d + s;

                M(row, r) = A(l, s, r);
            }
        }
    }

    auto svd = complexSvd(M);
    return svd.S;
}

double entropyFromSchimdtValues(const std::vector<double> &lambda)
{
    double norm2 = 0.0;

    for (double x : lambda)
    {
        norm2 += x * x;
    }

    if (norm2 <= 0.0)
    {
        return 0.0;
    }

    double entropy = 0.0;

    for (double x : lambda)
    {
        double p = x * x / norm2;

        if (p > 0.0)
        {
            entropy += -p * std::log(p);
        }
    }

    return entropy;
}

// Bipartite Entropy Helper
// bipartiteEntropyFirstEllSites(psi, 10) will give
// S( sites 0..9 | sites 10..N-1 )

double bipartiteEntropyFirstEllSites(const CMPS &psi, int ell)
{
    if (ell <= 0 || ell > psi.N)
    {
        throw std::runtime_error("bipartiteEntropyFirstEllSites: invalid bipartition");
    }

    int bond = ell - 1;
    auto lambda = schmidtValuesAcrossBond(psi, bond);

    return entropyFromSchimdtValues(lambda);
}

std::vector<double> bipartiteEntropiesFirstEllSites(CMPS psi, const std::vector<int> &ells)
{
    std::vector<double> entropies;
    entropies.reserve(ells.size());

    for (int ell : ells)
    {
        if (ell <= 0 || ell >= psi.N)
        {
            throw std::runtime_error("invalid ell in  bipartite entropy sweep");
        }

        int bond = ell - 1;

        moveOrthogonalityCenter(psi, bond);

        const auto &A = psi.A[bond];

        ComplexMatrix_col M(A.Dl * A.d, A.Dr);

        for (int l = 0; l < A.Dl; ++l)
        {
            for (int s = 0; s < A.d; ++s)
            {
                for (int r = 0; r < A.Dr; ++r)
                {
                    int row = l * A.d + s;
                    M(row, r) = A(l, s, r);
                }
            }
        }
        auto svd = complexSvd(M);
        entropies.push_back(entropyFromSchimdtValues(svd.S));
    }

    return entropies;
}

// h_m expectation dispatcher
// expectationFFDTerm(psi, 0)  // <X_0>
// expectationFFDTerm(psi, 1)  // <Z_0 X_1>
// expectationFFDTerm(psi, 2)  // <Z_0 Z_1 X_2>

Complex expectationFFDTerm(const CMPS &psi, int m)
{
    if (m < 0 || m >= psi.N)
    {
        throw std::runtime_error("expectationFFDTerm: invalid site number");
    }

    if (m == 0)
    {
        return expectationOneSiteOperator(psi, 0, makeFFDOneSiteTermOperator());
    }

    if (m == 1)
    {
        return expectationTwoSiteOperator(psi, 0, makeFFDTwoSiteTermOperator());
    }

    return expectationThreeSiteOperator(psi, m - 2, makeFFDThreeSiteTermOperator());
}

// Hamiltonian density dispatcher

Complex expectationCoupledFFDTerm(const CMPS &psi, const std::vector<double> &b, int m)
{
    return b[m] * expectationFFDTerm(psi, m);
}

// ----------------- END TEBD HELPERS ---------------

// ----------------- SANITY CHECKS -----------------------

void runPauliAndGateSanityChecks(bool check)
{

    if (check)
    {
        std::cout << " Pauli matrices and Gates sanity check";
        auto X = makePauliX();
        auto Y = makePauliY();
        auto Z = makePauliZ();

        auto ZX = complexKron(Z, X);
        auto ZY = complexKron(Z, Y);
        printComplexMatrix_col(ZX, "ZX");
        printComplexMatrix_col(ZY, "ZY");

        // Three site strings
        auto ZZX = complexKron(complexKron(Z, Z), X);
        auto XZZ = complexKron(complexKron(X, Z), Z);
        // printComplexMatrix_col(ZZX, "ZZX");
        // printComplexMatrix_col(XZZ, "XZZ");

        ComplexMatrix_col h3 = makeFFD3SiteHamiltonian(1.0, 1.0);
        printComplexMatrix_col(h3, "3 site bond Hamiltonian h3");

        // Paulis hould square to Identity
        auto I8 = makeComplexIdentity(8);
        auto ZZX2 = cmatmul(ZZX, ZZX);
        auto XZZ2 = cmatmul(XZZ, XZZ);
        auto comm = cmatmul(ZZX, XZZ);
        auto comm_reverse = cmatmul(XZZ, ZZX);

        std::cout << "ZZX^2 error: " << maxAbsDiffComplex(ZZX2, I8) << "\n";
        std::cout << "XZZ^2 error: " << maxAbsDiffComplex(XZZ2, I8) << "\n";
        std::cout << "[ZZX,XZZ] erro: " << maxAbsDiffComplex(comm, comm_reverse) << "\n";
    }
}

void runThreeSiteSplitSanityCheck(bool check)
{
    if (check)
    {
        std::cout << " Three site split sanity check";
        auto I8 = makeComplexIdentity(8);
        double dt = 0.01;
        auto G3 = makeFFD3SiteGate(1.0, 1.0, dt);
        auto G3_unitarity = cmatmul(G3, conjTranspose(G3));
        int chi = 4;

        std::cout << "G3 unitary error: "
                  << maxAbsDiffComplex(G3_unitarity, I8)
                  << "\n";

        // printComplexMatrix_col(G3, "G3");

        // 1. Build  a 3-site complex product state |101>

        // row = l * d + s0 = s0
        // col = 2*s1 + s2
        // col 0 -> (s1,s2) = (0,0)
        // col 1 -> (s1,s2) = (0,1)
        // col 2 -> (s1,s2) = (1,0)
        // col 3 -> (s1,s2) = (1,1)
        // we expect:
        // |101> Theta(2x4)
        // 0 0 0 0
        // 0 1 0 0

        auto ThreeSite_c = makeComplexProductState({1, 0, 1});

        auto Theta_ThreeSite_c = ComplexBuildThreeSiteTheta(
            ThreeSite_c.A[0],
            ThreeSite_c.A[1],
            ThreeSite_c.A[2]);

        printComplexMatrix_col(Theta_ThreeSite_c, "|101> Theta");

        // Apply Gate check (Identity)
        auto Theta_id = ComplexApplyThreeSiteGateToTheta(1, 1, 2, Theta_ThreeSite_c, I8);
        printComplexMatrix_col(Theta_id, "Theta after I8");
        std::cout << "Identity 3-site error: "
                  << maxAbsDiffComplex(Theta_id, Theta_ThreeSite_c) << "\n";
        auto Theta_dt = ComplexApplyThreeSiteGateToTheta(1, 1, 2, Theta_ThreeSite_c, G3);
        printComplexMatrix_col(Theta_dt, "Updated Theta");

        auto split3 = ComplexSplitStateThreeSiteThetaLeftToRight(Theta_dt,
                                                                 ThreeSite_c.A[0].Dl,
                                                                 ThreeSite_c.A[0].d,
                                                                 ThreeSite_c.A[2].Dr,
                                                                 chi);

        auto Theta_rebuilt = ComplexBuildThreeSiteTheta(
            split3.A0,
            split3.A1,
            split3.A2);
        printComplexMatrix_col(Theta_rebuilt, "Theta rebuilt after 2 SVDs");

        std::cout << "3-site split/rebuild error: "
                  << maxAbsDiffComplex(Theta_rebuilt, Theta_dt)
                  << "\n";

        std::cout << "chi01 = " << split3.schmidt01.size()
                  << ", chi12 = " << split3.schmidt12.size()
                  << "\n";

        std::cout << "discarded weight 01 = " << split3.discarded_weight01
                  << ", discarded weight 12 = " << split3.discarded_weight12
                  << "\n";
    }
}

void runThreeSiteUpdateSanityCheck(bool check)
{
    if (!check)
    {
        return;
    }
    std::cout << "Three Site update Sanity Check ";
    double dt = 0.01;
    int chi = 4;
    auto G3 = makeFFD3SiteGate(1.0, 1.0, dt);

    // 1. Build  a 3-site complex product state |101>

    auto ThreeSite_c = makeComplexProductState({1, 0, 1});

    auto Theta_ThreeSite_c = ComplexBuildThreeSiteTheta(
        ThreeSite_c.A[0],
        ThreeSite_c.A[1],
        ThreeSite_c.A[2]);

    // row = l * d + s0 = s0
    // col = 2*s1 + s2
    // col 0 -> (s1,s2) = (0,0)
    // col 1 -> (s1,s2) = (0,1)
    // col 2 -> (s1,s2) = (1,0)
    // col 3 -> (s1,s2) = (1,1)
    // we expect:
    // |101> Theta(2x4)
    // 0 0 0 0
    // 0 1 0 0

    printComplexMatrix_col(Theta_ThreeSite_c, "|101> Theta");

    // Apply Gate

    auto Theta_dt = ComplexApplyThreeSiteGateToTheta(1, 1, 2, Theta_ThreeSite_c, G3);
    printComplexMatrix_col(Theta_dt, "Updated Theta");

    auto update_report = ComplexUpdateThreeSite(ThreeSite_c,
                                                0, G3, chi);

    auto Theta_after_update = ComplexBuildThreeSiteTheta(
        ThreeSite_c.A[0],
        ThreeSite_c.A[1],
        ThreeSite_c.A[2]);

    std::cout << "CMPS update/rebuild error: "
              << maxAbsDiffComplex(Theta_after_update, Theta_dt) << "\n";

    std::cout << "CMPS update chi01 = " << update_report.schmidt01.size()
              << ", chi12 = " << update_report.schmidt12.size()
              << "\n";
}

void runThreeSiteLayerSanityCheck(bool check)
{
    if (check)
    {
        std::cout << "Three site layer Sanity Check ";
        double dt = 0.01;
        auto G3 = makeFFD3SiteGate(1.0, 1.0, dt);
        int chi = 4;

        // sweep check
        // CMPS layer_a;
        // layer_a.N = 6;
        // layer_a.d = 2;

        // CMPS layer_b;
        // layer_b.N = 6;
        // layer_b.d = 2;
        // set |psi> = |101010>

        // for (int j = 0; j < 6; ++j)
        // {
        //     ComplexSiteTensor T(1, 2, 1);
        //     int s = (j % 2 == 0) ? 1 : 0;
        //     T(0, s, 0) = Complex{1.0, 0.0};

        //     layer_a.A.push_back(T);
        //     layer_b.A.push_back(T);
        // }
        auto layer_a = makeComplexProductState({1, 0, 1, 0, 1, 0});
        auto layer_b = makeComplexProductState({1, 0, 1, 0, 1, 0});

        // Apply two methods
        ComplexApplyThreeSiteLayer(layer_a, G3, chi, 0);
        ComplexUpdateThreeSite(layer_b, 0, G3, chi);
        ComplexUpdateThreeSite(layer_b, 3, G3, chi);

        // Compare results
        auto theta_a_012 = ComplexBuildThreeSiteTheta(layer_a.A[0], layer_a.A[1], layer_a.A[2]);

        auto theta_b_012 = ComplexBuildThreeSiteTheta(layer_b.A[0], layer_b.A[1], layer_b.A[2]);

        auto theta_a_345 = ComplexBuildThreeSiteTheta(layer_a.A[3], layer_a.A[4], layer_a.A[5]);

        auto theta_b_345 = ComplexBuildThreeSiteTheta(layer_b.A[3], layer_b.A[4], layer_b.A[5]);

        std::cout << "layer start=0 triple 012 error: "
                  << maxAbsDiffComplex(theta_a_012, theta_b_012)
                  << "\n";

        std::cout << "layer start=0 triple 345 error: "
                  << maxAbsDiffComplex(theta_a_345, theta_b_345)
                  << "\n";
    }
}

void runFendleyGateSanityCheck(bool check)
{
    if (!check)
    {
        return;
    }

    std::cout << "Fendley Gate Sanity Check \n";

    double dt = 0.01;
    double b = 2.0;

    auto G1 = makeFFDOneSiteTermGate(b, dt);
    auto G2 = makeFFDTwoSiteTermGate(b, dt);
    auto G3 = makeFFDThreeSiteTermGate(b, dt);

    auto I2 = makeComplexIdentity(2);
    auto I4 = makeComplexIdentity(4);
    auto I8 = makeComplexIdentity(8);

    std::cout << "FFD one-site gate unitary error: "
              << maxAbsDiffComplex(cmatmul(G1, conjTranspose(G1)), I2)
              << "\n";

    std::cout << "FFD two-site gate unitary error: "
              << maxAbsDiffComplex(cmatmul(G2, conjTranspose(G2)), I4)
              << "\n";

    std::cout << "FFD three-site gate unitary error: "
              << maxAbsDiffComplex(cmatmul(G3, conjTranspose(G3)), I8)
              << "\n";
}

void runFendleyBulkLayerSanityCheck(bool check)
{
    if (!check)
    {
        return;
    }

    std::cout << "Fendley Bulk Layer Sanity Check \n";
    auto psi_a = makeComplexProductState({1, 0, 1, 0, 1, 0});
    auto psi_b = makeComplexProductState({1, 0, 1, 0, 1, 0});

    std::vector<double> b = {1.0, 2.0, 3.0, 1.0, 2.0, 3.0};
    double dt = 0.01;
    int chi = 4;

    ComplexApplyFFDThreeSiteBulkLayer(psi_a, b, dt, chi, 2);

    auto G2 = makeFFDThreeSiteTermGate(b[2], dt);
    auto G5 = makeFFDThreeSiteTermGate(b[5], dt);

    ComplexUpdateThreeSite(psi_b, 0, G2, chi);
    ComplexUpdateThreeSite(psi_b, 3, G5, chi);

    // Compare results
    auto theta_a_012 = ComplexBuildThreeSiteTheta(psi_a.A[0], psi_a.A[1], psi_a.A[2]);

    auto theta_b_012 = ComplexBuildThreeSiteTheta(psi_b.A[0], psi_b.A[1], psi_b.A[2]);

    auto theta_a_345 = ComplexBuildThreeSiteTheta(psi_a.A[3], psi_a.A[4], psi_a.A[5]);

    auto theta_b_345 = ComplexBuildThreeSiteTheta(psi_b.A[3], psi_b.A[4], psi_b.A[5]);

    std::cout << "FFD layer start=2 triple 012 error: "
              << maxAbsDiffComplex(theta_a_012, theta_b_012)
              << "\n";

    std::cout << "FFD layer start=2 triple 345 error: "
              << maxAbsDiffComplex(theta_a_345, theta_b_345)
              << "\n";
}

void runOneSiteGateSanityCheck(bool check)
{
    if (!check)
    {
        return;
    }

    std::cout << "run One Site Gate Sanity Check\n";

    auto psi = makeComplexProductState({1});
    double dt = 0.01;
    double b0 = 2.0;

    auto G1 = makeFFDOneSiteTermGate(b0, dt);

    ComplexApplyOneSiteGate(psi, 0, G1);

    ComplexMatrix_col v(2, 1);
    v(1, 0) = Complex{1.0, 0.0};
    // manual multiplication
    auto expected = cmatmul(G1, v);

    ComplexMatrix_col got(2, 1);
    got(0, 0) = psi.A[0](0, 0, 0);
    got(1, 0) = psi.A[0](0, 1, 0);

    // print result
    std::cout << "one-site gate update error: "
              << maxAbsDiffComplex(got, expected)
              << "\n";
}

void runTwoSiteUpdateSanityCheck(bool check)
{
    if (!check)
    {
        return;
    }

    std::cout << "Two-site update sanity check\n";

    double dt = 0.01;
    double b1 = 2.0;
    int chi = 4;
    double cutoff = 1e-10;

    auto psi = makeComplexProductState({1, 0});
    auto G2 = makeFFDTwoSiteTermGate(b1, dt);

    auto Theta = ComplexBuildTwoSiteTheta(psi.A[0], psi.A[1]);
    auto Theta_expected = ComplexApplyTwoSiteGateToTheta(1, 1, 2, Theta, G2);

    // truncated MPS
    auto report = ComplexUpdateTwoSite(psi, 0, G2, chi, cutoff);

    printComplexSplitReport(report, "two-site split");

    auto Theta_after = ComplexBuildTwoSiteTheta(psi.A[0], psi.A[1]);

    std::cout << "two-site update/rebuild error: "
              << maxAbsDiffComplex(Theta_after, Theta_expected)
              << "\n";
}

void runFFDStepSmallChainSanityCheck(bool check)
{
    if (!check)
    {
        return;
    }

    std::vector<double> b = {1.0, 2.0, 3.0};
    double dt = 0.01;
    int chi = 4;

    auto psi = makeComplexProductState({1, 0, 1});
    auto psi_test = makeComplexProductState({1, 0, 1});

    ComplexApplyFFDStep(psi, b, dt, chi);

    auto G0 = makeFFDOneSiteTermGate(b[0], dt);
    auto G1 = makeFFDTwoSiteTermGate(b[1], dt);

    ComplexApplyOneSiteGate(psi_test, 0, G0);
    ComplexUpdateTwoSite(psi_test, 0, G1, chi);

    // Apply bulk layers
    ComplexApplyFFDThreeSiteBulkLayer(psi_test, b, dt, chi, 2);
    ComplexApplyFFDThreeSiteBulkLayer(psi_test, b, dt, chi, 0);
    ComplexApplyFFDThreeSiteBulkLayer(psi_test, b, dt, chi, 1);

    auto Theta_layer = ComplexBuildThreeSiteTheta(psi.A[0], psi.A[1], psi.A[2]);
    auto Theta_test_layer = ComplexBuildThreeSiteTheta(psi_test.A[0], psi_test.A[1], psi_test.A[2]);

    std::cout << "FFD step small-chain error:"
              << maxAbsDiffComplex(Theta_layer, Theta_test_layer)
              << "\n";
}

void runInputStateSanityCheck(bool check)
{
    if (!check)
    {
        return;
    }
    double theta = std::acos(-1.0) / 8.0;
    auto psi = makeTiltedProductState(2, theta);
    std::cout << "tilted site amp |0>" << psi.A[0](0, 0, 0)
              << "expected" << std::sin(theta) << "\n";

    std::cout << "tilted site amp |1>" << psi.A[0](0, 1, 0)
              << "expected" << std::cos(theta) << "\n";

    auto b = makeStaggeredSqrtCouplings(6, 1.0, 2.0, 3.0);
    printDoubleVector(b, "sqrt of staggered couplings");
}

void runCanonicalizeSiteCheck(bool check, bool LR)
{
    if (!check)
    {
        return;
    }

    auto psi = makeTiltedProductState(4, std::acos(-1.0) / 8.0);
    ComplexApplyFFDStep(psi, makeStaggeredSqrtCouplings(4, 1, 2, 3), 0.01, 16, -1.0);
    if (LR)
    {
        leftCanonicalizeSite(psi, 0);

        // Left isometry check
        ComplexMatrix_col gram(psi.A[0].Dr, psi.A[0].Dr);

        for (int a = 0; a < psi.A[0].Dr; ++a)
        {
            for (int b = 0; b < psi.A[0].Dr; ++b)
            {
                Complex val = 0.0;
                for (int l = 0; l < psi.A[0].Dl; ++l)
                {
                    for (int s = 0; s < psi.A[0].d; ++s)
                    {
                        val += std::conj(psi.A[0](l, s, a)) * psi.A[0](l, s, b);
                    }
                }

                gram(a, b) = val;
            }
        }

        auto I = makeComplexIdentity(psi.A[0].Dr);

        std::cout << "left canonical site 0 error: "
                  << maxAbsDiffComplex(gram, I)
                  << "\n";
    }
    else
    {
        // Right isometry check
        rightCanonicalizeSite(psi, 1);

        ComplexMatrix_col gram(psi.A[1].Dl, psi.A[1].Dl);

        for (int a = 0; a < psi.A[1].Dl; ++a)
        {
            for (int b = 0; b < psi.A[1].Dl; ++b)
            {
                Complex val = 0.0;
                for (int r = 0; r < psi.A[1].Dr; ++r)
                {
                    for (int s = 0; s < psi.A[1].d; ++s)
                    {
                        val += psi.A[1](a, s, r) * std::conj(psi.A[1](b, s, r));
                    }
                }

                gram(a, b) = val;
            }
        }

        auto I = makeComplexIdentity(psi.A[1].Dl);

        std::cout << "right canonical site 1 error: "
                  << maxAbsDiffComplex(gram, I)
                  << "\n";
    }
}

void runBlockLevelCanonicalizationCheck(bool check, int order)
{

    if (!check)
    {
        return;
    }

    int N = 5;
    double dt = 0.01;
    int nsteps = 10;
    int chi_max = 32;
    double cutoff = 1e-10;
    double alpha = 1.0;
    double beta = 2.0;
    double gamma = 3.0;

    auto b = makeStaggeredSqrtCouplings(N, alpha, beta, gamma);

    // Initial Tilted state
    double theta = M_PI / 8.0;
    auto psi = makeTiltedProductState(N, theta);

    // loop over steps
    for (int step = 0; step <= nsteps; ++step)
    {
        if (step < nsteps)
        {
            if (order < 2)
            {
                ComplexApplyFFDStep(psi, b, dt, chi_max, cutoff);
            }
            else
            {
                ComplexApplyFFDStepSecondOrder(psi, b, dt, chi_max, cutoff);
            }
        }
    }

    prepareCanonicalBlock(psi, 1, 3);

    // Left isometry check
    ComplexMatrix_col gram_L(psi.A[0].Dr, psi.A[0].Dr);

    for (int a = 0; a < psi.A[0].Dr; ++a)
    {
        for (int b = 0; b < psi.A[0].Dr; ++b)
        {
            Complex val = 0.0;
            for (int l = 0; l < psi.A[0].Dl; ++l)
            {
                for (int s = 0; s < psi.A[0].d; ++s)
                {
                    val += std::conj(psi.A[0](l, s, a)) * psi.A[0](l, s, b);
                }
            }

            gram_L(a, b) = val;
        }
    }

    auto I_L = makeComplexIdentity(psi.A[0].Dr);

    std::cout << "block prep left canonical site 0 error: "
              << maxAbsDiffComplex(gram_L, I_L)
              << "\n";

    ComplexMatrix_col gram_R(psi.A[4].Dl, psi.A[4].Dl);

    for (int a = 0; a < psi.A[4].Dl; ++a)
    {
        for (int b = 0; b < psi.A[4].Dl; ++b)
        {
            Complex val = 0.0;
            for (int r = 0; r < psi.A[4].Dr; ++r)
            {
                for (int s = 0; s < psi.A[4].d; ++s)
                {
                    val += psi.A[4](a, s, r) * std::conj(psi.A[4](b, s, r));
                }
            }

            gram_R(a, b) = val;
        }
    }

    auto I_R = makeComplexIdentity(psi.A[4].Dl);

    std::cout << "block prep right canonical site 4 error: "
              << maxAbsDiffComplex(gram_R, I_R)
              << "\n";
}

void runExpectationValueCheck(bool check)
{
    if (!check)
    {
        return;
    }

    int N = 5;
    double theta = M_PI / 8.0;

    auto psi = makeTiltedProductState(N, theta);

    auto h_0 = expectationFFDTerm(psi, 0);
    auto h_1 = expectationFFDTerm(psi, 1);
    auto h_2 = expectationFFDTerm(psi, 2);

    std::cout << "h_0 =" << h_0 << "\n"
              << "h_1 =" << h_1 << "\n"
              << "h_2 =" << h_2 << "\n";
}

void runBipartiteEntropySanityCheck(bool check)
{

    if (!check)
    {
        return;
    }

    auto product = makeTiltedProductState(6, M_PI / 8.0);

    std::cout << "product S ell=3: "
              << bipartiteEntropyFirstEllSites(product, 3)
              << "\n";

    auto bell = makeComplexProductState({0, 0});

    double inv_sq = 1.0 / std::sqrt(2.0);
    bell.A[0] = ComplexSiteTensor(1, 2, 2);
    bell.A[1] = ComplexSiteTensor(2, 2, 1);

    bell.A[0](0, 0, 0) = inv_sq;
    bell.A[0](0, 1, 1) = inv_sq;
    bell.A[1](0, 0, 0) = 1.0;
    bell.A[1](1, 1, 0) = 1.0;

    bell.center = 0;

    std::cout << "Bell S ell =1"
              << bipartiteEntropyFirstEllSites(bell, 1)
              << " expected " << std::log(2.0)
              << "\n";
}

// ----------- DEMO FUNCTIONS -------------

void runFFDShortTimeEvolutionDemo(bool check, int order)
{
    if (!check)
    {
        return;
    }

    int N = 12;
    double dt = 0.01;
    int nsteps = 10;
    int chi_max = 16;
    double cutoff = 1e-10;
    double alpha = 1.0;
    double beta = 2.0;
    double gamma = 3.0;

    auto b = makeStaggeredSqrtCouplings(N, alpha, beta, gamma);

    // Initial state: Neel

    // std::vector<int> config(N);
    // for (int j = 0; j < N; ++j)
    // {
    //     config[j] = (j % 2 == 0) ? 1 : 0;
    // }
    // auto psi = makeComplexProductState(config);

    // Initial Tilted state
    double theta = M_PI / 8.0;
    auto psi = makeTiltedProductState(N, theta);

    // loop over steps
    for (int step = 0; step <= nsteps; ++step)
    {
        std::cout << "before step" << step << " ";
        printIntVector(ComplexBondDimensions(psi), "bond dims");
        // Bond dimension: the Schmidt rank across a cut. For a spin-1/2 chain, the exact maximum rank across a cut is:
        //  min(2^number_of_sites_left, 2^number_of_sites_right)

        if (step < nsteps)
        {
            if (order < 2)
            {
                auto report = ComplexApplyFFDStep(psi, b, dt, chi_max, cutoff);
                // print block decimaiton update
                std::cout << "after step" << step + 1
                          << " max discarded fraction: "
                          << report.max_discarded_fraction
                          << " total discarded weight: "
                          << report.total_discarded_weight
                          << " max bond dim: "
                          << report.max_bond_dim
                          << "\n";
            }
            else
            {
                auto report = ComplexApplyFFDStepSecondOrder(psi, b, dt, chi_max, cutoff);
                // print block decimaiton update
                std::cout << "after step" << step + 1
                          << " max discarded fraction: "
                          << report.max_discarded_fraction
                          << " total discarded weight: "
                          << report.total_discarded_weight
                          << " max bond dim: "
                          << report.max_bond_dim
                          << "\n";
            }
        }
    }
}

// <h_m(t)> simulation

void runFFDObservableTimeSeriesDemo(bool check, int order)
{
    if (!check)
    {
        return;
    }

    bool write_hm = false;
    bool write_entropy = true;

    int hm_measure_every = 5;
    int entropy_measure_every = 5;

    // Site number
    int N = 40;
    double dt = 0.1;
    // double dt = 0.02;
    // int n_steps = 100; // t_max= 2
    int n_steps = 50; // t_max =5
    int chi_max = 128;
    double cutoff = 1e-10;

    double theta = M_PI / 8.0;
    double alpha = 1.0;
    double beta = 2.0;
    double gamma = 3.0;

    if (!write_hm && !write_entropy)
    {
        throw std::runtime_error("No output spected: write_hm and write_entropy are both false");
    }

    auto b = makeStaggeredSqrtCouplings(N, alpha, beta, gamma);
    // Initial Tilted state
    auto psi = makeTiltedProductState(N, theta);

    // probes

    std::vector<int> probes = {
        N / 2 - 2,
        N / 2 - 1,
        N / 2,
        N / 2 + 1,
        N / 2 + 2};

    // entanglement cuts

    // std::vector<int> ells;
    // int ell_max = N / 2;
    // int n_ell = 10;

    // for (int k = 1; k <= n_ell; ++k)
    // {
    //     int ell = k * ell_max / n_ell;
    //     ells.push_back(ell);
    // }

    std::vector<int> ells = {5, 10, 15, 20};

    int dt_milli = static_cast<int>(std::round(1000.0 * dt));

    std::string mode_label = (write_hm && write_entropy) ? "hm_entropy" : write_hm ? "hm"
                                                                                   : "entropy";

    std::string run_label =
        "fast_" + mode_label +
        "_N" + std::to_string(N) +
        "_chi" + std::to_string(chi_max) +
        "_steps" + std::to_string(n_steps) +
        "_dtm" + std::to_string(dt_milli) +
        "_order" + std::to_string(order) +
        "_hmeas" + std::to_string(hm_measure_every) +
        "_Smeas" + std::to_string(entropy_measure_every);

    std::string filename =
        "tutorials/ffd_manual/ffd_hm_" + run_label + ".csv";

    std::string entropy_filename = "tutorials/ffd_manual/ffd_entropy_" + run_label + ".csv";

    std::ofstream out;
    std::ofstream entropy_out;

    if (write_hm)
    {
        if (std::ifstream(filename))
        {
            throw std::runtime_error("Output file already exists: " + filename);
        }

        out.open(filename);

        if (!out)
        {
            throw std::runtime_error("Could not open output file");
        }

        out << "step,t,m,h_real,h_imag,bh_real,bh_imag,"
            << "max_discarded_fraction,total_discarded_weight,max_bond_dim\n";
    }

    if (write_entropy)
    {
        if (std::ifstream(entropy_filename))
        {
            throw std::runtime_error("Output file already exists: " + entropy_filename);
        }

        entropy_out.open(entropy_filename);

        if (!entropy_out)
        {
            throw std::runtime_error("Could not open entropy output file");
        }

        entropy_out << "step,t,ell,S,max_discarded_fraction,total_discarded_weight,max_bond_dim\n";
    }

    FFDStepReport last_report;

    for (int step = 0; step <= n_steps; ++step)
    {
        double t = step * dt;

        if (step % 10 == 0 || step == n_steps)
        {
            std::cout << "step: " << step
                      << " / " << n_steps
                      << ", t = " << t
                      << "\n";
        }

        if (write_hm && (step % hm_measure_every == 0 || step == n_steps))
        {
            for (int m : probes)
            {
                auto hm = expectationFFDTerm(psi, m);
                auto bhm = expectationCoupledFFDTerm(psi, b, m);

                out << step << ","
                    << t << ","
                    << m << ","
                    << hm.real() << ","
                    << hm.imag() << ","
                    << bhm.real() << ","
                    << bhm.imag() << ","
                    << last_report.max_discarded_fraction << ","
                    << last_report.total_discarded_weight << ","
                    << last_report.max_bond_dim << "\n";
            }
        }

        if (write_entropy && (step % entropy_measure_every == 0 || step == n_steps))
        {
            auto entropies = bipartiteEntropiesFirstEllSites(psi, ells);

            for (std::size_t i = 0; i < ells.size(); ++i)
            {
                entropy_out << step << ","
                            << t << ","
                            << ells[i] << ","
                            << entropies[i] << ","
                            << last_report.max_discarded_fraction << ","
                            << last_report.total_discarded_weight << ","
                            << last_report.max_bond_dim << "\n";
            }
        }

        if (step < n_steps)
        {
            if (order < 2)
            {
                last_report = ComplexApplyFFDStep(psi, b, dt, chi_max, cutoff);
            }
            else
            {
                last_report = ComplexApplyFFDStepSecondOrder(psi, b, dt, chi_max, cutoff);
            }
        }
    }
    if (write_hm)
    {
        std::cout << "saved observable data to " << filename << "\n";
    }

    if (write_entropy)
    {
        std::cout << "saved entropy data to " << entropy_filename << "\n";
    }
}

int main()
{
    std::cout << "FFD manual TEBD tutorial placeholder.\n";
    bool PauliAndGate = false;
    bool ThreeSiteSplit = false;
    bool ThreeSiteUpdate = false;
    bool ThreeSiteLayer = false;
    bool FendleyGates = false;
    bool FendleyBulkLayer = false;
    bool OneSiteGate = false;
    bool TwoSiteUpdate = false;
    bool FFDStepSmallChain = false;
    bool FFDShortTimeEvolution = false;
    int Order = 2;
    bool InputStateSanityCheck = false;
    bool CanonicalizeSite = false;
    bool LeftRight = false;
    bool BlockLevelCanonicalization = false;
    bool ExpectationValue = false;
    bool FFDObservableTimeSeries = true;

    runPauliAndGateSanityChecks(PauliAndGate);
    runThreeSiteSplitSanityCheck(ThreeSiteSplit);
    runThreeSiteUpdateSanityCheck(ThreeSiteUpdate);
    runThreeSiteLayerSanityCheck(ThreeSiteLayer);
    runFendleyGateSanityCheck(FendleyGates);
    runFendleyBulkLayerSanityCheck(FendleyBulkLayer);
    runOneSiteGateSanityCheck(OneSiteGate);
    runTwoSiteUpdateSanityCheck(TwoSiteUpdate);
    runFFDStepSmallChainSanityCheck(FFDStepSmallChain);
    runFFDShortTimeEvolutionDemo(FFDShortTimeEvolution, Order);
    runInputStateSanityCheck(InputStateSanityCheck);
    runCanonicalizeSiteCheck(CanonicalizeSite, LeftRight);
    runCanonicalizeSiteCheck(CanonicalizeSite, true);
    runBlockLevelCanonicalizationCheck(BlockLevelCanonicalization, Order);
    runExpectationValueCheck(ExpectationValue);
    runFFDObservableTimeSeriesDemo(FFDObservableTimeSeries, Order);

    return 0;
}
