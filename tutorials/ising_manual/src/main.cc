#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#ifndef ACCELERATE_NEW_LAPACK
#define ACCELERATE_NEW_LAPACK
#endif
#include <Accelerate/Accelerate.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <string>
#include <vector>

using Complex = std::complex<double>;

struct IsingParams
{
  int N = 24;
  double J = 1.0;
  double h = 1.0;
  int nsteps = 50;
  double dt = 0.01;
};

// Matrix structure, save data stacking columns
struct Matrix_col
{
  int rows = 0;
  int cols = 0;
  std::vector<double> data;

  Matrix_col() = default;

  Matrix_col(int r, int c)
      : rows(r), cols(c), data(r * c, 0.0)
  {
  }

  double &operator()(int i, int j)
  {
    if (i < 0 || i >= rows || j < 0 || j >= cols)
    {
      throw std::out_of_range("Matrix index out of range");
    }
    return data[j * rows + i];
  }

  double operator()(int i, int j) const
  {
    if (i < 0 || i >= rows || j < 0 || j >= cols)
    {
      throw std::out_of_range("Matrix index out of range");
    }
    return data[j * rows + i];
  }
};

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

struct SiteTensor
{
  int Dl = 1;
  int d = 2;
  int Dr = 1;
  std::vector<double> data;

  SiteTensor() = default;

  SiteTensor(int Dl_, int d_, int Dr_)
      : Dl(Dl_), d(d_), Dr(Dr_), data(Dl_ * d_ * Dr_, 0.0)
  {
  }
  // Dim= Dr Dl d = (Dr-1)dD_l +(d-1)D_L +(D_l-1)
  double &operator()(int l, int s, int r)
  {
    if (l < 0 || l >= Dl || s < 0 || s >= d || r < 0 || r >= Dr)
    {
      throw std::out_of_range("SiteTensor index out of range");
    }
    return data[(r * d + s) * Dl + l];
  }

  double operator()(int l, int s, int r) const
  {
    if (l < 0 || l >= Dl || s < 0 || s >= d || r < 0 || r >= Dr)
    {
      throw std::out_of_range("SiteTensor index out of range");
    }
    return data[(r * d + s) * Dl + l];
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

struct MPS
{
  int N = 0;
  // Spin 1/2 states
  int d = 2;
  std::vector<SiteTensor> A;
};

struct CMPS
{
  int N = 0;
  // Spin 1/2 states
  int d = 2;
  std::vector<ComplexSiteTensor> A;
};

MPS makeProductState(const std::vector<int> &local_states, int d = 2)
{
  if (d <= 0)
  {
    throw std::runtime_error("makeProductState expects positive local dimension");
  }
  if (local_states.empty())
  {
    throw std::runtime_error("makeProductState expects at least one site");
  }

  MPS psi;
  psi.N = static_cast<int>(local_states.size());
  psi.d = d;
  psi.A.reserve(psi.N);

  for (int s : local_states)
  {
    if (s < 0 || s >= d)
    {
      throw std::runtime_error("makeProductState found invalid local state");
    }
    SiteTensor Aj(1, d, 1);
    Aj(0, s, 0) = 1.0;
    psi.A.push_back(Aj);
  }

  return psi;
}

CMPS makeComplexProductState(const std::vector<int> &local_states, int d = 2)
{
  if (d <= 0)
  {
    throw std::runtime_error("makeProductState expects positive local dimension");
  }
  if (local_states.empty())
  {
    throw std::runtime_error("makeProductState expects at least one site");
  }

  CMPS psi;
  psi.N = static_cast<int>(local_states.size());
  psi.d = d;
  psi.A.reserve(psi.N);

  for (int s : local_states)
  {
    if (s < 0 || s >= d)
    {
      throw std::runtime_error("makeProductState found invalid local state");
    }
    ComplexSiteTensor Aj(1, d, 1);
    Aj(0, s, 0) = 1.0;
    psi.A.push_back(Aj);
  }

  return psi;
}

MPS makeNeelState(int N)
{
  if (N <= 0)
  {
    throw std::runtime_error("makeNeelState expects N > 0");
  }

  std::vector<int> local_states(N, 0);
  for (int j = 0; j < N; ++j)
  {
    local_states[j] = (j % 2 == 0) ? 0 : 1;
  }

  return makeProductState(local_states, 2);
}

CMPS makeComplexNeelState(int N)
{
  if (N <= 0)
  {
    throw std::runtime_error("makeNeelState expects N > 0");
  }

  std::vector<int> local_states(N, 0);
  for (int j = 0; j < N; ++j)
  {
    local_states[j] = (j % 2 == 0) ? 0 : 1;
  }

  return makeComplexProductState(local_states, 2);
}

MPS makeDomainWallState(int N, int cut)
{
  if (N <= 0)
  {
    throw std::runtime_error("makeDomainWallState expects N > 0");
  }
  if (cut < 0 || cut > N)
  {
    throw std::runtime_error("makeDomainWallState expects 0 <= cut <= N");
  }

  std::vector<int> local_states(N, 1);
  for (int j = 0; j < cut; ++j)
  {
    local_states[j] = 0;
  }

  return makeProductState(local_states, 2);
}

std::vector<int> productConfiguration(const MPS &psi)
{
  // answer to : Which local basis state is occupied at each site?

  if (static_cast<int>(psi.A.size()) != psi.N)
  {
    throw std::runtime_error("productConfiguration found inconsistent MPS length");
  }

  std::vector<int> config(psi.N, 0);
  for (int j = 0; j < psi.N; ++j)
  {
    auto const &Aj = psi.A.at(j);
    if (Aj.Dl != 1 || Aj.Dr != 1)
    {
      throw std::runtime_error("productConfiguration only supports bond-dimension-1 states");
    }

    int best_s = 0;
    double best_amp = std::abs(Aj(0, 0, 0));
    for (int s = 1; s < Aj.d; ++s)
    {
      double amp = std::abs(Aj(0, s, 0));
      if (amp > best_amp)
      {
        best_amp = amp;
        best_s = s;
      }
    }
    config[j] = best_s;
  }

  return config;
}

void printProductConfiguration(const MPS &psi, const std::string &name)
{
  auto config = productConfiguration(psi);
  std::cout << name << " (N=" << psi.N << "):";
  for (int s : config)
  {
    std::cout << " " << s;
  }
  std::cout << "\n";
}

Matrix_col add(const Matrix_col &A, const Matrix_col &B)
{

  if (A.rows != B.rows || A.cols != B.cols)
  {
    throw std::runtime_error(" dimension mismatch");
  }

  Matrix_col C(A.rows, A.cols);
  for (int i = 0; i < A.rows; ++i)
  {
    for (int j = 0; j < A.cols; ++j)
    {
      C(i, j) = A(i, j) + B(i, j);
    }
  }
  return C;
}

Matrix_col scale(double alpha, const Matrix_col &A)
{
  Matrix_col B(A.rows, A.cols);
  for (int i = 0; i < A.rows; ++i)
  {
    for (int j = 0; j < A.cols; ++j)
    {
      B(i, j) = alpha * A(i, j);
    }
  }
  return B;
}

// MATRIX TRANSPOSE
Matrix_col transpose(const Matrix_col &A)
{

  Matrix_col A_T(A.cols, A.rows);

  if (A.rows < 0 || A.cols < 0)
  {
    throw std::out_of_range("Matrix index out of range");
  }
  for (int i = 0; i < A.rows; ++i)
  {
    for (int j = 0; j < A.cols; ++j)
    {
      A_T(j, i) = A(i, j);
    }
  }

  return A_T;
}

// REAL Matrix eigen decomposition

struct EigenResult
{
  std::vector<double> evals;
  Matrix_col evecs;
};

EigenResult eighSymmetric(const Matrix_col &A)
{
  if (A.rows != A.cols)
  {
    throw std::runtime_error("eighSymmetric requires square matrices");
  }

  using LapackInt = int;

  LapackInt n = A.rows;
  LapackInt lda = n;
  LapackInt info = 0;

  Matrix_col Acopy = A;
  std::vector<double> w(n);

  char jobz = 'V';
  char uplo = 'U';

  double work_query = 0.0;
  LapackInt lwork = -1;

  // Ask for memory
  dsyev_(
      &jobz,
      &uplo,
      &n,
      Acopy.data.data(),
      &lda,
      w.data(),
      &work_query,
      &lwork,
      &info);

  if (info != 0)
  {
    throw std::runtime_error("dsyev workspace query failed");
  }

  // allocate temporary memory
  lwork = static_cast<LapackInt>(work_query);
  if (lwork < 1)
  {
    lwork = 1;
  }

  std::vector<double> work(lwork);

  dsyev_(
      &jobz,
      &uplo,
      &n,
      Acopy.data.data(),
      &lda,
      w.data(),
      work.data(),
      &lwork,
      &info);

  if (info != 0)
  {
    throw std::runtime_error("dysev failed");
  }

  return {w, Acopy};
}

// Reshave two sites state into a 2 by 2 matrix,
// |00>, |01>, |10>, |11>

Matrix_col reshapeTwoSitesState(const Matrix_col &psi_vec)
{
  if (psi_vec.rows != 4 || psi_vec.cols != 1)
  {
    throw std::runtime_error("reshapeTwoSitesState expects a 4x1 vector");
  }

  Matrix_col Theta(2, 2);

  for (int s1 = 0; s1 < 2; ++s1)
  {
    for (int s2 = 0; s2 < 2; ++s2)
    {
      int flat = 2 * s1 + s2;
      Theta(s1, s2) = psi_vec(flat, 0);
    }
  }

  return Theta;
}

// Matrix structure, save data stacking rows
// struct Matrix_row
// {
//   int rows = 0;
//   int cols = 0;
//   std::vector<double> data;

//   Matrix_row() = default;

//   Matrix_row(int r, int c)
//       : rows(r), cols(c), data(r * c, 0.0)
//   {
//   }

//   double &operator()(int i, int j)
//   {
//     if (i < 0 || i >= rows || j < 0 || j >= cols)
//     {
//       throw std::out_of_range("Matrix index out of range");
//     }
//     return data[i * cols + j];
//   }

//   double operator()(int i, int j) const
//   {
//     if (i < 0 || i >= rows || j < 0 || j >= cols)
//     {
//       throw std::out_of_range("Matrix index out of range");
//     }
//     return data[i * cols + j];
//   }
// };

void printMatrix_col(const Matrix_col &A, const std::string &name)
{
  std::cout << name << "(" << A.rows << "x" << A.cols
            << ")\n";
  for (int i = 0; i < A.rows; ++i)
  {
    for (int j = 0; j < A.cols; ++j)
    {
      std::cout << std::setw(10) << A(i, j) << " ";
    }
    std::cout << "\n";
  }
}

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

ComplexMatrix_col toComplexMatrix(const Matrix_col &A)
{
  ComplexMatrix_col out(A.rows, A.cols);
  for (int i = 0; i < A.rows; ++i)
  {
    for (int j = 0; j < A.cols; ++j)
    {
      out(i, j) = Complex{A(i, j), 0.0};
    }
  }
  return out;
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

// void printMatrix_row(const Matrix_row &A, const std::string &name)
// {
//   std::cout << name << "(" << A.rows << "x" << A.cols
//             << ")\n";
//   for (int i = 0; i < A.rows; ++i)
//   {
//     for (int j = 0; j < A.cols; ++j)
//     {
//       std::cout << std::setw(10) << A(i, j) << " ";
//     }
//     std::cout << "\n";
//   }
// }

// difference between row and col
void printVector(const std::vector<double> &v, const std::string &name)
{
  std::cout << name << ":";
  for (double x : v)
    std::cout << x << " ";
  std::cout << "\n";
}

void printIntVector(const std::vector<int> &v, const std::string &name)
{
  std::cout << name << ":";
  for (int x : v)
  {
    std::cout << " " << x;
  }
  std::cout << "\n";
}

// Matrix multiplication col
Matrix_col matmul(const Matrix_col &A, const Matrix_col &B)
{
  if (A.cols != B.rows)
  {
    throw std::runtime_error("matmul dimension mismatch");
  }

  Matrix_col C(A.rows, B.cols);

  // C = alpha * A * B + beta * C
  cblas_dgemm(
      CblasColMajor,
      CblasNoTrans,
      CblasNoTrans,
      A.rows, // m
      B.cols, // n
      A.cols, // k
      1.0,    // alpha
      A.data.data(),
      A.rows, // lda
      B.data.data(),
      B.rows, // ldb
      0.0,    // beta
      C.data.data(),
      C.rows // ldc
  );

  return C;
}

ComplexMatrix_col cmatmul(const ComplexMatrix_col &A, const ComplexMatrix_col &B)
{
  if (A.cols != B.rows)
  {
    throw std::runtime_error("cmatmul dimension mismatch");
  }

  ComplexMatrix_col C(A.rows, B.cols);

  // TODO(juan): replace this  with cblas_zgemm.
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

struct SVDResult
{
  Matrix_col U;
  std::vector<double> S;
  Matrix_col VT;
};

struct ComplexSVDResult
{
  ComplexMatrix_col U;
  std::vector<double> S;
  ComplexMatrix_col VH;
};

Matrix_col makeDiagonal(const std::vector<double> &s)
{
  int k = static_cast<int>(s.size());
  Matrix_col D(k, k);
  for (int i = 0; i < k; ++i)
  {
    D(i, i) = s[i];
  }
  return D;
}

ComplexMatrix_col makeComplexDiagonal(const std::vector<Complex> &s)
{
  int k = static_cast<int>(s.size());
  ComplexMatrix_col D(k, k);
  for (int i = 0; i < k; ++i)
  {
    D(i, i) = s[i];
  }
  return D;
}

ComplexMatrix_col makeComplexDiagonalFromReal(const std::vector<double> &s)
{
  int k = static_cast<int>(s.size());
  ComplexMatrix_col D(k, k);
  for (int i = 0; i < k; ++i)
  {
    D(i, i) = Complex{s[i], 0.0};
  }
  return D;
}
Matrix_col makeSqrtDiagonal(const std::vector<double> &s)
{
  int k = static_cast<int>(s.size());
  Matrix_col sqrtD(k, k);
  for (int i = 0; i < k; ++i)
  {
    sqrtD(i, i) = std::sqrt(s[i]);
  }

  return sqrtD;
}

double maxAbsDiff(const Matrix_col &A, const Matrix_col &B)
{
  if (A.rows != B.rows || A.cols != B.cols)
  {
    throw std::runtime_error("maxAbsDiff dimension mismatch");
  }
  // element-wise difference
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

double bondEntropyFromSchmidt(const std::vector<double> &schmidt)
{
  double norm2 = 0.0;
  for (double lambda : schmidt)
  {
    norm2 += lambda * lambda;
  }
  if (norm2 <= 0.0)
  {
    return 0.0;
  }

  double entropy = 0.0;
  for (double lambda : schmidt)
  {
    double p = (lambda * lambda) / norm2;
    if (p > 0.0)
    {
      entropy -= p * std::log(p);
    }
  }
  return entropy;
}

SVDResult svd(const Matrix_col &A)
{
  using LapackInt = int;

  LapackInt m = A.rows;
  LapackInt n = A.cols;
  LapackInt k = std::min(m, n);

  Matrix_col Acopy = A;
  Matrix_col U(m, k);
  Matrix_col VT(k, n);
  std::vector<double> S(k);

  char jobu = 'S';
  char jobvt = 'S';

  LapackInt lda = m;
  LapackInt ldu = m;
  LapackInt ldvt = k;
  LapackInt info = 0;

  double work_query = 0.0;
  LapackInt lwork = -1;

  // Ask for memory
  dgesvd_(
      &jobu,
      &jobvt,
      &m,
      &n,
      Acopy.data.data(),
      &lda,
      S.data(),
      U.data.data(),
      &ldu,
      VT.data.data(),
      &ldvt,
      &work_query,
      &lwork,
      &info);

  if (info != 0)
  {
    throw std::runtime_error("SVD workspace query failed");
  }

  // allocate temporary memory
  lwork = static_cast<LapackInt>(work_query);
  if (lwork < 1)
  {
    lwork = 1;
  }

  std::vector<double> work(lwork);

  dgesvd_(
      &jobu,
      &jobvt,
      &m,
      &n,
      Acopy.data.data(),
      &lda,
      S.data(),
      U.data.data(),
      &ldu,
      VT.data.data(),
      &ldvt,
      work.data(),
      &lwork,
      &info);

  if (info != 0)
  {
    throw std::runtime_error("SVD failed");
  }

  return {U, S, VT};
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

  //  real work space for zgesvd
  std::vector<double> rwork(5 * k);
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

SVDResult truncateSVD(const SVDResult &full, int chi)
{
  int kfull = static_cast<int>(full.S.size());

  if (chi < 0)
  {
    throw std::runtime_error("truncateSVD: chi must be positive");
  }
  // bond dimension
  int keep = std::min(chi, kfull);

  SVDResult out{
      Matrix_col(full.U.rows, keep),
      std::vector<double>(keep),
      Matrix_col(keep, full.VT.cols)};
  // fill truncated S vector
  for (int a = 0; a < keep; ++a)
  {
    out.S[a] = full.S[a];
  }

  //  fill truncated cols of U
  for (int i = 0; i < full.U.rows; ++i)
  {
    for (int a = 0; a < keep; ++a)
    {
      out.U(i, a) = full.U(i, a);
    }
  }

  //  fill truncated rows of VT
  for (int a = 0; a < keep; ++a)
  {
    for (int j = 0; j < full.VT.cols; ++j)
    {
      out.VT(a, j) = full.VT(a, j);
    }
  }

  return out;
}

// norm for imaginary time evolution

double schmidtNorm2(const std::vector<double> &s)
{
  double n2 = 0.0;
  for (double x : s)
  {
    n2 += x * x;
  }
  return n2;
}

// save discarded S entries sum

double discardedWeight(const SVDResult &full, int chi)
{
  int kfull = static_cast<int>(full.S.size());
  int keep = std::min(std::max(chi, 0), kfull);

  double w = 0.0;
  for (int a = keep; a < kfull; ++a)
  {
    w += full.S[a] * full.S[a];
  }
  return w;
}
//  Save discarded fraction
double discardedFraction(const SVDResult &full, int chi)
{
  double total = schmidtNorm2(full.S);
  if (total <= 0.0)
  {
    return 0.0;
  }
  return discardedWeight(full, chi) / total;
}

Matrix_col kron(const Matrix_col &A, const Matrix_col &B)
{
  int r = A.rows * B.rows;
  int c = A.cols * B.cols;
  Matrix_col kron_AB(r, c);
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
          kron_AB(row, col) = A(iA, jA) * B(iB, jB);
        }
      }
    }
  }
  return kron_AB;
}

// Complex tensor product

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

// Two site state update + truncation
// Bond Hamiltonian

Matrix_col makeBondHamiltonian(double J, double h)
{

  // one site Pauli Matrices
  Matrix_col I(2, 2);
  I(0, 0) = 1.0;
  I(1, 1) = 1.0;

  Matrix_col X(2, 2);
  X(0, 1) = 1.0;
  X(1, 0) = 1.0;

  Matrix_col Z(2, 2);
  Z(0, 0) = 1.0;
  Z(1, 1) = -1.0;

  // Tensor product check
  Matrix_col ZZ = kron(Z, Z);
  Matrix_col XI = kron(X, I);
  Matrix_col IX = kron(I, X);
  // 0.5 h to validate against exactly solvable Hamiltonian J=0

  Matrix_col h_bond = add(scale(-J, ZZ), add(scale(-0.5 * h, XI), scale(-0.5 * h, IX)));

  return h_bond;
}

Matrix_col makeBondHamiltonianForBond(int j, int N, double J, double h)
{

  // one site Pauli Matrices
  Matrix_col I(2, 2);
  I(0, 0) = 1.0;
  I(1, 1) = 1.0;

  Matrix_col X(2, 2);
  X(0, 1) = 1.0;
  X(1, 0) = 1.0;

  Matrix_col Z(2, 2);
  Z(0, 0) = 1.0;
  Z(1, 1) = -1.0;

  // Tensor product check
  Matrix_col ZZ = kron(Z, Z);
  Matrix_col XI = kron(X, I);
  Matrix_col IX = kron(I, X);

  // Normalize Hamiltonian
  double wL;
  double wR;
  if (j == 0)
  {
    wL = 1.0;
  }
  else
  {
    wL = 0.5;
  }
  if (j + 1 == N - 1)
  {
    wR = 1.0;
  }
  else
  {
    wR = 0.5;
  }

  Matrix_col h_bond = add(scale(-J, ZZ), add(scale(-wL * h, XI), scale(-wR * h, IX)));

  return h_bond;
}

// J=0 Hamiltonian
Matrix_col makeFullHamiltonianJ0(int N, double h)
{
  if (N <= 0)
  {
    throw std::runtime_error("makeFullHamiltonianJ0 expects N>0");
  }
  Matrix_col H(1 << N, 1 << N);

  Matrix_col I(2, 2);
  I(0, 0) = 1.0;
  I(1, 1) = 1.0;

  Matrix_col X(2, 2);
  X(0, 1) = 1.0;
  X(1, 0) = 1.0;

  for (int site = 0; site < N; ++site)
  {
    Matrix_col O_site(1, 1);
    O_site(0, 0) = 1.0;

    for (int i = 0; i < N; ++i)
    {
      // O = factor_0 ⊗ factor_1 ⊗ ... ⊗ factor_{N-1}
      O_site = kron(O_site, (i == site) ? X : I);
    }
    H = add(H, scale(-h, O_site));
  }

  return H;
}

Matrix_col makeImagTimeGate(const Matrix_col &h_bond, double dt)
{
  // bond exponentiation for imaginary time
  auto eig = eighSymmetric(h_bond);

  // DIAGONAL MATRIX EXPONENTIATION

  std::vector<double> exp_evals(eig.evals.size());

  for (int i = 0; i < static_cast<int>(eig.evals.size()); ++i)
  {
    exp_evals[i] = std::exp(-dt * eig.evals[i]);
  }
  Matrix_col ExpD = makeDiagonal(exp_evals);

  // Assemble the gate

  Matrix_col VT = transpose(eig.evecs);
  Matrix_col G = matmul(matmul(eig.evecs, ExpD), VT);

  return G;
}

ComplexMatrix_col makeRealTimeGate(const Matrix_col &h_bond, double dt)
{
  // bond exponentiation for real time
  auto eig = eighSymmetric(h_bond);

  // DIAGONAL MATRIX EXPONENTIATION

  std::vector<Complex> exp_evals(eig.evals.size());

  for (int i = 0; i < static_cast<int>(eig.evals.size()); ++i)
  {
    exp_evals[i] = std::exp(Complex{0.0, -dt * eig.evals[i]});
  }
  ComplexMatrix_col ExpD = makeComplexDiagonal(exp_evals);

  // Assemble the gate
  ComplexMatrix_col U = toComplexMatrix(eig.evecs);

  ComplexMatrix_col VT = toComplexMatrix(transpose(eig.evecs));
  ComplexMatrix_col G = cmatmul(cmatmul(U, ExpD), VT);

  // sanity check: real-time gate should be unitary
  ComplexMatrix_col Gdag = conjTranspose(G);
  ComplexMatrix_col check = cmatmul(G, Gdag);
  ComplexMatrix_col I_c = makeComplexIdentity(G.rows);
  if (maxAbsDiffComplex(check, I_c) > 1e-12)
  {
    throw std::runtime_error("makeRealTimeGate failed unitary check");
  }

  return G;
}

// Exact state evolution
ComplexMatrix_col evolveExactOneStep(const ComplexMatrix_col &psi,
                                     const Matrix_col &H,
                                     double dt)
{
  if (H.cols != psi.rows)
  {
    throw std::runtime_error("evolveExactOneStep Hamiltonian and state dim mismatch ");
  }
  ComplexMatrix_col U = makeRealTimeGate(H, dt);
  return cmatmul(U, psi);
}

std::vector<ComplexMatrix_col> makeRealTimeGates(int N, double J, double h, double dt)
{
  std::vector<ComplexMatrix_col> gates;
  for (int j = 0; j < N - 1; j++)
  {
    Matrix_col h_j = makeBondHamiltonianForBond(j, N, J, h);
    gates.push_back(makeRealTimeGate(h_j, dt));
  }
  return gates;
}

// Apply the two body gate

Matrix_col applyGate(const Matrix_col &G, const Matrix_col &psi)
{
  Matrix_col psi_new = matmul(G, psi);
  return psi_new;
}

ComplexMatrix_col ComplexApplyGate(const ComplexMatrix_col &G, const ComplexMatrix_col &psi)
{
  ComplexMatrix_col psi_new = cmatmul(G, psi);
  return psi_new;
}

// Split state for SVD
struct SplitStateResult
{
  Matrix_col left;
  Matrix_col right;
  std::vector<double> schmidt;
  double discarded_weight = 0.0;
  double discarded_fraction = 0.0;
};

SplitStateResult splitState(const Matrix_col &Theta, int chi, bool left_to_right)
{
  // Bond split \theta = U\sqrt{S_theta} \sqrt{S_theta}V^T
  auto svdTheta = svd(Theta);

  // truncate A and B
  // Truncated svd
  // Bond dimension
  auto svdTheta_trunc = truncateSVD(svdTheta, chi);

  double discarded_weight = discardedWeight(svdTheta, chi);
  double discarded_fraction = discardedFraction(svdTheta, chi);

  std::vector<double> S_trunc = svdTheta_trunc.S;
  Matrix_col S = makeDiagonal(S_trunc);
  // Symmetric splitting
  // Matrix_col Ssqrt_trunc = makeSqrtDiagonal(svdTheta_trunc.S);
  // Matrix_col left_trunc = matmul(svdTheta_trunc.U, Ssqrt_trunc);
  // Matrix_col right_trunc = matmul(Ssqrt_trunc, svdTheta_trunc.VT);
  Matrix_col left_trunc;
  Matrix_col right_trunc;

  if (left_to_right == true)
  {
    left_trunc = svdTheta_trunc.U;
    right_trunc = matmul(S, svdTheta_trunc.VT);
  }
  else
  {

    left_trunc = matmul(svdTheta_trunc.U, S);
    right_trunc = svdTheta_trunc.VT;
  }

  return {left_trunc, right_trunc, S_trunc, discarded_weight, discarded_fraction};
}

// Complex split state
struct ComplexSplitStateResult
{
  ComplexMatrix_col left;
  ComplexMatrix_col right;
  std::vector<double> schmidt;
  double discarded_weight = 0.0;
  double discarded_fraction = 0.0;
};

ComplexSplitStateResult ComplexsplitState(const ComplexMatrix_col &Theta, int chi, bool left_to_right)
{
  auto complexSvdTheta = complexSvd(Theta);
  int kfull = static_cast<int>(complexSvdTheta.S.size());
  if (chi < 0)
  {
    throw std::runtime_error("ComplexsplitState: chi must be nonnegative");
  }
  int keep = std::min(chi, kfull);

  ComplexMatrix_col U_trunc(complexSvdTheta.U.rows, keep);
  ComplexMatrix_col VH_trunc(keep, complexSvdTheta.VH.cols);
  std::vector<double> S_trunc(keep);

  for (int a = 0; a < keep; ++a)
  {
    S_trunc[a] = complexSvdTheta.S[a];
  }
  for (int i = 0; i < complexSvdTheta.U.rows; ++i)
  {
    for (int a = 0; a < keep; ++a)
    {
      U_trunc(i, a) = complexSvdTheta.U(i, a);
    }
  }
  for (int a = 0; a < keep; ++a)
  {
    for (int j = 0; j < complexSvdTheta.VH.cols; ++j)
    {
      VH_trunc(a, j) = complexSvdTheta.VH(a, j);
    }
  }

  double discarded_weight = 0.0;
  for (int a = keep; a < kfull; ++a)
  {
    discarded_weight += complexSvdTheta.S[a] * complexSvdTheta.S[a];
  }

  double total_weight = schmidtNorm2(complexSvdTheta.S);
  double discarded_fraction = (total_weight > 0.0) ? (discarded_weight / total_weight) : 0.0;

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

// TODO(COMPLEX-1): Implement ComplexWriteBondFromSplit just below this comment.
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

// TODO(COMPLEX-2): Implement ComplexUpdateBond after ComplexWriteBondFromSplit.
ComplexMatrix_col ComplexBuildBondTheta(const ComplexSiteTensor &Aj, const ComplexSiteTensor &Aj1)
// required bond dimension >= rank(Theta)
// bond dimension = Schmidt rank = number of nonzero singular values
// So the bond dimension : size of the internal summed index in the factorization
{
  // left side (l,s1)
  // right side (s2,r)

  ComplexMatrix_col Theta(Aj.Dl * Aj.d, Aj1.d * Aj1.Dr);

  // Check
  if (Aj.Dr != Aj1.Dl)
  {
    throw std::runtime_error("buildBondTheta bond mismatch");
  }

  int S1 = Aj.d;
  int S2 = Aj1.d;
  // sum index for Theta(l,s1,s2,r) = sum_m Aj(l,s1,m) * Aj1(m,s2,r)
  int Dm = Aj.Dr;

  for (int l = 0; l < Aj.Dl; ++l)
  {
    for (int r = 0; r < Aj1.Dr; ++r)
    {
      for (int s1 = 0; s1 < S1; ++s1)
      {
        for (int s2 = 0; s2 < S2; ++s2)
        {
          Complex val = 0.0;
          for (int m = 0; m < Dm; ++m)
          {
            val += Aj(l, s1, m) * Aj1(m, s2, r);
          }
          int row = l * S1 + s1;
          int col = s2 * Aj1.Dr + r;
          Theta(row, col) = val;
        }
      }
    }
  }
  return Theta;
}

ComplexMatrix_col ComplexApplyGateToTheta(int Dl, int Dr, int d, const ComplexMatrix_col &Theta, const ComplexMatrix_col &G)
{
  if (Theta.rows != Dl * d || Theta.cols != d * Dr)
  {
    throw std::runtime_error("applyGateToTheta shape mismatch");
  }
  if (G.rows != d * d || G.cols != d * d)
  {
    throw std::runtime_error("applyGateToTheta expects a (d*d)x(d*d) gate");
  }

  ComplexMatrix_col Theta_new(Dl * d, d * Dr);

  for (int l = 0; l < Dl; ++l)
  {
    for (int r = 0; r < Dr; ++r)
    {
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

      ComplexMatrix_col theta_lr_new = ComplexApplyGate(G, theta_lr);

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

ComplexSplitStateResult ComplexUpdateBond(CMPS &psi, int j, const ComplexMatrix_col &G, int chi, bool left_to_right)
{
  ComplexMatrix_col Theta = ComplexBuildBondTheta(psi.A[j], psi.A[j + 1]);
  int Dl = psi.A[j].Dl;
  int d = psi.A[j].d;
  int Dr = psi.A[j + 1].Dr;

  // SVD with chi bond dimension
  ComplexMatrix_col Theta_gate = ComplexApplyGateToTheta(Dl, Dr, d, Theta, G);
  auto split = ComplexsplitState(Theta_gate, chi, left_to_right);
  ComplexWriteBondFromSplit(psi, j, split);
  return split;
}

// Pipeline should be:
//   Theta      = ComplexBuildBondTheta(psi.A[j], psi.A[j + 1])
//   Theta_gate = ComplexApplyGateToTheta(Dl, Dr, d, Theta, G)
//   split      = ComplexsplitState(Theta_gate, chi, left_to_right)
//   ComplexWriteBondFromSplit(psi, j, split)
//   return split

// TODO(COMPLEX-3): Implement ComplexApplyBondSweep after Costruct SweepReport
struct SweepReport
{
  std::vector<double> discarded_weights;
  std::vector<double> discarded_fractions;
  std::vector<double> bond_entropies;
};

SweepReport ComplexApplyBondSweep(CMPS &psi, const ComplexMatrix_col &G, int chi, int start_bond, bool left_to_right)
{
  if (start_bond < 0 || start_bond > 1)
  {
    throw std::runtime_error("applyBondSweep expects start_bond to be 0 or 1");
  }

  SweepReport report;

  if (left_to_right)
  {
    for (int j = start_bond; j + 1 < psi.N; j += 2)
    {
      auto split = ComplexUpdateBond(psi, j, G, chi, left_to_right);
      report.discarded_weights.push_back(split.discarded_weight);
      report.discarded_fractions.push_back(split.discarded_fraction);
      report.bond_entropies.push_back(bondEntropyFromSchmidt(split.schmidt));
    }
  }
  else
  {
    int j0 = psi.N - 2;
    while (j0 >= 0 && (j0 % 2) != start_bond)
    {
      --j0;
    }
    for (int j = j0; j >= 0; j -= 2)
    {
      auto split = ComplexUpdateBond(psi, j, G, chi, left_to_right);
      report.discarded_weights.push_back(split.discarded_weight);
      report.discarded_fractions.push_back(split.discarded_fraction);
      report.bond_entropies.push_back(bondEntropyFromSchmidt(split.schmidt));
    }
  }

  return report;
}

SweepReport ComplexApplyBondSweepGates(
    CMPS &psi,
    const std::vector<ComplexMatrix_col> &gates,
    int chi,
    int start_bond,
    bool left_to_right)
{

  if (static_cast<int>(gates.size()) != psi.N - 1)
  {
    throw std::runtime_error("ComplexApplyBondSweepGates expects one gate per bond");
  }

  if (start_bond < 0 || start_bond > 1)
  {
    throw std::runtime_error("applyBondSweep expects start_bond to be 0 or 1");
  }

  SweepReport report;

  if (left_to_right)
  {
    for (int j = start_bond; j + 1 < psi.N; j += 2)
    {
      // important for implementing free Hamiltonian
      auto split = ComplexUpdateBond(psi, j, gates[j], chi, left_to_right);
      report.discarded_weights.push_back(split.discarded_weight);
      report.discarded_fractions.push_back(split.discarded_fraction);
      report.bond_entropies.push_back(bondEntropyFromSchmidt(split.schmidt));
    }
  }
  else
  {
    int j0 = psi.N - 2;
    while (j0 >= 0 && (j0 % 2) != start_bond)
    {
      --j0;
    }
    for (int j = j0; j >= 0; j -= 2)
    {
      auto split = ComplexUpdateBond(psi, j, gates[j], chi, left_to_right);
      report.discarded_weights.push_back(split.discarded_weight);
      report.discarded_fractions.push_back(split.discarded_fraction);
      report.bond_entropies.push_back(bondEntropyFromSchmidt(split.schmidt));
    }
  }

  return report;
}

ComplexMatrix_col mpsToFullState(const CMPS &psi)
{

  if (static_cast<int>(psi.A.size()) != psi.N)
  {
    throw std::runtime_error("mpsToFullState found inconsistent MPS length");
  }
  if (psi.N < 0)
  {
    throw std::runtime_error("mpsToFullState expects a positive N");
  }
  if (psi.d != 2)
  {
    throw std::runtime_error("mpsToFullState assumes local dimension =2");
  }
  // << :“shift bits to the left”. dim = 2^(psi.N) 1 << 1 = 2
  int dim = 1 << psi.N;
  ComplexMatrix_col state(dim, 1);

  for (int b = 0; b < dim; ++b)
  {
    // start with left boundary vector of size 1 (sweep from left to right -->)
    std::vector<Complex> v(1, Complex{1.0, 0.0});

    for (int j = 0; j < psi.N; ++j)
    {
      auto const &Aj = psi.A[j];
      if (static_cast<int>(v.size()) != Aj.Dl)
      {
        throw std::runtime_error("mpsToFullState found bond dimension mismatch while contracting");
      }
      //  b= 5 = 0101
      // s = (5 >> 2) & 1 -> 0001= 1 &1 = 1
      // &1 extracts the site bit moved to the far right s3
      // We keep states as
      // b   binary   s0 s1   ket
      // 0   00       0  0    |00>
      // 1   01       0  1    |01>
      // 2   10       1  0    |10>
      // 3   11       1  1    |11> ....
      int s = (b >> (psi.N - 1 - j)) & 1;

      // propagate from left bond to right bond space
      std::vector<Complex> next(Aj.Dr, Complex{0.0, 0.0});

      for (int l = 0; l < Aj.Dl; ++l)
      {
        for (int r = 0; r < Aj.Dr; ++r)
        {
          next[r] += v[l] * Aj(l, s, r);
        }
      }
      v = next;
    }

    if (static_cast<int>(v.size()) != 1)
    {
      throw std::runtime_error("mpsToFullState expected a scalar at the right boundary");
    }
    state(b, 0) = v[0];
  }

  return state;
}

// One-site expectation value
Complex oneSiteExpectationFullState(
    const ComplexMatrix_col &psi,
    int site,
    int N,
    const Matrix_col &op)
{
  if (psi.cols != 1 || psi.rows != (1 << N))
  {
    throw std::runtime_error("oneSiteExpectationFullState state shape mismatch");
  }
  if (site < 0 || site >= N)
  {
    throw std::runtime_error("oneSiteExpectationFullState site is out of range");
  }
  if (op.rows != 2 || op.cols != 2)
  {
    throw std::runtime_error("oneSiteExpectationFullState expects 2x2 operator (qubit)");
  }

  ComplexMatrix_col I_c = makeComplexIdentity(2);
  ComplexMatrix_col O_c(1, 1);
  O_c(0, 0) = Complex{1.0, 0.0};
  // O = factor_0 ⊗ factor_1 ⊗ ... ⊗ factor_{N-1}
  for (int i = 0; i < N; ++i)
  {
    ComplexMatrix_col factor = (i == site) ? toComplexMatrix(op) : makeComplexIdentity(2);
    O_c = complexKron(O_c, factor);
  }
  // O |\psi\rangle
  ComplexMatrix_col Opsi = cmatmul(O_c, psi);

  Complex numerator{0.0, 0.0};
  Complex denominator{0.0, 0.0};

  for (int i = 0; i < psi.rows; ++i)
  {
    numerator += std::conj(psi(i, 0)) * Opsi(i, 0);
    denominator += std::conj(psi(i, 0)) * psi(i, 0); // \langle \psi| \psi \rangle
  }

  if (std::abs(denominator) == 0.0)
  {
    throw std::runtime_error("oneSiteExpectationFullState yielded zero-norm state!");
  }

  return numerator / denominator;
}

// General reshape for MPS bond

Matrix_col reshapeStateVectorToMatrix(const Matrix_col &psi_vec, int left_dim, int right_dim)
{
  if (psi_vec.rows != left_dim * right_dim || psi_vec.cols != 1)
  {
    throw std::runtime_error("reshapeStateVectorToMatrix expects a (left_dim*right_dim) x 1 vector");
  }

  Matrix_col Theta(left_dim, right_dim);

  for (int s_l = 0; s_l < left_dim; ++s_l)
  {
    for (int s_r = 0; s_r < right_dim; ++s_r)
    {
      int flat = right_dim * s_l + s_r;
      Theta(s_l, s_r) = psi_vec(flat, 0);
    }
  }

  return Theta;
}

Matrix_col buildBondTheta(const SiteTensor &Aj, const SiteTensor &Aj1)
// required bond dimension >= rank(Theta)
// bond dimension = Schmidt rank = number of nonzero singular values
// So the bond dimension : size of the internal summed index in the factorization
{
  // left side (l,s1)
  // right side (s2,r)

  Matrix_col Theta(Aj.Dl * Aj.d, Aj1.d * Aj1.Dr);

  // Check
  if (Aj.Dr != Aj1.Dl)
  {
    throw std::runtime_error("buildBondTheta bond mismatch");
  }

  int S1 = Aj.d;
  int S2 = Aj1.d;
  // sum index for Theta(l,s1,s2,r) = sum_m Aj(l,s1,m) * Aj1(m,s2,r)
  int Dm = Aj.Dr;

  for (int l = 0; l < Aj.Dl; ++l)
  {
    for (int r = 0; r < Aj1.Dr; ++r)
    {
      for (int s1 = 0; s1 < S1; ++s1)
      {
        for (int s2 = 0; s2 < S2; ++s2)
        {
          double val = 0.0;
          for (int m = 0; m < Dm; ++m)
          {
            val += Aj(l, s1, m) * Aj1(m, s2, r);
          }
          int row = l * S1 + s1;
          int col = s2 * Aj1.Dr + r;
          Theta(row, col) = val;
        }
      }
    }
  }
  return Theta;
}

void writeBondFromSplit(MPS &psi, int j, const SplitStateResult &split)
{
  int Dl = psi.A[j].Dl;
  int d = psi.A[j].d;
  int Dr = psi.A[j + 1].Dr;
  int chi = split.schmidt.size(); // BOND DIMENSION!

  // updated tensors
  SiteTensor Aj_new(Dl, d, chi);
  SiteTensor Aj1_new(chi, d, Dr);

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

Matrix_col applyGateToTheta(int Dl, int Dr, int d, const Matrix_col &Theta, const Matrix_col &G)
{
  if (Theta.rows != Dl * d || Theta.cols != d * Dr)
  {
    throw std::runtime_error("applyGateToTheta shape mismatch");
  }
  if (G.rows != d * d || G.cols != d * d)
  {
    throw std::runtime_error("applyGateToTheta expects a (d*d)x(d*d) gate");
  }

  Matrix_col Theta_new(Dl * d, d * Dr);

  for (int l = 0; l < Dl; ++l)
  {
    for (int r = 0; r < Dr; ++r)
    {
      Matrix_col theta_lr(d * d, 1);

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

      Matrix_col theta_lr_new = applyGate(G, theta_lr);

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

SplitStateResult updateBond(MPS &psi, int j, const Matrix_col &G, int chi, bool left_to_right)
{
  Matrix_col Theta = buildBondTheta(psi.A[j], psi.A[j + 1]);
  int Dl = psi.A[j].Dl;
  int d = psi.A[j].d;
  int Dr = psi.A[j + 1].Dr;

  // SVD with chi bond dimension
  Matrix_col Theta_gate = applyGateToTheta(Dl, Dr, d, Theta, G);
  auto split = splitState(Theta_gate, chi, left_to_right);
  writeBondFromSplit(psi, j, split);
  return split;
}

std::vector<int> bondDimensions(const MPS &psi)
{
  if (static_cast<int>(psi.A.size()) != psi.N)
  {
    throw std::runtime_error("bondDimensions found inconsistent MPS length");
  }

  std::vector<int> dims;
  if (psi.N <= 1)
  {
    return dims;
  }

  dims.reserve(psi.N - 1);
  for (int j = 0; j + 1 < psi.N; ++j)
  {
    if (psi.A[j].Dr != psi.A[j + 1].Dl)
    {
      throw std::runtime_error("bondDimensions found mismatched neighboring bond dimensions");
    }
    dims.push_back(psi.A[j].Dr);
  }
  return dims;
}

std::vector<int> ComplexBondDimensions(const CMPS &psi)
{
  if (static_cast<int>(psi.A.size()) != psi.N)
  {
    throw std::runtime_error("bondDimensions found inconsistent MPS length");
  }

  std::vector<int> dims;
  if (psi.N <= 1)
  {
    return dims;
  }

  dims.reserve(psi.N - 1);
  for (int j = 0; j + 1 < psi.N; ++j)
  {
    if (psi.A[j].Dr != psi.A[j + 1].Dl)
    {
      throw std::runtime_error("bondDimensions found mismatched neighboring bond dimensions");
    }
    dims.push_back(psi.A[j].Dr);
  }
  return dims;
}

void printBondDimensions(const MPS &psi, const std::string &name)
{
  printIntVector(bondDimensions(psi), name);
}

void printComplexBondDimensions(const CMPS &psi, const std::string &name)
{
  printIntVector(ComplexBondDimensions(psi), name);
}

// UpdateBond across the sites of the MPS psi
SweepReport applyBondSweep(MPS &psi, const Matrix_col &G, int chi, int start_bond, bool left_to_right)
{
  if (start_bond < 0 || start_bond > 1)
  {
    throw std::runtime_error("applyBondSweep expects start_bond to be 0 or 1");
  }

  SweepReport report;

  if (left_to_right)
  {
    for (int j = start_bond; j + 1 < psi.N; j += 2)
    {
      auto split = updateBond(psi, j, G, chi, left_to_right);
      report.discarded_weights.push_back(split.discarded_weight);
      report.discarded_fractions.push_back(split.discarded_fraction);
      report.bond_entropies.push_back(bondEntropyFromSchmidt(split.schmidt));
    }
  }
  else
  {
    int j0 = psi.N - 2;
    while (j0 >= 0 && (j0 % 2) != start_bond)
    {
      --j0;
    }
    for (int j = j0; j >= 0; j -= 2)
    {
      auto split = updateBond(psi, j, G, chi, left_to_right);
      report.discarded_weights.push_back(split.discarded_weight);
      report.discarded_fractions.push_back(split.discarded_fraction);
      report.bond_entropies.push_back(bondEntropyFromSchmidt(split.schmidt));
    }
  }

  return report;
}

void runLegacySanityChecks(const Matrix_col &h_bond, const Matrix_col &G, const ComplexMatrix_col &Urt)
{
  Matrix_col psi(4, 1);
  psi(0, 0) = 1.0;
  psi(1, 0) = 0.0;
  psi(2, 0) = 0.0;
  psi(3, 0) = 0.0;

  Matrix_col hpsi = matmul(h_bond, psi);
  printMatrix_col(psi, "psi");
  printMatrix_col(hpsi, "hpsi");

  auto eig = eighSymmetric(h_bond);
  printVector(eig.evals, "evals");
  printMatrix_col(eig.evecs, "evecs");

  Matrix_col psi_dt = applyGate(G, psi);
  printMatrix_col(psi_dt, "psi_dt");

  Matrix_col Theta = reshapeStateVectorToMatrix(psi_dt, 2, 2);
  printMatrix_col(Theta, "Theta");

  auto svdTheta = svd(Theta);
  printMatrix_col(svdTheta.U, "U_theta");
  printVector(svdTheta.S, "S_theta");
  printMatrix_col(svdTheta.VT, "VT_theta");

  Matrix_col Sigma_theta = makeDiagonal(svdTheta.S);
  Matrix_col Theta_rec = matmul(matmul(svdTheta.U, Sigma_theta), svdTheta.VT);
  printMatrix_col(Theta_rec, "Theta_rec");
  std::cout << "Theta reconstruction error:"
            << maxAbsDiff(Theta, Theta_rec) << "\n";

  int chi_theta = 1;
  auto split = splitState(Theta, chi_theta, true);

  Matrix_col Theta_split = matmul(split.left, split.right);
  printMatrix_col(split.left, "Left_trunc");
  printMatrix_col(split.right, "Right_trunc");
  printVector(split.schmidt, "S_trunc");
  printMatrix_col(Theta_split, "Theta_split");
  std::cout << "max truncation error: "
            << maxAbsDiff(Theta_split, Theta) << "\n";
  std::cout << "discarded weight: "
            << split.discarded_weight << "\n";
  std::cout << "discarded weight fraction: "
            << split.discarded_fraction << "\n";
  std::cout << "full bond entropy: "
            << bondEntropyFromSchmidt(svdTheta.S) << "\n";
  std::cout << "truncated bond entropy: "
            << bondEntropyFromSchmidt(split.schmidt) << "\n";

  SiteTensor P1(1, 2, 1);
  SiteTensor P2(1, 2, 1);

  P1(0, 0, 0) = 1.0;
  P2(0, 1, 0) = 1.0;

  Matrix_col Theta_prod = buildBondTheta(P1, P2);
  printMatrix_col(Theta_prod, "Theta_product_01");

  SiteTensor B1(1, 2, 2);
  SiteTensor B2(2, 2, 1);

  double invsqrt2 = 1.0 / std::sqrt(2.0);

  B1(0, 0, 0) = invsqrt2;
  B1(0, 1, 1) = invsqrt2;

  B2(0, 0, 0) = 1.0;
  B2(1, 1, 0) = 1.0;

  Matrix_col Theta_bell = buildBondTheta(B1, B2);
  printMatrix_col(Theta_bell, "Theta_bell");

  MPS bell_mps;
  bell_mps.N = 2;
  bell_mps.d = 2;
  bell_mps.A.push_back(B1);
  bell_mps.A.push_back(B2);

  Matrix_col Theta_before = buildBondTheta(bell_mps.A[0], bell_mps.A[1]);
  printMatrix_col(Theta_before, "Theta_before_writeback");

  auto split_bell = splitState(Theta_before, 2, true);
  writeBondFromSplit(bell_mps, 0, split_bell);

  Matrix_col Theta_after = buildBondTheta(bell_mps.A[0], bell_mps.A[1]);
  printMatrix_col(Theta_after, "Theta_after_writeback");

  std::cout << "Bell writeback error: "
            << maxAbsDiff(Theta_before, Theta_after) << "\n";

  Matrix_col I4(4, 4);
  for (int i = 0; i < 4; ++i)
  {
    I4(i, i) = 1.0;
  }

  MPS bell_mps_id;
  bell_mps_id.N = 2;
  bell_mps_id.d = 2;
  bell_mps_id.A.push_back(B1);
  bell_mps_id.A.push_back(B2);

  Matrix_col Theta_id_before = buildBondTheta(bell_mps_id.A[0], bell_mps_id.A[1]);
  updateBond(bell_mps_id, 0, I4, 2, false);
  Matrix_col Theta_id_after = buildBondTheta(bell_mps_id.A[0], bell_mps_id.A[1]);

  printMatrix_col(Theta_id_before, "Theta_id_before");
  printMatrix_col(Theta_id_after, "Theta_id_after");
  std::cout << "Identity-gate update error: "
            << maxAbsDiff(Theta_id_before, Theta_id_after) << "\n";

  MPS bell_mps_h;
  bell_mps_h.N = 2;
  bell_mps_h.d = 2;
  bell_mps_h.A.push_back(B1);
  bell_mps_h.A.push_back(B2);

  MPS bell_mps_h_trunc;
  bell_mps_h_trunc.N = 2;
  bell_mps_h_trunc.d = 2;
  bell_mps_h_trunc.A.push_back(B1);
  bell_mps_h_trunc.A.push_back(B2);

  auto split_h = updateBond(bell_mps_h, 0, G, 2, false);
  Matrix_col Theta_h_after = buildBondTheta(bell_mps_h.A[0], bell_mps_h.A[1]);

  auto split_h_trunc = updateBond(bell_mps_h_trunc, 0, G, 1, false);
  Matrix_col Theta_h_trunc_after = buildBondTheta(bell_mps_h_trunc.A[0], bell_mps_h_trunc.A[1]);

  printMatrix_col(Theta_h_after, "Theta_h_after");
  printMatrix_col(Theta_h_trunc_after, "Theta_h__trunc_after");
  std::cout << "bell state gate difference: "
            << maxAbsDiff(Theta_h_trunc_after, Theta_h_after) << "\n";

  printVector(split_h.schmidt, "bell gate Schmidt");
  std::cout << "bell gate discarded weight: "
            << split_h.discarded_weight << "\n";
  printVector(split_h_trunc.schmidt, "bell gate Schmidt trunc");
  std::cout << "bell gate discarded weight trunc: "
            << split_h_trunc.discarded_weight << "\n";

  // Complex Linear Algebra checks
  Matrix_col D_real(2, 2);
  D_real(1, 0) = -1;
  D_real(0, 1) = 1;

  Matrix_col E_real(2, 2);
  E_real(1, 0) = 1;
  E_real(0, 1) = 1;

  ComplexMatrix_col D_c = toComplexMatrix(D_real);
  ComplexMatrix_col E_c = toComplexMatrix(E_real);
  ComplexMatrix_col F = cmatmul(D_c, E_c);
  printComplexMatrix_col(F, "F matrix");

  ComplexMatrix_col check = cmatmul(Urt, conjTranspose(Urt));
  printComplexMatrix_col(check, "Urt Urt^dag");
  std::cout << "unitarity error: "
            << maxAbsDiffComplex(check, makeComplexIdentity(Urt.rows)) << "\n";
}

void runNeelSweepSanityCheck(const IsingParams &p, int chi)
{
  // imaginary time TEBD sweep with Neel state

  MPS neel_chain = makeNeelState(8);

  printProductConfiguration(neel_chain, "neel chain before sweep");
  printBondDimensions(neel_chain, "Bond dims before sweep");

  // even sites sweep

  // 1. build bond Hamiltonian and gate
  Matrix_col h_bond = makeBondHamiltonian(p.J, p.h);
  Matrix_col G = makeImagTimeGate(h_bond, p.dt);

  auto even_report = applyBondSweep(neel_chain, G, chi, 0, true);
  printBondDimensions(neel_chain, "Bond dims after even sweep");
  printVector(even_report.discarded_weights, "Even discarded weights");
  printVector(even_report.discarded_fractions, "Even discarded fractions");
  printVector(even_report.bond_entropies, "Even bond entropies");

  // odd sites sweep

  auto odd_report = applyBondSweep(neel_chain, G, chi, 1, false);
  printBondDimensions(neel_chain, "Bond dims after odd sweep");
  printVector(odd_report.discarded_weights, "Odd discarded weights");
  printVector(odd_report.discarded_fractions, "Odd discarded fractions");
  printVector(odd_report.bond_entropies, "Odd bond entropies");
}

void runComplexNeelSweepSanityCheck(const IsingParams &p, int chi)
{
  // imaginary time TEBD sweep with Neel state

  CMPS neel_chain = makeComplexNeelState(8);

  // printProductConfiguration(neel_chain, "neel chain before sweep");
  printComplexBondDimensions(neel_chain, "Bond dims before sweep");

  // even sites sweep

  // 1. build bond Hamiltonian and gate
  bool GorGates = false;
  if (GorGates)
  {
    Matrix_col h_bond = makeBondHamiltonian(p.J, p.h);
    ComplexMatrix_col G = makeRealTimeGate(h_bond, p.dt);

    auto even_report = ComplexApplyBondSweep(neel_chain, G, chi, 0, true);
    printComplexBondDimensions(neel_chain, "Bond dims after even sweep");
    printVector(even_report.discarded_weights, "Even discarded weights");
    printVector(even_report.discarded_fractions, "Even discarded fractions");
    printVector(even_report.bond_entropies, "Even bond entropies");

    // odd sites sweep

    auto odd_report = ComplexApplyBondSweep(neel_chain, G, chi, 1, false);
    printComplexBondDimensions(neel_chain, "Bond dims after odd sweep");
    printVector(odd_report.discarded_weights, "Odd discarded weights");
    printVector(odd_report.discarded_fractions, "Odd discarded fractions");
    printVector(odd_report.bond_entropies, "Odd bond entropies");
  }
  else
  {
    auto gates = makeRealTimeGates(neel_chain.N, p.J, p.h, p.dt);
    auto even_c = ComplexApplyBondSweepGates(neel_chain, gates, chi, 0, true);
    printComplexBondDimensions(neel_chain, "Bond dims after even sweep");
    printVector(even_c.discarded_weights, "Even discarded weights");
    printVector(even_c.discarded_fractions, "Even discarded fractions");
    printVector(even_c.bond_entropies, "Even bond entropies");

    auto odd_c = ComplexApplyBondSweepGates(neel_chain, gates, chi, 1, false);
    printComplexBondDimensions(neel_chain, "Bond dims after odd sweep");
    printVector(odd_c.discarded_weights, "Odd discarded weights");
    printVector(odd_c.discarded_fractions, "Odd discarded fractions");
    printVector(odd_c.bond_entropies, "Odd bond entropies");
  }
}

double sumEntries(const std::vector<double> &v)
{
  double sum = 0.0;
  for (double v_ent : v)
  {
    sum += v_ent;
  }
  return sum;
}

double maxEntry(const std::vector<double> &v)
{
  double ans = 0.0;

  for (double j : v)
  {
    ans = std::max(ans, j);
  }
  return ans;
}

struct StepSummary
{
  int step = 0;
  double tau = 0.0;
  double discarded_even = 0.0;
  double discarded_odd = 0.0;
  double max_entropy_even = 0.0;
  double max_entropy_odd = 0.0;
  double max_discarded_frac_even = 0.0;
  double max_discarded_frac_odd = 0.0;

  std::vector<int> bond_dims;
};

struct ImagTimeRun
{
  MPS psi;
  std::vector<StepSummary> history;
};

ImagTimeRun runImagTimeIsing(const IsingParams &p, int chi)
{
  // 1. build bond Hamiltonian and gate
  Matrix_col h_bond = makeBondHamiltonian(p.J, p.h);
  Matrix_col G = makeImagTimeGate(h_bond, p.dt);

  // 2. initialize the state
  MPS psi = makeNeelState(p.N);

  // 3. prepare output container
  ImagTimeRun run;

  // 4. time-step loop
  for (int step = 0; step < p.nsteps; ++step)
  {
    // even sweep
    auto even_report = applyBondSweep(psi, G, chi, 0, true);

    // odd sweep
    auto odd_report = applyBondSweep(psi, G, chi, 1, false);

    StepSummary rec;
    rec.step = step + 1;
    rec.tau = (step + 1) * p.dt;
    rec.max_discarded_frac_even = maxEntry(even_report.discarded_fractions);
    rec.max_discarded_frac_odd = maxEntry(odd_report.discarded_fractions);

    // TODO: fill rec.discarded_even
    // hint: sum over even_report.discarded_weights

    rec.discarded_even = sumEntries(even_report.discarded_weights);

    // TODO: fill rec.discarded_odd
    // hint: sum over odd_report.discarded_weights
    rec.discarded_odd = sumEntries(odd_report.discarded_weights);

    // TODO: fill rec.max_entropy_even
    // hint: max over even_report.bond_entropies
    rec.max_entropy_even = maxEntry(even_report.bond_entropies);

    // TODO: fill rec.max_entropy_odd
    // hint: max over odd_report.bond_entropies
    rec.max_entropy_odd = maxEntry(odd_report.bond_entropies);
    // TODO: fill rec.bond_dims
    // hint: use bondDimensions(psi)
    rec.bond_dims = bondDimensions(psi);
    // TODO: push rec into run.history
    run.history.push_back(rec);
  }

  // 5. store final state
  // TODO: set run.psi = psi
  run.psi = psi;

  return run;
}

// CSV benchmark driver
void runJ0BenchmarkCsv(
    int N,
    double h,
    const std::vector<double> &dts,
    double tmax,
    int chi,
    const std::string &csv_path)
{
  // open std::ofstream
  std::ofstream out(csv_path);
  if (!out)
  {
    throw std::runtime_error("failed to open CSV file");
  }

  out << std::setprecision(17);
  out << "dt,step,t,site,z_tebd,z_exact,z_analytic,abs_err_exact,abs_err_analytic\n";

  Matrix_col Z(2, 2);
  Z(0, 0) = 1.0;
  Z(1, 1) = -1.0;

  Matrix_col H = makeFullHamiltonianJ0(N, h);

  for (double dt : dts)
  {
    int nsteps = static_cast<int>(std::round(tmax / dt));

    CMPS psi_tebd = makeComplexNeelState(N);
    ComplexMatrix_col psi_exact = mpsToFullState(psi_tebd);

    auto gates = makeRealTimeGates(N, 0.0, h, dt);

    for (int step = 0; step <= nsteps; ++step)
    {
      double t = step * dt;

      ComplexMatrix_col psi_tebd_full = mpsToFullState(psi_tebd);

      for (int site = 0; site < N; ++site)
      {
        Complex z_tebd_c = oneSiteExpectationFullState(psi_tebd_full, site, N, Z);
        Complex z_exact_c = oneSiteExpectationFullState(psi_exact, site, N, Z);
        double z_tebd = z_tebd_c.real();
        double z_exact = z_exact_c.real();

        // analytic result
        double sigma = (site % 2 == 0) ? 1.0 : -1.0;
        double z_analytic = sigma * std::cos(2.0 * h * t);

        // CSV read
        double abs_err_exact = std::abs(z_tebd - z_exact);
        double abs_err_analytic = std::abs(z_tebd - z_analytic);

        out << dt << ","
            << step << ","
            << t << ","
            << site << ","
            << z_tebd << ","
            << z_exact << ","
            << z_analytic << ","
            << abs_err_exact << ","
            << abs_err_analytic << "\n";
      }
      if (step == nsteps)
      {
        break;
      }
      //  exact update
      psi_exact = evolveExactOneStep(psi_exact, H, dt);

      // TEBD update
      // even sweep
      ComplexApplyBondSweepGates(psi_tebd, gates, chi, 0, true);
      // odd sweep
      ComplexApplyBondSweepGates(psi_tebd, gates, chi, 1, false);
    }
  }
}

double twoSiteExpectation(const MPS &psi, int j, const Matrix_col &op2)
{
  (void)psi;
  (void)j;
  (void)op2;
  throw std::runtime_error("TODO: implement nearest-neighbor expectation value <O_{j,j+1}>");
}

int main()
{
  auto p = IsingParams{};
  std::cout << "Non-ITensor manual track placeholder.\n";
  std::cout << "Goal: implement the tensor-network machinery directly on BLAS/LAPACK.\n";
  std::cout << "Current defaults: N=" << p.N
            << ", J=" << p.J
            << ", h=" << p.h
            << ", nsteps=" << p.nsteps
            << ", dt=" << p.dt
            << "\n";

  // bond hamiltonian
  // Matrix addition
  double J = p.J;
  double hx = p.h;
  double dt = p.dt;

  Matrix_col h_bond = makeBondHamiltonian(J, hx);
  // print
  printMatrix_col(h_bond, "h_bond");

  Matrix_col G = makeImagTimeGate(h_bond, dt);
  printMatrix_col(G, "G");

  ComplexMatrix_col Urt = makeRealTimeGate(h_bond, dt);

  bool SanityChecks = false;
  if (SanityChecks)
  {
    runLegacySanityChecks(h_bond, G, Urt);
  }
  bool SanityNeelSweep = false;
  if (SanityNeelSweep)
  {
    runNeelSweepSanityCheck(p, 2);
  }

  // Implement runImagTimeIsing
  bool SanityNeelImagTime = false;
  if (SanityNeelImagTime)
  {
    auto run = runImagTimeIsing(p, 2);
    for (auto const &rec : run.history)
    {
      std::cout << "step " << rec.step
                << " even_dw=" << rec.discarded_even
                << " odd_dw=" << rec.discarded_odd
                << " even_Smax=" << rec.max_entropy_even
                << " odd_Smax=" << rec.max_entropy_odd
                << ", tau=" << rec.tau
                << " even_df=" << rec.max_discarded_frac_even
                << " odd_df=" << rec.max_discarded_frac_odd
                << "\n";
      printIntVector(rec.bond_dims, "bond dims");
    }
  }

  // Complex SVD sanity check

  bool SanityComplexSVD = false;
  if (SanityComplexSVD)
  {
    ComplexMatrix_col A(2, 2);
    A(0, 0) = Complex{1.0, 0.0};
    A(0, 1) = Complex{0.0, 1.0};
    A(1, 0) = Complex{1.0, -1.0};
    A(1, 1) = Complex{-1.0, 0.5};

    auto svdA = complexSvd(A);
    ComplexMatrix_col S = makeComplexDiagonalFromReal(svdA.S);
    ComplexMatrix_col Arec = cmatmul(cmatmul(svdA.U, S), svdA.VH);

    printComplexMatrix_col(A, "A_complex");
    printComplexMatrix_col(Arec, "Arec_complex");
    std::cout << "complex SVD reconstruction error: "
              << maxAbsDiffComplex(A, Arec) << "\n";

    ComplexMatrix_col Ucheck = cmatmul(conjTranspose(svdA.U), svdA.U);
    ComplexMatrix_col Vcheck = cmatmul(svdA.VH, conjTranspose(svdA.VH));

    std::cout << "U unitary error: "
              << maxAbsDiffComplex(Ucheck, makeComplexIdentity(svdA.U.cols)) << "\n";
    std::cout << "VH unitary error: "
              << maxAbsDiffComplex(Vcheck, makeComplexIdentity(svdA.VH.rows)) << "\n";
  }

  // Complex tEBD check
  bool SanityComplexBond = false;
  if (SanityComplexBond)
  {
    // 1. Build  a2-site complex bell state
    CMPS bell_c;
    bell_c.N = 2;
    bell_c.d = 2;

    ComplexSiteTensor B1c(1, 2, 2);
    ComplexSiteTensor B2c(2, 2, 1);

    double invsqrt2 = 1.0 / std::sqrt(2.0);

    B1c(0, 0, 0) = Complex{invsqrt2, 0.0};
    B1c(0, 1, 1) = Complex{invsqrt2, 0.0};

    B2c(0, 0, 0) = Complex{1.0, 0.0};
    B2c(1, 1, 0) = Complex{1.0, 0.0};

    bell_c.A.push_back(B1c);
    bell_c.A.push_back(B2c);

    // 2. Combine into a two-site Theta
    ComplexMatrix_col Theta_before = ComplexBuildBondTheta(bell_c.A[0], bell_c.A[1]);

    // 3. Apply one complex bond update
    ComplexSplitStateResult split_c = ComplexUpdateBond(bell_c, 0, Urt, 2, false);

    // 4. Rebuild Theta after update
    ComplexMatrix_col Theta_after = ComplexBuildBondTheta(bell_c.A[0], bell_c.A[1]);

    // Check with no truncation
    int Dl = bell_c.A[0].Dl;
    int d = bell_c.A[0].d;
    int Dr = bell_c.A[1].Dr;

    ComplexMatrix_col Theta_direct = ComplexApplyGateToTheta(Dl, Dr, d, Theta_before, Urt);

    // 5. Print diagnostics

    printComplexMatrix_col(Theta_before, "Theta_before");
    printComplexMatrix_col(Theta_after, "Theta_after");
    std::cout << " Truncation error: "
              << maxAbsDiffComplex(Theta_after, Theta_direct) << "\n";

    printVector(split_c.schmidt, "complex bond Schmidt");
    std::cout << "complex discarded weight:"
              << split_c.discarded_weight << "\n";
    std::cout << "complex discarded fraction:"
              << split_c.discarded_fraction << "\n";
  }

  // Complex Neel Sweep
  bool ComplexSanityNeelSweep = false;
  if (ComplexSanityNeelSweep)
  {
    runComplexNeelSweepSanityCheck(p, 2);
  }

  // mpsToFullState

  bool MPSToFull = false;
  if (MPSToFull)
  {
    //     ( |00> + |11> ) / sqrt(2)
    // =
    // [1/sqrt(2), 0, 0, 1/sqrt(2)]^T
    // b = 0 -> |00>
    // b = 1 -> |01>
    // b = 2 -> |10>
    // b = 3 -> |11>

    CMPS bell_c;
    bell_c.N = 2;
    bell_c.d = 2;

    ComplexSiteTensor B1c(1, 2, 2);
    ComplexSiteTensor B2c(2, 2, 1);

    double invsqrt2 = 1.0 / std::sqrt(2.0);

    B1c(0, 0, 0) = Complex{invsqrt2, 0.0};
    B1c(0, 1, 1) = Complex{invsqrt2, 0.0};

    B2c(0, 0, 0) = Complex{1.0, 0.0};
    B2c(1, 1, 0) = Complex{1.0, 0.0};

    bell_c.A.push_back(B1c);
    bell_c.A.push_back(B2c);

    ComplexMatrix_col bell_c_state = mpsToFullState(bell_c);
    printComplexMatrix_col(bell_c_state, "bell mps to state");
  }

  bool runExactCheck = true;
  if (runExactCheck)
  {
    int N_bench = 6;
    double tmax = 2.0;
    std::vector<double> dts{0.2, 0.1, 0.05, 0.02, 0.01};
    runJ0BenchmarkCsv(
        N_bench,
        p.h,
        dts,
        tmax,
        16,
        "/Users/juan/Desktop/Git/FFD_GGE/tutorials/ising_manual/j0_benchmark.csv");
  }

  return 0;
}
