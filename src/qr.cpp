#include <Rcpp.h>
#include <R_ext/Lapack.h>

using Rcpp::NumericMatrix;
using Rcpp::List;


// X = QR, X[m x n], m >= n
// Overwrite X and an additional vector tau [n] with an internal representation of Q and R
void qr_inplace_impl(double* x, double* tau, const int m, const int n)
{
  double work_query;
  int lwork = -1;
  int info;
  // Query workspace size
  F77_CALL(dgeqrf)(&m, &n, x, &m, tau, &work_query, &lwork, &info);

  // QR decomposition
  lwork = (int)(work_query);
  double* work = new double[lwork];
  F77_CALL(dgeqrf)(&m, &n, x, &m, tau, work, &lwork, &info);
  delete[] work;
}

void qr_Q_inplace_impl(double* qr, const double* tau, const int m, const int n)
{
  double work_query;
  int lwork = -1;
  int info;
  // Query workspace size
  F77_CALL(dorgqr)(&m, &n, &n, qr, &m, tau, &work_query, &lwork, &info);

  // Extract Q
  lwork = (int)(work_query);
  double* work = new double[lwork];
  F77_CALL(dorgqr)(&m, &n, &n, qr, &m, tau, work, &lwork, &info);
  delete[] work;
}

// In-place QR decomposition
// [[Rcpp::export]]
void qr_Q_inplace(NumericMatrix x)
{
  const int m = x.nrow();
  const int n = x.ncol();
  if(m < n)
    Rcpp::stop("nrow(x) must be greater than or equal to ncol(x)");

  double* tau = new double[n];
  qr_inplace_impl(x.begin(), tau, m, n);
  qr_Q_inplace_impl(x.begin(), tau, m, n);

  delete [] tau;
}
