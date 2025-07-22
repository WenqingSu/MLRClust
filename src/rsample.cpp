#include <RcppEigen.h>

using Rcpp::NumericVector;
using Rcpp::IntegerVector;

using SpMat = Eigen::SparseMatrix<double>;
using MapSpMat = Eigen::Map<SpMat>;

// [[Rcpp::export]]
Rcpp::S4 rsample(Rcpp::S4 A, double p)
{
  // A simpler but less efficient implementation
  /* SpMat A2 = Rcpp::as<SpMat>(A);
   const int nnz = A2.nonZeros();
   double* xptr = A2.valuePtr();
   for(int i = 0; i < nnz; i++)
   {
   xptr[i] = double(R::unif_rand() < p);
   }

   A2.prune(0.0);
   A2.makeCompressed();
   return Rcpp::wrap(A2); */

  MapSpMat A2 = Rcpp::as<MapSpMat>(A);
  const int n = A2.cols();
  const int nnz = A2.nonZeros();
  const int* inner = A2.innerIndexPtr();
  const int* outer = A2.outerIndexPtr();

  IntegerVector slotp(n + 1);  // will be initialized to 0
  int* new_outer = slotp.begin();
  std::vector<int> new_inner;
  new_inner.reserve(nnz * p * 1.2);

  for(int j = 0; j < n; j++)
  {
    const int* Ai_start = inner + outer[j];
    const int* Ai_end = inner + outer[j + 1];
    for(; Ai_start < Ai_end; Ai_start++)
    {
      if(R::unif_rand() <= p)
      {
        new_inner.emplace_back(*Ai_start);
        (new_outer[j + 1])++;
      }
    }
    new_outer[j + 1] += new_outer[j];
  }

  Rcpp::S4 res("dgCMatrix");
  NumericVector slotx(new_inner.size(), 1.0);

  res.slot("i") = Rcpp::wrap(new_inner);
  res.slot("p") = slotp;
  res.slot("Dim") = IntegerVector::create(n, n);
  res.slot("x") = slotx;

  return res;
}

// [[Rcpp::export]]
Rcpp::S4 rsample_sym(Rcpp::S4 A, double p, bool use_lower = true)
{
  SpMat A2;
  if(use_lower)
    A2 = Rcpp::as<SpMat>(A).triangularView<Eigen::StrictlyLower>();
  else
    A2 = Rcpp::as<SpMat>(A).triangularView<Eigen::StrictlyUpper>().transpose();
  const int nnz = A2.nonZeros();
  double* xptr = A2.valuePtr();
  for(int i = 0; i < nnz; i++)
  {
    xptr[i] = double(R::unif_rand() < p);
  }

  A2.prune(0.0);
  A2.makeCompressed();

  return Rcpp::wrap(A2);
}
