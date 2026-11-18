#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Number of threads actually used: 1 when the package was built without OpenMP.
static int n_threads(int ncores) {
#ifdef _OPENMP
  return ncores < 1 ? 1 : ncores;
#else
  (void) ncores;
  return 1;
#endif
}

// Shared worker: y is a q x n matrix (one observation per column, as y <- t(x) in R).
// Returns f[i] = sum_{j != i} exp( s * (<y_i, y_j> - 1) ) * com.
// Subtracting 1 keeps the exponent <= 0, so exp() cannot overflow for large s;
// com must therefore use the exponentially scaled Bessel function.
// The dot products are symmetric, so each pair (i, j) is computed once
// and added to both f[i] and f[j]: n(n-1)/2 exp() calls instead of n(n-1).
// With several threads each thread accumulates into its own vector of length n
// (no data races), and the vectors are summed at the end. Rows i near the top
// of the triangle carry more work, hence the dynamic schedule.
static NumericVector kernel_loop(const NumericMatrix& y, const double s, const double com,
                                 const int ncores) {
  const int q = y.nrow(), n = y.ncol();
  const double* Y = y.begin();          // column-major: column i starts at Y + i * q
  const int nt = n_threads(ncores);
  std::vector<double> acc( (std::size_t)n * nt, 0.0 );

#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(dynamic, 16)
#endif
  for (int i = 0; i < n - 1; ++i) {
#ifdef _OPENMP
    double* a = acc.data() + (std::size_t)omp_get_thread_num() * n;
#else
    double* a = acc.data();
#endif
    const double* yi = Y + (std::size_t)i * q;
    double ai = 0.0;
    for (int j = i + 1; j < n; ++j) {
      const double* yj = Y + (std::size_t)j * q;
      double d = 0.0;
      for (int k = 0; k < q; ++k) d += yi[k] * yj[k];
      const double e = std::exp(s * (d - 1.0));
      ai += e;
      a[j] += e;
    }
    a[i] += ai;
  }

  NumericVector f(n);
  for (int i = 0; i < n; ++i) {
    double t = 0.0;
    for (int r = 0; r < nt; ++r) t += acc[(std::size_t)r * n + i];
    f[i] = t * com;
  }
  return f;
}

// vm.kernel(): replaces
//   for (i in 1:n) {
//     disa <- Rfast::eachcol.apply(y, y[, i])[-i]
//     f[i] <- sum( exp(disa)^h ) / com
//   }
// y is the 2 x n matrix rbind(cos(u), sin(u)); exp(disa)^h = exp(h * disa).
// com = (n - 1) * 2 * pi * besselI(h, 0, expon.scaled = TRUE).
// [[Rcpp::export]]
NumericVector vm_kernel_cpp(NumericMatrix y, double h, double com, int ncores = 1) {
  return kernel_loop(y, h, 1.0 / com, ncores);
}

// vmf.kernel(): replaces
//   for ( i in 1:n ) {
//     disa <- Rfast::eachcol.apply(y, y[, i])[-i] / h^2
//     f[i] <- sum( exp(disa) ) * com
//   }
// y is the p x n matrix t(x). com uses besselI(1/h^2, p/2 - 1, expon.scaled = TRUE).
// [[Rcpp::export]]
NumericVector vmf_kernel_cpp(NumericMatrix y, double h, double com, int ncores = 1) {
  return kernel_loop(y, 1.0 / (h * h), com, ncores);
}

// Kernel discriminant analysis scores for vm.kernel.da() and vmf.kernel.da().
// Replaces, for every class j,
//   a <- kappa[j] * tcrossprod( xnew, x[ina == j, , drop = FALSE] )
//   mx <- Rfast::rowMaxs(a, value = TRUE)
//   score[, j] <- cons[j] + mx + log( Rfast::rowsums( exp(a - mx) ) )
// without ever forming the m x n_j matrix a: each new observation is processed
// on its own with a streaming (one-pass) log-sum-exp, so the extra memory is
// O(k) per row, whatever the number of new observations m.
// The rows of ynew are independent, so they are split across ncores threads.
// They are processed in blocks so that user interrupts can be checked (from
// the main thread only) between blocks.
//
// ynew : p x m matrix, one new observation per column (t(xnew) in R)
// y    : p x n matrix, one training observation per column (t(x) in R)
// ina  : integer labels 1..k, length n
// kappa: concentration per class (h for vm, 1/h^2 for vmF)
// cons : constant per class, log(prior) - log(n_j) + log C(kappa_j)
// Returns list(class, score); class is the first column attaining the row
// maximum, as in Rfast::rowMaxs().
// [[Rcpp::export]]
List kernel_da_cpp(NumericMatrix ynew, NumericMatrix y, IntegerVector ina,
                   NumericVector kappa, NumericVector cons, int ncores = 1) {
  const int p = y.nrow(), n = y.ncol(), m = ynew.ncol(), k = kappa.size();
  if ( ynew.nrow() != p ) stop("ynew and y must have the same dimension");
  if ( ina.size() != n ) stop("ina must have one label per training observation");
  if ( cons.size() != k ) stop("kappa and cons must have the same length");

  // Training data grouped by class, stored contiguously so the inner loop
  // walks memory sequentially: Yc holds the columns of class 0, then 1, ...
  std::vector<int> start(k + 1, 0);
  for (int l = 0; l < n; ++l) {
    const int g = ina[l];
    if ( g < 1 || g > k ) stop("labels in ina must lie in 1..k");
    ++start[g];
  }
  for (int j = 0; j < k; ++j) start[j + 1] += start[j];
  std::vector<double> Yc( (std::size_t)n * p );
  {
    std::vector<int> pos(start.begin(), start.end() - 1);
    const double* Y = y.begin();
    for (int l = 0; l < n; ++l) {
      const int g = ina[l] - 1;
      std::copy(Y + (std::size_t)l * p, Y + (std::size_t)(l + 1) * p,
                Yc.begin() + (std::size_t)pos[g]++ * p);
    }
  }

  // Plain C++ copies/pointers: no R API calls inside the parallel region.
  const std::vector<double> kap(kappa.begin(), kappa.end());
  const std::vector<double> cst(cons.begin(), cons.end());
  const double* Ynew = ynew.begin();
  const double* Yp = Yc.data();
  const int* st = start.data();
  NumericMatrix score(m, k);
  IntegerVector cl(m);
  double* S = score.begin();            // column-major: S[i + j * m]
  int* C = cl.begin();
  const int nt = n_threads(ncores);
  const int block = 4096;

  for (int b0 = 0; b0 < m; b0 += block) {
    checkUserInterrupt();
    const int b1 = std::min(m, b0 + block);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static)
#endif
    for (int i = b0; i < b1; ++i) {
      const double* yi = Ynew + (std::size_t)i * p;
      int best = 0;
      double bestval = R_NegInf;
      for (int j = 0; j < k; ++j) {
        double mx = R_NegInf, s = 0.0;
        for (int l = st[j]; l < st[j + 1]; ++l) {
          const double* yl = Yp + (std::size_t)l * p;
          double d = 0.0;
          for (int r = 0; r < p; ++r) d += yi[r] * yl[r];
          const double a = kap[j] * d;
          if ( a > mx ) {            // new maximum: rescale the running sum
            s = s * std::exp(mx - a) + 1.0;
            mx = a;
          } else {
            s += std::exp(a - mx);
          }
        }
        const double sc = cst[j] + mx + std::log(s);
        S[i + (std::size_t)j * m] = sc;
        if ( sc > bestval ) { bestval = sc;  best = j; }
      }
      C[i] = best + 1;
    }
  }
  return List::create(Named("class") = cl, Named("score") = score);
}
