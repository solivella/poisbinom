#include <Rcpp.h>
#include <R_ext/Utils.h> // for findInterval();
#include <fftw3.h>


void dft_pmf(fftw_complex* out, int m,  Rcpp::NumericVector& pp);


/*
  Probability mass function
*/
// [[Rcpp::export]]
Rcpp::NumericVector dpoisbinom(Rcpp::IntegerVector& x,
			       Rcpp::NumericVector& pp,
			       bool log_d = false)
{
  //Check bounds
  if( Rcpp::is_true(Rcpp::any(pp > 1)) || Rcpp::is_true(Rcpp::any(pp < 0)) ){
    Rcpp::stop("Values in pp must be between 0 and 1.");
  }

  int m = pp.size() + 1;
  int nn = x.size();
  fftw_complex* out;
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m);

  //dft
  dft_pmf(out, m, pp);

  //form return object
  Rcpp::NumericVector res(nn);
  double scale = 1.0 / m;
  int kk;
  for(int k = 0; k < nn; ++k)
    {
      kk = x[k];
      res[k] = out[kk][0] * scale;
    }

  //Destroy dft object
  fftw_free(out);

  if(log_d)
    return(log(res));
  else
    return(res);
}


/*
  Raw function that computes the cumulative probabilities for first 'max_q' values
  of the support, i.e., 0, 1, ..., max_q - 1.
*/
Rcpp::NumericVector ppoisbinom_raw(int max_q,
				   Rcpp::NumericVector& pp)
{

  int m = pp.size() + 1;

  fftw_complex* out;
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m);

  //dft
  dft_pmf(out, m, pp);

  // form cumulative probabilities
  Rcpp::NumericVector csum(max_q);
  double scale = 1.0 / m;
  csum[0] = out[0][0] * scale;
  int k = 1;
  do
    {
      csum[k] = out[k][0] * scale + csum[k - 1];
      ++k;
    }
  while(k < max_q);

  //Destroy dft object
  fftw_free(out);

  return csum;
}


/*
  C++ wrapper for the cumulative probability function
*/
// [[Rcpp::export]]
Rcpp::NumericVector ppoisbinom(Rcpp::IntegerVector& q,
			       Rcpp::NumericVector& pp,
			       bool lower_tail = true,
			       bool log_p = false)
{
  //Check bounds
  if( Rcpp::is_true(Rcpp::any(pp > 1)) || Rcpp::is_true(Rcpp::any(pp < 0)) ){
    Rcpp::stop("Values in pp must be between 0 and 1.");
  }

  int max_q = max(q) + 1; // maximum of quantiles plus one

  Rcpp::NumericVector csum = ppoisbinom_raw(max_q, pp);

  //form return object
  int nn = q.size();
  Rcpp::NumericVector res(nn);
  int kk;
  for(int k = 0; k < nn; ++k)
    {
      kk = q[k];
      res[k] = csum[kk];
    }

  if(!lower_tail)
    res = 1.0 - res;

  if(log_p)
    return(log(res));
  else
    return(res);
}


/*
  Subroutine to find interval in cdf
*/
Rcpp::IntegerVector find_from_cdf(Rcpp::NumericVector& csum,
				  Rcpp::NumericVector& s_invec,
				  Rcpp::IntegerVector& order,
				  int n,
				  int t_res)
{
  Rcpp::IntegerVector res(n);
  int kk;
  int flag;
  for(int k = 0; k < n; ++k)
    {
      t_res = findInterval(&csum[0], csum.size(), s_invec[k], FALSE, FALSE, t_res, &flag);
      kk = order[k];
      res[kk-1] = t_res;
    }

  return res;

}


/*
  Quantile function
*/
// [[Rcpp::export]]
Rcpp::IntegerVector qpoisbinom(Rcpp::NumericVector& p,
			       Rcpp::NumericVector& pp,
			       bool lower_tail = true,
			       bool log_p = false)
{
  //Check bounds
  if( Rcpp::is_true(Rcpp::any(pp > 1)) || Rcpp::is_true(Rcpp::any(pp < 0)) ){
    Rcpp::stop("Values in pp must be between 0 and 1.");
  }
  if( Rcpp::is_true(Rcpp::any(p > 1)) || Rcpp::is_true(Rcpp::any(p < 0)) ){
    Rcpp::stop("Values in p must be between 0 and 1.");
  }

  if (log_p) p = exp(p);

  Rcpp::NumericVector csum = ppoisbinom_raw(pp.size() + 1, pp);

  //sort keeping track of original order
  int nn = p.size();
  Rcpp::NumericVector s_invec = Rcpp::clone(p).sort();
  Rcpp::IntegerVector order = Rcpp::match(s_invec, p);

  //find interval on sorted vector, and form return object
  int t_res = std::floor(Rcpp::sum(pp));
  Rcpp::IntegerVector res = find_from_cdf(csum, s_invec, order, nn, t_res);

  return(res);

}


/*
  Random number generation
*/
// [[Rcpp::export]]
Rcpp::IntegerVector rpoisbinom(int n,
			       Rcpp::NumericVector& pp)
{
  //Check bounds
  if( Rcpp::is_true(Rcpp::any(pp > 1)) || Rcpp::is_true(Rcpp::any(pp < 0)) ){
    Rcpp::stop("Values in pp must be between 0 and 1.");
  }
  // generate random number from uniform(0, 1) and sort the output
  Rcpp::NumericVector u = Rcpp::runif(n);
  Rcpp::NumericVector sorted_u = Rcpp::clone(u).sort();
  Rcpp::IntegerVector order = Rcpp::match(sorted_u, u);

  // obtain cdf
  Rcpp::NumericVector csum = ppoisbinom_raw(pp.size() + 1, pp);

  // inverse-cdf method: find interval on sorted vector, and form return object
  int t_res = std::floor(Rcpp::sum(pp));
  Rcpp::IntegerVector res = find_from_cdf(csum, sorted_u, order, n, t_res);

  return(res);

}

/*
  Probability mass function using DFT
*/
void dft_pmf(fftw_complex* out, int m,  Rcpp::NumericVector& pp)
{
  int n = m - 1;
  fftw_complex* in;
  fftw_plan p;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m);
  std::complex<double> C(0.0,2);
  C = exp(C * 3.1415926535897 / ((double)m));
  double C_real = C.real();
  double C_imag = C.imag();


  //build input vector
  double tmp_real;
  double tmp_imag;
  std::complex<double> temp;
  std::complex<double> f(1.,0.0);
  in[0][0] = 1.0;
  in[0][1] = 0.0;


  int halfn = n / 2 + 1;
  for (int i = 1; i <= halfn; ++i){
    temp = 1.;
    tmp_real = f.real();
    tmp_imag = f.imag();
    f.real(tmp_real * C_real - tmp_imag * C_imag);
    f.imag(tmp_imag * C_real + tmp_real * C_imag);
    for(int j = 0; j < n; ++j){
      temp *= (1. + (f - 1.) * pp[j]);
    }

    in[i][0] = temp.real();
    in[i][1] = temp.imag();
    in[m-i][0] = temp.real();
    in[m-i][1] = - temp.imag();
  }


  //dft
  p = fftw_plan_dft_1d(m, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  fftw_free(in);
}
