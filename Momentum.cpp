#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Momentum(): computes the (scaled) momentum statistic using three arguments.
//    Nsim : the number of posterior simulations we compute each statistic for
//    N: the number of observations in the random walk
//    theSim : a matrix of (simulated) posterior values 

// [[Rcpp::export]]
NumericVector Momentum(double Nsim, double N, NumericMatrix theSim) {
  // initialize variables 
  NumericVector posterior_samples(Nsim);    // one entry for each posterior sample
  
  // calculate the number of positive differences 
  for (int ii = 0; ii < Nsim; ii++) { // ii = posterior sample 
    double total_contribution = 0; 
    double counter = 0;
    for (int jj = 0; jj < N-1; jj++) { 
      for (int kk = jj+1; kk < N; kk++) { 
        if ( theSim(kk,ii) - theSim(jj,ii) > 0 ) { total_contribution += 1; }
        if ( theSim(kk,ii) - theSim(jj,ii) == 0 ) { total_contribution += 0.5; }
        if ( theSim(kk,ii) - theSim(jj,ii) < 0 ) { total_contribution += 0; }
        counter += 1; 
      }
    }
    posterior_samples[ii] = total_contribution / counter;
  }
  return posterior_samples; 
}
