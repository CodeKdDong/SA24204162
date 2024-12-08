#include <RcppArmadillo.h>
#include <random>
using namespace Rcpp;

//' @title Gibbs sampling algorithm for sampling from the joint distribution of 
//' two variables x and y. 
//' @description The function performs sampling from the conditionaldistributions p(x | y) and p(y | x) iteratively.
//' @param n Integer. The total number of trials for the binomial distribution in each step.
//' @param a Double. The first shape parameter for the Beta distribution.
//' @param b Double. The second shape parameter for the Beta distribution.
//' @param num_samples Integer. The number of posterior samples to draw after burn-in.
//' @param burn_in Integer. The number of initial iterations to discard as burn-in.
//' @return NumericMatrix A matrix of size (num_samples, 2), where each row contains
//' the sampled values of x and y.
//' @examples
//' \dontrun{
//' samples <- gibbs_sampler(10, 2, 2, 1000, 200)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbs_sampler(int n, double a, double b, int num_samples, int burn_in) {
  // Initialize the matrix to store the sampled values of x and y
  NumericMatrix samples(num_samples, 2);
  
  // Initialize the variables x and y with arbitrary starting values
  int x = n / 2;  // Choose an initial value for x
  double y = 0.5; // Choose an initial value for y within the range (0, 1)
  
  // Set up the random number generator
  std::random_device rd;
  std::mt19937 gen(rd());
  
  // Main loop for Gibbs sampling
  for (int i = 0; i < num_samples + burn_in; ++i) {
    // Step 1: Sample x from the conditional distribution x | y ~ Binomial(n, y)
    x = R::rbinom(n, y);
    
    // Step 2: Sample y from the conditional distribution y | x ~ Beta(x + a, n - x + b)
    y = R::rbeta(x + a, n - x + b);
    
    // After the burn-in period, store the sampled values of x and y
    if (i >= burn_in) {
      samples(i - burn_in, 0) = x;
      samples(i - burn_in, 1) = y;
    }
  }
  
  return samples;
}
