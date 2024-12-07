#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath> // For sqrt function
using namespace Rcpp;

// [[Rcpp::export]]
double customSign(double x) {
  if (x > 0) return 1.0;
  if (x < 0) return -1.0;
  return 0.0;
}

// [[Rcpp::export]]
/*
Implements the group Lasso optimization algorithm.
*/
NumericMatrix groupLasso(NumericMatrix X, NumericMatrix Y, double lambda, int T, int n, int p) {
  NumericMatrix beta(p, T);  // Coefficients for each task
  NumericMatrix grad(p, T); // Gradient matrix
  
  double tolerance = 1e-6;  // Convergence tolerance
  double step_size = 0.01;  // Step size for gradient descent
  int max_iter = 1000;      // Maximum number of iterations
  
  // Initialize beta matrix to 0
  for (int t = 0; t < T; t++) {
    for (int j = 0; j < p; j++) {
      beta(j, t) = 0.0;
    }
  }
  
  // Main optimization loop
  for (int iter = 0; iter < max_iter; iter++) {
    // Compute gradients
    for (int t = 0; t < T; t++) {
      for (int j = 0; j < p; j++) {
        double grad_val = 0.0;
        
        // Accumulate gradient over samples
        for (int i = 0; i < n; i++) {
          grad_val += X(i + t * n, j) * (X(i + t * n, j) * beta(j, t) - Y(i, t));
        }
        
        // Add L1 regularization term
        grad(j, t) = grad_val + lambda * customSign(beta(j, t));
      }
    }
    
    // Update beta
    for (int t = 0; t < T; t++) {
      for (int j = 0; j < p; j++) {
        beta(j, t) -= step_size * grad(j, t);
      }
    }
    
    // Check convergence: calculate L2 norm of gradient
    double norm_grad = 0.0;
    for (int t = 0; t < T; t++) {
      for (int j = 0; j < p; j++) {
        norm_grad += grad(j, t) * grad(j, t);
      }
    }
    norm_grad = sqrt(norm_grad);
    if (norm_grad < tolerance) { // Stop if gradient norm is below tolerance
      break;
    }
  }
  
  // Return the optimized coefficient matrix
  return beta;
}
