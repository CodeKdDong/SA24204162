#include <RcppArmadillo.h>
using namespace Rcpp;
//' @title Group Lasso Optimization Algorithm
//' 
//' @description This function implements the Group Lasso optimization algorithm to learn coefficients for a multi-task 
//' regression problem. It uses gradient descent with L1 regularization on each group of coefficients (each task).
//' The optimization minimizes the residual sum of squares subject to the group lasso penalty.
//'
//' @param X A numeric matrix of size \code{n * p * T}, where \code{n} is the number of samples, 
//'        \code{p} is the number of features, and \code{T} is the number of tasks.
//' @param Y A numeric matrix of size \code{n * T}, where \code{n} is the number of samples and 
//'        \code{T} is the number of tasks. Each column represents the response for a particular task.
//' @param lambda A regularization parameter that controls the strength of the group lasso penalty.
//' @param T The number of tasks. This defines the number of columns in the output \code{beta} matrix.
//' @param n The number of samples.
//' @param p The number of features in the input matrix.
//' 
//' @return A numeric matrix \code{beta} of size \code{p * T} containing the estimated coefficients for each task.
//'         Each column represents the estimated coefficients for a particular task, and each row corresponds 
//'         to the coefficients of a specific feature.
//'
//' @examples
//' \dontrun{
//'     // Example usage:
//'     beta_hat <- groupLasso(X = X_data, Y = Y_data, lambda = 0.1, T = 5, n = 100, p = 20)
//'     // beta_hat will contain the optimized coefficients for each task
//' }
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix groupLasso(NumericMatrix X, NumericMatrix Y, double lambda, int T, int n, int p) {
  NumericMatrix beta(p, T);  // Coefficients for each task
  NumericMatrix grad(p, T); // Gradient matrix
  
  double tolerance = 1e-6;  // Convergence tolerance
  double step_size = 0.01;  // Step size for gradient descent
  int max_iter = 1000;      // Maximum number of iterations
  
  // Define customSign function as a lambda
  auto customSign = [](double x) -> double {
    if (x > 0) return 1.0;
    if (x < 0) return -1.0;
    return 0.0;
  };
  
  
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
