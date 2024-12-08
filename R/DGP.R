#' @title data generation process given some information using R
#' @description data generation process given some information using R
#' @param n sample size of our task
#' @param p sample dimension of our task
#' @param T task size for our MTL
#' @param r intrinsic dimension for linear coefficients
#' @param h similarity concentration level for representation matrix A
#' @param epsilon_bar projection quantile level for Spectral 
#' Method(eg:0.1 represent projection to 0.9-quantile level)
#' @return data contains sample and some parameter information
#' @examples
#' \dontrun{
#'     data <- generate_data(n, p, T, r, h, epsilon_bar)
#'     beta_true <- data$beta  # True beta values from the model
#'     beta_final_spectral <- Spectral_Method_kr(data$X, data$Y, 
#'     T, r, epsilon_bar, gamma = 1)
#' }
#' @export

generate_data <- function(n, p, T, r, h, epsilon_bar) {
  # Generate X: Input array (n, p, T)
  X <- array(rnorm(n * p * T), dim = c(n, p, T))
  
  # Generate C: Random p x r matrix with standard normal entries
  C <- matrix(rnorm(p * r), nrow = p, ncol = r)
  
  # Generate A_hat: Left singular matrix of C
  svd_C <- svd(C)
  A_hat <- svd_C$u[, 1:r]
  
  # Generate a(t) for each time t in [T]
  a <- runif(T, min = -h, max = h)
  
  # Generate A_t for each task, ensuring correct dimension
  A_t <- lapply(1:T, function(t) {
    # Perturb A_hat with a(t) by adding random perturbation
    perturbation <- a[t] * matrix(rnorm(p * r), nrow = p, ncol = r)
    A_t_val <- A_hat + perturbation
    return(A_t_val)
  })
  
  # Generate theta_star: Linear coefficients for each task t
  theta_star <- lapply(1:T, function(t) {
    runif(r, min = -2, max = 2)
  })
  
  
  # Generate Y: Response variable (n, T)
  Y <- matrix(0, nrow = n, ncol = T)
  for (t in 1:T) {
    Y[, t] <- X[, , t] %*% A_t[[t]] %*% theta_star[[t]] + rnorm(n)
  }
  
  beta <- list()
  for(t in 1:T){
    beta[[t]] <- A_t[[t]] %*% theta_star[[t]]
  }
  
  
  # Return generated data
  return(list(X = X, Y = Y, A_hat = A_hat, theta_star = theta_star,beta=beta))
}



