single_task <- function(Y,X,T){
  # Initialize the storage for estimated betas
  beta_hat <- list()
  for (t in 1:T) {
    beta_hat[[t]] <- lm(Y[, t] ~ X[, , t] - 1)$coefficients  # Linear regression (no intercept)
  }
  return(beta_hat)
}