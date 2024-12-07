projection <- function(T,beta_hat,X,epsilon_bar){
  p <- dim(X)[2]
  n <- dim(X)[1]
  # Step 2: Projection and concatenation
  beta_norms <- sapply(1:T, function(t) norm(beta_hat[[t]], type = "2"))  #         Calculate the 2-norms
  quantile_val <- quantile(beta_norms, 1 - epsilon_bar)  # Upper quantile of the   norms
  
  # Projection operator QR onto the ℓ2-ball
  B_hat <- matrix(0, nrow = p, ncol = T)  # Matrix to store projected betas
  
  for (t in 1:T) {
    beta_t <- beta_hat[[t]]
    norm_beta_t <- norm(beta_t, type = "2")
    
    # Perform decomposition to project onto ℓ2-ball of radius quantile_val
    if (norm_beta_t > quantile_val) {
      # Scale the vector to fit inside the ℓ2-ball
      scaling_factor <- quantile_val / norm_beta_t
      B_hat[, t] <- scaling_factor * beta_t
    } else {
      B_hat[, t] <- beta_t  # No scaling needed
    }
  }
  return(B_hat)
}
