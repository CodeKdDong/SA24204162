#' @title Projection onto the ℓ2-ball
#' @description This function performs a projection of a sequence of vectors (betas) onto the ℓ2-ball 
#' based on the upper quantile of the norms of the vectors. The projection scales vectors
#' with a norm greater than the quantile value to fit within the ℓ2-ball of radius defined 
#' by the quantile value.
#'
#' @param T An integer representing the number of iterations or the length of `beta_hat`.
#' @param beta_hat A list of vectors (length T), where each element is a vector of coefficients.
#' @param X A matrix of input data, where rows represent observations and columns represent features.
#' @param epsilon_bar A numeric value (between 0 and 1) used to determine the quantile value 
#'        for projection based on the upper quantile of the norms.
#' 
#' @return A matrix of projected coefficients, with each column representing the projected 
#'         vector corresponding to each element in `beta_hat`.
#' 
#' @examples{
#' \dontrun{
#' # Example usage:
#' projection_result <- projection(100, beta_hat_list, X_matrix, 0.05)
#' # Visualizing the result
#' plot(projection_result[, 1], type = 'l')
#' plot(projection_result[, 2], type = 'l')
#' }
#' }
#' @export
projection <- function(T, beta_hat, X, epsilon_bar) {
  p <- dim(X)[2]
  n <- dim(X)[1]
  
  # Step 2: Projection and concatenation
  beta_norms <- sapply(1:T, function(t) norm(beta_hat[[t]], type = "2"))  # Calculate the 2-norms
  quantile_val <- quantile(beta_norms, 1 - epsilon_bar)  # Upper quantile of the norms
  
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
