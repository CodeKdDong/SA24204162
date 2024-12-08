#' @title Single Task Linear Regression
#' @description This function performs linear regression for each task
#'  in a multi-task setting.
#' @param Y A matrix of response variables of size \code{n} x \code{T}, where each column represents 
#'        the responses for a particular task.
#' @param X A 3D array of input features of size \code{n} x \code{p} x \code{T}, where each slice corresponds 
#'        to the feature matrix for a particular task.
#' @param T The number of tasks. This determines how many times the regression will be performed.
#'
#' @return A list containing the estimated coefficients for each task
#' @examples
#' \dontrun{
#'     # Example usage:
#'     beta_hat <- single_task(Y = data$Y, X = data$X, T = 10)
#'     # beta_hat contains the estimated coefficients for each task
#' }
#' 
#' @export
single_task <- function(Y, X, T) {
  # Initialize the storage for estimated betas
  beta_hat <- list()
  
  for (t in 1:T) {
    # Perform linear regression (no intercept) for each task t
    beta_hat[[t]] <- lm(Y[, t] ~ X[, , t] - 1)$coefficients
  }
  
  return(beta_hat)
}
