#' @title Spectral Method for Representation Multitask Learning(RMTL) 
#' when unknow intrinsic dimension r
#' @description Spectral Method for Representation Multitask Learning
#'  when unknow intrinsic dimension r
#' @param X Explation variable array with dimension 
#' n(sample size),p(dimension),T(task size) 
#' @param Y Response variable with dimension 
#' n(sample size) and T(task size)
#' @param T task size for MTL 
#' @param epsilon_bar projection quantile level for Spectral Method
#' @param gamma regularization parameter for bias regularization step
#' @param T1 thresholding penalty parameter for sqrt((p+logT)/n),default 0.5
#' @param T2 thresholding penalty parameter for R*epsilon_bar(default 0.25)
#' @return linear coefficients for T tasks using RMTL method
#' @examples
#' \dontrun{
#'     X <- array(runif(5 * 100 * 10), dim = c(100, 10, 5))
#'     Y <- array(runif(5 * 100), dim = c(100, 5))
#'     beta_final <- Spectral_Method_ur(X,Y,5,0.1,1)
#'     print(beta_final)
#' }
#' @export



Spectral_Method_ur <- function(X,Y,T,epsilon_bar,gamma,T1=0.5,T2=0.25){
  # Step 4: θ(t) optimization
  n <- dim(X)[1]
  r <- r_estimation(Y,X,T,epsilon_bar,T1,T2)
  theta_hat <- list()
  beta_hat <- single_task(Y,X,T)
  B_hat <- projection(T,beta_hat,X,epsilon_bar)
  A_hat <- SVD(B_hat,r)
  for (t in 1:T) {
    # Define the objective function for θ(t)
    objective_function_theta <- function(theta) {
      # Squared loss
      f_t_theta <- sum((Y[, t] - X[, , t] %*% (A_hat %*% theta))^2)
      return(f_t_theta)
    }
    
    # Initial guess for θ(t)
    theta_initial <- rep(0, r)  # Start with zero vector of dimension r
    
    # Use optim to find the optimal θ(t)
    theta_hat[[t]] <- optim(theta_initial, objective_function_theta)$par
  }
  
  # Step 5: Biased regularization
  beta_final <- list()
  
  for (t in 1:T) {
    # Initial guess for beta_t
    beta_t <- as.vector(beta_hat[[t]])  # Ensure it's a column vector
    
    # Define the objective function for β(t)
    objective_function_beta <- function(beta) {
      # Squared loss
      loss <- sum((Y[, t] - X[, , t] %*% beta)^2)
      
      # Regularization term
      regularization <- gamma / sqrt(n) * sum((beta - A_hat %*% theta_hat[[t]])^2)
      
      # Return the combined objective
      return(loss + regularization)
    }
    
    # Use optim to minimize the objective function for β(t)
    beta_final[[t]] <- optim(beta_t, objective_function_beta)$par
  }
  return(beta_final)
}


