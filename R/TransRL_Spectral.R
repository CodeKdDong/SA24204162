#' @title Spectral Method for Representation Multitask Learning(RMTL) 
#' applies to new transfer learning task
#' @description Spectral Method for Representation Multitask Learning
#'  applies to new transfer learning task
#' @param X_source Explation source variable array with dimension 
#' n(sample size),p(dimension),T(task size) 
#' @param Y_source Response source variable with dimension n(sample size) 
#' and T(task size)
#' @param X_target New explation variable with dimension n0(sample size) and p(variable dimension) 
#' @param Y_target New response variable with dimension n0
#' @param T task size for MTL 
#' @param epsilon_bar projection quantile level for Spectral Method
#' @param gamma regularization parameter for bias regularization step
#' @param T1 thresholding penalty parameter for sqrt((p+logT)/n),default 0.5
#' @param T2 thresholding penalty parameter for R*epsilon_bar(default 0.25)
#' @return linear coefficients for new task
#' @examples
#' \dontrun{
#'     X_source <- array(runif(5 * 100 * 10), dim = c(100, 10, 5))
#'     Y_source <- array(runif(5 * 100), dim = c(100, 5))
#'     Y_target  <- runif(n)
#'     X_target <- matrix(runif(n*p),ncol = p,nrow = n)
#'     beta0 <- TransRL_Spectral(X_source,Y_source,X_target,Y_target,T,epsilon_bar,
#'     gamma,T1=0.5,T2=0.25)
#'     print(beta0)
#' }
#' @export


TransRL_Spectral <- function(X_source,Y_source,X_target,Y_target,T,epsilon_bar,gamma,T1=0.5,T2=0.25){
  
  n0 <- dim(X_target)[1]
  r <- r_estimation(Y_source,X_source,T,epsilon_bar,T1,T2)
  theta_hat <- list()
  beta_hat <- single_task(Y_source,X_source,T)
  B_hat <- projection(T,beta_hat,X_source,epsilon_bar)
  A_hat <- SVD(B_hat,r)
  
  # Define the objective function for θ(t)
  objective_function_theta <- function(theta) {
    # Squared loss
    
    product <- as.vector(A_hat %*% theta)
    f_t_theta <- sum((Y_target - X_target %*% product)^2)
    return(f_t_theta)
  }
  # Initial guess for θ(t)
  theta_initial <- rep(0, r)  # Start with zero vector of dimension r
  
  # Use optim to find the optimal θ(t)
  theta_hat0 <- optim(theta_initial, objective_function_theta)$par
  
  
  # Step 5: Biased regularization
  beta_0 <- list()
  
  # Initial guess for beta_t
  beta0_hat <- lm(Y_target ~ X_target - 1)$coefficients
  beta_0 <- as.vector(beta0_hat)  # Ensure it's a column vector
  
  # Define the objective function for β(t)
  objective_function_beta <- function(beta) {
    # Squared loss
    loss <- sum((Y_target - X_target %*% beta)^2)
    
    # Regularization term
    regularization <- gamma / sqrt(n0) * sum((beta - A_hat %*% theta_hat0)^2)
    
    # Return the combined objective
    return(loss + regularization)
  }
  
  # Use optim to minimize the objective function for β(t)
  beta_final <- optim(beta_0, objective_function_beta)$par
  
  return(beta_final)
}















