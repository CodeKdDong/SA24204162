#' @title An estimation of intrinsic dimension of linear coefficients
#' @description An estimation of intrinsic dimension of linear coefficients
#' @param X Explation variable array with dimension n(sample size),
#' p(dimension),T(task size) 
#' @param Y Response variable with dimension n(sample size) and T(task size)
#' @param T task size for MTL 
#' @param T1 thresholding penalty parameter for sqrt((p+logT)/n),default 0.5
#' @param T2 thresholding penalty parameter for R*epsilon_bar(default 0.25)
#' @param epsilon_bar projection quantile level for Spectral 
#' Method(eg:0.1 represent projection to 0.9-quantile level)
#' @return intrinsic dimension estimation for linear coefficient
#' @examples
#' \dontrun{
#'     X <- array(runif(5 * 100 * 10), dim = c(100, 10, 5))
#'     Y <- array(runif(5 * 100), dim = c(100, 5))
#'     r_hat <- r_estimation(Y,X,5,0.1)
#'     print(r_hat)
#' }
#' @export


r_estimation <- function(Y,X,T,epsilon_bar,T1=0.5,T2=0.25){
  p <- dim(X)[2]
  n <- dim(X)[1]
  beta_hat <- single_task(Y,X,T)
  B_hat <- projection(T,beta_hat,X,epsilon_bar)
  svd_B <- svd(B_hat)
  singular_values <- svd_B$d  # Singular values of B
  # Step 2: Projection and concatenation
  beta_norms <- sapply(1:T, function(t) norm(beta_hat[[t]], type = "2"))  #         Calculate the 2-norms
  quantile_val <- quantile(beta_norms, 1 - epsilon_bar)  # Upper quantile of the   norms
  threshold <- T1 * sqrt((p + log(T)) / n) + T2 * quantile_val * sqrt(epsilon_bar)
  
  # Estimate the intrinsic dimension r
  r_hat <- max(which(singular_values >= threshold))
  
  # Output the estimated r
  return(r_hat)
}



