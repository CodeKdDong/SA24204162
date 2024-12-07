#' Group Lasso Optimization
#'
#' @param X Design matrix, dimensions (n * T) x p
#' @param Y Response matrix, dimensions n x T
#' @param lambda Regularization parameter for L1 penalty
#' @param T Number of tasks
#' @param n Number of samples per task
#' @param p Number of features per task
#'
#' @return A matrix of coefficients (p x T).
#' @export
