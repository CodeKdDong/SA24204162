#' @title Singular Value Decomposition (SVD) for Matrix Decomposition
#' 
#' @description This function performs Singular Value Decomposition (SVD) on a given matrix \code{B_hat}.
#' It then extracts the first \code{r} columns of the left singular vector matrix \code{U}, 
#' which corresponds to the intrinsic dimension of the matrix.
#'
#' @param B_hat A matrix of size \code{n} x \code{p} to undergo singular value decomposition.
#' @param r The intrinsic dimension for the decomposition, i.e., the number of left singular vectors to extract.
#' 
#' @return A matrix \code{A_hat} containing the first \code{r} columns of the left singular vector matrix \code{U}
#'         from the SVD of \code{B_hat}.
#'
#' @examples
#' \dontrun{
#'     # Example usage:
#'     A_hat <- SVD(B_hat = some_matrix, r = 5)
#'     # A_hat will contain the first 5 left singular vectors of B_hat
#' }
#' 
#' @export
SVD <- function(B_hat, r) {
  # Step 3: Singular Value Decomposition (SVD)
  svd_B <- svd(B_hat)
  U_hat <- svd_B$u  # U matrix from SVD
  A_hat <- U_hat[, 1:r]  # Take the first r columns of U (intrinsic dimension)
  
  return(A_hat)
}
