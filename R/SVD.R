SVD <- function(B_hat,r){
  # Step 3: Singular Value Decomposition (SVD)
  svd_B <- svd(B_hat)
  U_hat <- svd_B$u  # U matrix from SVD
  A_hat <- U_hat[, 1:r]  # Take the first r columns of U (intrinsic dimension)
  return(A_hat)
}