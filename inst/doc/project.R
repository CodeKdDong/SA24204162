## -----------------------------------------------------------------------------
library(SA24204162)
library(Rcpp)
library(RcppArmadillo)

## -----------------------------------------------------------------------------
dir_cpp <- '../src/'

# Can create source file in Rstudio
sourceCpp(paste0(dir_cpp,"Grouplasso.cpp"))


## -----------------------------------------------------------------------------
# parametric settings
n <- 150
p <- 80
T <- 10
r <- 10
h <- 0.5
epsilon_bar <- 0.1
gamma = 1
data <- generate_data(n, p, T, r, h, epsilon_bar)
data_target <- generate_data(n, p, 1, r, h, epsilon_bar)

# Spectral Method
beta_final_spectral <- Spectral_Method_kr(data$X, data$Y, T, r, epsilon_bar, gamma = 1)


# 调用 Rcpp 函数计算 beta
# 将 X 重塑为二维矩阵 (n * T, p)
X_matrix <- matrix(NA, nrow = n * T, ncol = p)
for (t in 1:T) {
  X_matrix[((t - 1) * n + 1):(t * n), ] <- data$X[, , t]  # 这里要确保X_matrix正确赋值
}


# Group Lasso (using the custom C++ function)
lambda <- 0.1  # Regularization parameter for Group Lasso
beta_group_lasso <- groupLasso(X_matrix, data$Y, lambda, T, n, p)

# Compute L2 Loss for Spectral Method
beta_true <- data$beta  # True beta values from the model



## -----------------------------------------------------------------------------
r_hat <- r_estimation(data$Y,data$X,T,epsilon_bar)
r_hat


## -----------------------------------------------------------------------------
l2_loss <- function(x,y){
  return(min(mean((x-y)^2),1))
}

## -----------------------------------------------------------------------------
loss_beta <- 0
for(t in 1:T){
  beta_true[[t]] <- as.numeric(beta_true[[t]])
  loss_beta <- l2_loss(beta_true[[t]],beta_final_spectral[[t]]) + loss_beta
}
loss_beta/T

## -----------------------------------------------------------------------------
loss_lasso <- 0
for(t in 1:T){
  beta_true[[t]] <- as.numeric(beta_true[[t]])
  loss_lasso <- l2_loss(beta_true[[t]],beta_group_lasso[,t]) + loss_lasso
}
loss_lasso/T

## -----------------------------------------------------------------------------
data <- generate_data(100, 10, T, r, h, epsilon_bar)
data_target <- generate_data(100, 10, 1, r, h, epsilon_bar)

## -----------------------------------------------------------------------------
X_target <- as.matrix(data_target$X[,,1])

Y_target <- as.vector(data_target$Y)
beta0 <- TransRL_Spectral(data$X,data$Y,X_target,Y_target,T,epsilon_bar,gamma,T1=0.5,T2=0.25)
l2_loss(beta0,data_target$beta[[1]])

## -----------------------------------------------------------------------------
beta_ls <- lm(Y_target ~ X_target - 1)$coefficients  
l2_loss(beta_ls,data_target$beta[[1]])

## -----------------------------------------------------------------------------
data <- generate_data(100, 80, T, r, h, epsilon_bar)
data_target <- generate_data(100, 80, 1, r, h, epsilon_bar)

## -----------------------------------------------------------------------------
X_target <- as.matrix(data_target$X[,,1])

Y_target <- as.vector(data_target$Y)
beta0 <- TransRL_Spectral(data$X,data$Y,X_target,Y_target,T,epsilon_bar,gamma,T1=0.5,T2=0.25)
l2_loss(beta0,data_target$beta[[1]])

## -----------------------------------------------------------------------------
beta_ls <- lm(Y_target ~ X_target - 1)$coefficients  
l2_loss(beta_ls,data_target$beta[[1]])

