## -----------------------------------------------------------------------------
Generate_Rayleign <- function(N,seed,sig){
set.seed(seed)
u <- runif(N)
x <- sqrt(-2*sig^2*log(1-u))
hist(x,probability = TRUE,main = paste('Empirical Density Graph of Rayleign Distribution--',sig))
y <- seq(0,8,0.01)
lines(y,y/(sig[1]^2) * exp(-y^2/(2*sig^2)) )
}

## -----------------------------------------------------------------------------
seed = 123
N = 10000
sig <- seq(0.1,2,0.4)
Generate_Rayleign(N,seed,sig[1])

## -----------------------------------------------------------------------------
Generate_Rayleign(N,seed,sig[2])

## -----------------------------------------------------------------------------
Generate_Rayleign(N,seed,sig[3])

## -----------------------------------------------------------------------------
Generate_Rayleign(N,seed,sig[4])

## -----------------------------------------------------------------------------
Generate_Rayleign(N,seed,sig[5])

## -----------------------------------------------------------------------------
theory_density <- function(p1,y){
  return(p1/sqrt(2*pi)*exp(-((y)^2)/(2))+(1-p1)/sqrt(2*pi)*exp(-((y-3)^2)/(2)))
}

## -----------------------------------------------------------------------------
set.seed(3)
N = 1000
p1 <- 0.75
x1 <- rnorm(N,0,1)
x2 <- rnorm(N,3,1)
p <- sample(c(1,0),1,prob = c(p1,1-p1))
x <- p*x1+(1-p)*x2
hist(x,probability = TRUE,main = paste('Mix-Gaussian Distribution--',p1))
y <- seq(-4,4,0.01)
lines(y,theory_density(p1,y))

## -----------------------------------------------------------------------------
Generate_MixGaussian<- function(p1,N,seed){
set.seed(seed)
x1 <- rnorm(N,0,1)
x2 <- rnorm(N,3,1)
p <- sample(c(1,0),N,prob = c(p1,1-p1),replace = TRUE)
x <- p*x1+(1-p)*x2
hist(x,probability = TRUE,main = paste('Mix-Gaussian Distribution--',p1))
y <- seq(-8,8,0.01)
lines(y,theory_density(p1,y))
}

## -----------------------------------------------------------------------------
p <- seq(0.1,0.9,0.1)

## -----------------------------------------------------------------------------
for(i in 1:9){
Generate_MixGaussian(p[i],N,123)  
}


## -----------------------------------------------------------------------------
# Function to simulate the compound Poisson-Gamma process
Generate_compound_poisson_gamma <- function(lambda, k, theta, t, N) {
results <- numeric(N)
  
for (i in 1:N) {
# Simulate N(t), the number of events in time t, from Poisson distribution
    N_t <- rpois(1, lambda * t)
    
    # Simulate Y_i from Gamma distribution, shape = k, rate = theta
    if (N_t > 0) {
      Y <- rgamma(N_t, shape = k, rate = theta)
      results[i] <- sum(Y)  # X(t) = sum of the Y's
    } else {
      results[i] <- 0  # If N(t) = 0, X(t) = 0
    }
  }
  
  return(results)
}

## -----------------------------------------------------------------------------
#parameters
t <- 10# Time point for X(t)
lambda <- c(1,3,5,10,20)# Rate of Poisson process
k <- c(1,5,10,20,100)# Shape parameter of Gamma distribution
theta <- c(1/2,1,3,10,20)# Rate parameter of Gamma distribution
N <- 10000# Number of simulations

## -----------------------------------------------------------------------------
#Simalation
set.seed(123)
cat("Parameters Combination","Mean(Predict)","Mean(True)","Variance(Predict)","Variance(True)","Mean Error","Variance Error","\n")
for(i in 1:length(lambda)){
  for(j in 1:length(k)){
    for(l in 1:length(theta)){
      xt <- numeric(N)
      xt <-Generate_compound_poisson_gamma(lambda[i],k[j], theta[l], t, N)
      mean_hat <- mean(xt)
      var_hat <- var(xt)
      mean_true <- lambda[i] * t *(k[j]/theta[l])
      var_true <- lambda[i] * t *(((k[j]+k[j]^2)/(theta[l]^2)))
      value <- c(mean_hat,mean_true,var_hat,var_true)
      relative_error<-c(abs((mean_hat-mean_true)/mean_true),abs((var_hat-var_true)/var_true))
      cat("(",lambda[i],",",k[j],",",theta[l],"):",value,relative_error,"\n")
    }
  }
}


## ----include=FALSE------------------------------------------------------------
#load packages
library(microbenchmark)
library(ggplot2)

## -----------------------------------------------------------------------------
#Quicksort Alogorithm Realization
Quicksort <- function(a){
  if (length(a)<=1) {
      return(a)
  }
  #select a bechmark element
  pivot <- a[1]
  left <- a[a<pivot]
  right <- a[a>pivot]
  #recursively use Quicksort function
  return(c(Quicksort(left),pivot,Quicksort(right)))
}

## -----------------------------------------------------------------------------
#Test on Quicksort function
set.seed(124)
v <- c(2,3,4,6,3,-1,4,9)
Quicksort(v)

## -----------------------------------------------------------------------------
#Construct a function to compute the compution time
measure_sort_time <- function(n,simulations = 100){
  times <- replicate(simulations,{
    # Generate a random permutation of numbers from 1 to n
    perm <- sample(1:n)
    # Measure the sorting time via Quicksort
    start_time <- Sys.time()
    Quicksort(perm)
    end_time <- Sys.time()
    time <- end_time - start_time
    time
  })
  return(mean(times))
}

## -----------------------------------------------------------------------------
n_values <- c(1,2,4,6,8)*10^4
# Calculate average sorting times
a_n <- sapply(n_values,measure_sort_time)
#Calculate tn
t_n <- n_values * log(n_values)

## -----------------------------------------------------------------------------
df <- data.frame(
  t_n = t_n, 
  a_n = a_n
)
fit <- lm(a_n~t_n,data= df)
summary(fit)
#Regress Plot
# Plot scatter plot and regression line
ggplot(df, aes(x = t_n, y = a_n)) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue") +
  labs(title = "Regression of Average Sorting Time on t_n",
       x = "t_n = n log(n)",
       y = "Average Sorting Time (a_n)") +
  theme_minimal()

## -----------------------------------------------------------------------------
#A function of estimate for Beta(3,3)
BetaF_estimate <- function(x,seed=123,N=10000){
  set.seed(seed)
  X <- numeric(N)
  X <- rbeta(N,3,3)
  return(sum(X<x)/N)
}

## -----------------------------------------------------------------------------
#Get an estimate via Monte Carlo
x_values <- seq(0.1,0.9,0.1)
F_estimate <- sapply(x_values,BetaF_estimate)
#Generate value via pbeta function
F_theory <- pbeta(x_values,3,3)
#Compare the theory and simulation result
plot(x_values,y = F_estimate,type = 'l',col = 'red',
     main = 'cdf value of Beta(3,3)',xlab = 'x',ylab = 'F(x)')
lines(x_values,y=F_theory,col = 'blue')
legend('topleft',col = c('red','blue'),lty = c(1,1),
       lwd = c(2,2),legend = c('Estimate','Theory'))

## -----------------------------------------------------------------------------
#density function of rayleigh 
f <- function(sig,x){
  return(x/(sig^2)*exp(-x^2/(2*sig^2)))
}


## -----------------------------------------------------------------------------
#Generate Rayleign distribution via Inverse transformation
Generate_Rayleign <- function(sig,u,seed = 123){
set.seed(seed)
u1 <- u[1:(length(u)/2)]
u2 <- u[(length(u)/2+1):(length(u))]
x1 <- sqrt(-2*sig^2*log(1-u1))
x2 <- sqrt(-2*sig^2*log(1-u2))
x <- (as.vector(x1) + as.vector(x2))/2
return(x)
}

## -----------------------------------------------------------------------------
#antithetic variables
N = 10000
U1 <- runif(N/2,0,1)
u <- c(U1,1-U1)
sig_values <- seq(0.1,2,0.4)
A <- matrix(0,nrow = 3,ncol = length(sig_values))
rownames(A) <- c('anti_variance','simple_variance','reduction percent')
colnames(A) <- sig_values
count <- 1
for(sig in sig_values){
  x <- Generate_Rayleign(sig,u)
  var_anti <- var(x)
  A[1,count] <- var_anti
  count <- count + 1
}

#simple variables
u <- runif(N)
count <- 1
for(sig in sig_values){
  x <- Generate_Rayleign(sig,u)
  var_simple <- var(x)
  A[2,count] <- var_simple
  count <- count + 1
}
A[3,] <- (A[2,]-A[1,])/A[2,]
A

## -----------------------------------------------------------------------------
g <- function(x){
  return(x^2/sqrt(2*pi)*exp(-x^2/2)) * (x>1)
}
f3 <- function(x) x^2

## -----------------------------------------------------------------------------
#Estimate the integral via importance sampling
set.seed(123)
N <- 10000
replica <- 1000
estimate_f1 <- numeric(replica)
estimate_f2 <- numeric(replica)
estimate_f3 <- numeric(replica)
for(i in 1:replica){
  sample_1 <- rnorm(N,0,1)
  sample_2 <- rnorm(N,3,0.5)
  u <- runif(N)
  sample_3 <- (3*u)^{1/3}
  estimate_f1[i] <- mean(g(sample_1)/dnorm(sample_1)*(sample_1>1))
  estimate_f2[i] <- mean(g(sample_2)/dnorm(sample_2,3,0.5)*(sample_2>1))
  estimate_f3[i] <- mean(g(sample_3)/f3(sample_3)*(sample_3>1))
}
cat('the mean and variance of f1:',mean(estimate_f1),",",var(estimate_f1),'\n')
cat('the mean and variance of f2:',mean(estimate_f2),",",var(estimate_f2),'\n')
cat('the mean and variance of f3:',mean(estimate_f3),",",var(estimate_f3),'\n')

## -----------------------------------------------------------------------------
integrand <- function(x) (x^2 / sqrt(2 * pi)) * exp(-x^2 / 2)

true_value <- integrate(integrand, 1, Inf)$value
true_value


## -----------------------------------------------------------------------------
library(e1071)
library(MASS)

## -----------------------------------------------------------------------------
#Generate data
Generate_skewness <- function(sample_size,mean,se){
  x <- rnorm(sample_size,mean,se)
  skewness_values <- skewness(x)
  return(skewness_values)
}

## -----------------------------------------------------------------------------
#statistical inference 
Skewness_quantiles <- function(skewness_values,quantile_value,mean,se){
  quantiles <- quantile(skewness_values, probs = quantile_value)
  variances <- quantile_value * (1-quantile_value)/(length(skewness_values)*dnorm(quantiles,mean,se))
  SE <- sqrt(variances)
  return(c(quantiles,SE))
}

## -----------------------------------------------------------------------------
#result reporting
Skewness_Result <- function(quantile.hat,quantile.se,quantile_values,sample_size){
  #Generate large sample approximation of quantile
  approx_quantiles <- qnorm(quantile_values, mean = 0, sd = sqrt(6/sample_size))
  table <- matrix(c(quantile.hat, quantile.se, approx_quantiles), nrow = 3, byrow = TRUE)

#rename the matrix
rownames(table) <- c("quantile(estimate)", "quantile(sd)", "quantile(theory)")
colnames(table) <- quantile_values

# output the result
print(table)
}

## -----------------------------------------------------------------------------
set.seed(123)
m <- 1e4
n <- 1e3
mean <- 0
se <- 1
quantile_values <- c(0.025,0.05,0.95,0.975)
Skw <- numeric(m)
quantile.hat <- quantile.se <- numeric(length(quantile_values))
for(i in 1:m){
  Skw[i] <- Generate_skewness(sample_size=n,mean=mean,se=se)
}

for(j in 1:length(quantile_values)){
  quantile.hat[j] <- Skewness_quantiles(Skw,quantile_values[j],mean,se)[1]
  quantile.se[j] <- Skewness_quantiles(Skw,quantile_values[j],mean,se)[2] 
}
Skewness_Result(quantile.hat,quantile.se,quantile_values,n)

## -----------------------------------------------------------------------------
set.seed(123)
m <- 1e4
n <- 10
mean <- 0
se <- 1
quantile_values <- c(0.025,0.05,0.95,0.975)
Skw <- numeric(m)
quantile.hat <- quantile.se <- numeric(length(quantile_values))
for(i in 1:m){
  Skw[i] <- Generate_skewness(sample_size=n,mean=mean,se=se)
}

for(j in 1:length(quantile_values)){
  quantile.hat[j] <- Skewness_quantiles(Skw,quantile_values[j],mean,se)[1]
  quantile.se[j] <- Skewness_quantiles(Skw,quantile_values[j],mean,se)[2] 
}
Skewness_Result(quantile.hat,quantile.se,quantile_values,n)

## -----------------------------------------------------------------------------
#Generate data
Generate_biogauss <- function(sample_size,mu,sigma){
  data_normal <- mvrnorm(sample_size, mu, sigma)
  x_normal <- data_normal[, 1]
  y_normal <- data_normal[, 2]
  return(list(x=x_normal,y=y_normal))
}

Generate_alternative <- function(sample_size){
  x_dep <- rnorm(sample_size)
  y_dep <- sqrt(abs(x_dep)) + rnorm(sample_size, sd = 0.1)  # create dependent relationship
  return(list(x=x_dep,y=y_dep))
  }

## -----------------------------------------------------------------------------
#Statistical Inference
Compare_Test <- function(x,y){
  pearson_result <- cor.test(x, y)
  spearman_result <- cor.test(x, y, method = "spearman")
  kendall_result <- cor.test(x, y, method = "kendall")
  return(list(p=pearson_result,s=spearman_result,k=kendall_result))
}

## -----------------------------------------------------------------------------
#result reporting
Correlation_Result <-function(normal,alternative){
  normal_value <- c(normal$p$p.value,normal$s$p.value,normal$k$p.value)
  alternative_value <- c(alternative$p$p.value,alternative$s$p.value,alternative$k$p.value)
table <- matrix(c(normal_value,alternative_value),byrow = T,nrow = 2)
#rename the matrix
rownames(table) <- c("normal condition", "other condition")
colnames(table) <- c('pearson','spearman','kendall')

# output the result
print(table)
}

## -----------------------------------------------------------------------------
#Experiment procedure
## parameter settings
set.seed(123)
n <- 10
mu <- c(0,0)
sigma <- matrix(c(1,0.8,0.8,1),2,2)
rnormal <- Generate_biogauss(n,mu,sigma)
ralternative <- Generate_alternative(n)
normal <- Compare_Test(rnormal$x,rnormal$y)
alternative <- Compare_Test(ralternative$x,ralternative$y)
Correlation_Result(normal,alternative)

## -----------------------------------------------------------------------------
#data generation
data_generation1 <- function(N0,N1){
  p_null <- runif(N0)                            # Null hypotheses
  p_alt <- rbeta(N1, 0.1, 1)                     # Alternative hypotheses
  p_values <- c(p_null, p_alt)
  return(p_values)
}

## -----------------------------------------------------------------------------
# Statistical Inference
statistical_inference1 <- function(p_values,alpha,N0,N1){
  # Bonferroni correction
  bonf_rejected <- (p.adjust(p_values, method = "bonferroni") < alpha)
  # Benjamini-Hochberg correction
  bh_rejected <- (p.adjust(p_values, method = "BH") < alpha)
  
  # Functions to calculate FWER, FDR, TPR
calc_metrics <- function(rejected, true_alt) {
  FWER <- as.numeric(sum(rejected[1:N0])>0)   # FWER: At least one false rejection
  FDR <- sum(rejected[1:N0]) / max(1, sum(rejected))
  TPR <- sum(rejected[(N0 + 1):(N0+N1)]) / N1
  return(c(FWER = FWER, FDR = FDR, TPR = TPR))
}
  
  # Update metrics for Bonferroni correction
  bonf_metrics <- calc_metrics(bonf_rejected, true_alt = p_values[(N0+1):(N0+N1)])
  
  
  # Update metrics for B-H correction
  bh_metrics <- calc_metrics(bh_rejected, true_alt = p_values[(N0+1):(N0+N1)])
  return(list(bonf_metrics=bonf_metrics,bh_metrics=bh_metrics))
}

## -----------------------------------------------------------------------------
#Result Reporting
result_reporting1 <- function(results,bonf_metrics,bh_metrics,m){
  
# Store averages of metrics
  results[,1 ] <- results[,1 ] + bonf_metrics/m 
  results[,2 ] <- results[,2] + bh_metrics /m
  return(results)
}

## -----------------------------------------------------------------------------
set.seed(123)

# Simulation parameters
N <- 1000      # Total number of hypotheses
N0 <- 950        # Number of null hypotheses
N1 <- N - N0     # Number of alternative hypotheses
m <- 10000       # Number of simulation replicates
alpha <- 0.1     # Nominal significance level
# Storing results
  results <- matrix(0, nrow = 3, ncol = 2)
  colnames(results) <- c("Bonferroni correction", "B-H correction")
  rownames(results) <- c("FWER", "FDR", "TPR")

for (i in 1:m) {
  p_values <- data_generation1(N0,N1)
  
  bonf_metrics <- statistical_inference1(p_values,alpha,N0,N1)$bonf_metrics
  bh_metrics <- statistical_inference1(p_values,alpha,N0,N1)$bh_metrics
  
  results <- result_reporting1(results,bonf_metrics,bh_metrics,m)
  
}

print(results)


## -----------------------------------------------------------------------------
#loading package
library(boot)

## -----------------------------------------------------------------------------
#data generation
data_generation2 <- function(){
  # Air-conditioning failure times data
  failure_times <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
  return(failure_times)
}

## -----------------------------------------------------------------------------
#statistical inference
statistical_inference2 <- function(data,Repetitions){
  # Function to compute MLE of lambda for bootstrap samples
  lambda_mle_fn <- function(data, indices) {
  resampled_data <- data[indices]
  return(1 / mean(resampled_data))
}

  # Perform bootstrap
  bootstrap_results <- boot(data = data, statistic = lambda_mle_fn, R = Repetitions)
  return(bootstrap_results)
}

## -----------------------------------------------------------------------------
#result reporting
result_reporting2 <- function(lambda_mle,bootstrap_results){
  # Extract bias and standard error
  lambda_mle_bias <- mean(bootstrap_results$t) - lambda_mle
  lambda_mle_se <- sd(bootstrap_results$t)
  # Output the results
  cat("MLE of lambda:", lambda_mle, "\n")
  cat("Bootstrap estimate of bias:", lambda_mle_bias, "\n")
  cat("Bootstrap estimate of standard error:", lambda_mle_se, "\n")
  cat('Relative Error of Estimate:',abs(lambda_mle_bias)/lambda_mle*100,'%')
}

## -----------------------------------------------------------------------------
#Experimental Settings
set.seed(123)
Repetitions <- 1000 
data <- data_generation2()
# MLE of lambda
lambda_mle <- 1 / mean(data)

bootstrap_results <- statistical_inference2(data,Repetitions)
result_reporting2(lambda_mle,bootstrap_results)

## -----------------------------------------------------------------------------
#statistical inference
statistical_inference3 <- function(data,Repetitions,alpha){
  # Function to compute MLE of 1/lambda for bootstrap samples
  lambda_mle_fn2 <- function(data, indices) {
  resampled_data <- data[indices]
  return(mean(resampled_data))
}

  # Perform bootstrap
  bootstrap_results <- boot(data = data, statistic = lambda_mle_fn2, R = Repetitions)
  ci <- boot.ci(bootstrap_results,type=c("norm","basic","perc","bca"),conf = alpha)
  return(ci)
}

## -----------------------------------------------------------------------------
#result reporting
result_reporting3 <- function(lambdad_mle,norm_ci,basic_ci,perc_ci,BCa_ci){
  # Output results
  mu <- 1 / lambda_mle
  cat("Mean time between failures:", mu, "\n")
cat("Normal CI:(",norm_ci[1],',',norm_ci[2],')', "\n")
cat("Basic CI:(",basic_ci[1],',',basic_ci[2],')', "\n")
cat("Percentile CI:(",perc_ci[1],',',perc_ci[2],')', "\n")
cat("BCa CI:(",BCa_ci[1],',',BCa_ci[2],')', "\n")

}

## -----------------------------------------------------------------------------
#Experimental Settings
alpha <- 0.95
m <- 1000
Repetitions <- 1000

ci.norm<-ci.basic<-ci.perc<-ci.bca<-matrix(NA,m,2)
for(i in 1:m){
ci <- statistical_inference3(data,Repetitions,alpha)
ci.norm[i,]<-ci$norm[2:3];ci.basic[i,]<-ci$basic[4:5]
ci.perc[i,]<-ci$percent[4:5];ci.bca[i,]<-ci$bca[4:5]
}

norm_ci <- c(mean(ci.norm[,1]),mean(ci.norm[,2]))
basic_ci <- c(mean(ci.basic[,1]),mean(ci.basic[,2]))
perc_ci <- c(mean(ci.perc[,1]),mean(ci.perc[,2]))
BCa_ci <- c(mean(ci.bca[,1]),mean(ci.bca[,2]))

result_reporting3(lambda_mle,norm_ci,basic_ci,perc_ci,BCa_ci)


## -----------------------------------------------------------------------------
library(bootstrap)
library(DAAG)

## -----------------------------------------------------------------------------
#generating data
scor_generate <- function(){
  data(scor)
  return(scor)
}

## -----------------------------------------------------------------------------
#statistical inference
jackknife_infer <- function(data){
  # Jackknife 方法
  n <- nrow(data)
  theta_jackknife <- numeric(n)

  for (i in 1:n) {
    # 排除第 i 个观测
    data_excluded <- data[-i, ]
    Sigma_hat_i <- cov(data_excluded)
    eigenvalues_i <- eigen(Sigma_hat_i)$values
    lambda_hat_i <- sort(eigenvalues_i, decreasing = TRUE)
  
    # 计算每次排除的theta
    theta_jackknife[i] <- lambda_hat_i[1] / sum(lambda_hat_i)
  }

  return(theta_jackknife)
}

## -----------------------------------------------------------------------------
#result reporting
jackknife_report <- function(data,theta_jackknife){
  n <- nrow(data)
  # 计算协方差矩阵的MLE
  Sigma_hat <- cov(data)

  # 计算协方差矩阵的特征值
  eigenvalues <- eigen(Sigma_hat)$values
  lambda_hat <- sort(eigenvalues, decreasing = TRUE)

  # 计算样本估计的theta
  hat_theta <- lambda_hat[1] / sum(lambda_hat)
  # 计算 Jackknife 偏差和标准误
  bias_theta <- (n - 1) * (mean(theta_jackknife) - hat_theta)
  se_theta <- sqrt((n - 1) / n * sum((theta_jackknife - mean(theta_jackknife))^2))
  cat('the jackknife estimates of theta is',(n-1) * mean(theta_jackknife),'\n')
  cat('the jackknife estimates of bias is',bias_theta,'\n')
  cat('the jackknife estimates of standard error is',se_theta,'\n')
}


## -----------------------------------------------------------------------------
data <- scor_generate()
theta_jackknife <- jackknife_infer(data)
jackknife_report(data,theta_jackknife)


## -----------------------------------------------------------------------------
#data generating
ironslag_generate <- function(){
  attach(ironslag)
 
  data <- list(chemical=chemical,magnetic=magnetic)
  detach(ironslag)
  return(data)
}

## -----------------------------------------------------------------------------
#Statistical inference
Cross_Validation <- function(chemical,magnetic){
  n <- length(magnetic)
  e1 <- e2 <- e3 <- e4 <- numeric(n)
  
  #for n-fold cross validation
  #fit models on leave-one-out samples
  for(k in 1:n){
    x <- chemical[-k]
    y <- magnetic[-k]
    
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
    e1[k] <- magnetic[k] - yhat1
    
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k]+J2$coef[3] * chemical[k]^2
    e2[k] <- magnetic[k] - yhat2
    
    J3 <- lm(log(y) ~ x )
    logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
    yhat3 <- exp(logyhat3)
    e3[k] <- magnetic[k] - yhat3
    
    J4 <- lm(y ~ poly(x, 3, raw = TRUE))
    yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k] + J4$coef[3] * chemical[k]^2+ J4$coef[4] * chemical[k]^3
    e4[k] <- magnetic[k] - yhat4
  }
  mse <- c(mean(e1^2),mean(e2^2),mean(e3^2),mean(e4^2))
  x <- chemical
  y <- magnetic
  J1 <- lm(y ~ x)
  J2 <- lm(y ~ x + I(x^2))
  J3 <- lm(log(y) ~ x )
  J4 <- lm(y ~ poly(x, 3, raw = TRUE))
  model_result <- list(mse = mse,model1=J1,model2=J2,model3=J3,model4=J4) 
  return(model_result)
}

## -----------------------------------------------------------------------------
#Result Reporting
CV_result <- function(model_result){
  adj_r <- c(summary(model_result$model1)$adj.r.squared,summary(model_result$model2)$adj.r.squared,summary(model_result$model3)$adj.r.squared,summary(model_result$model4)$adj.r.squared)
 
  pred_error <- matrix(c(model_result$mse,adj_r),nrow = 2,ncol = 4,byrow = T)
  colnames(pred_error) <- c('Linear','Quadratic','Exponential','Cubic')
  rownames(pred_error) <- c('prediction error','adj_Rsquare')
  print(pred_error)
}

## -----------------------------------------------------------------------------
data <- ironslag_generate()
chemical <- as.vector(data$chemical)
magnetic <- as.vector(data$magnetic)
model_result <- Cross_Validation(chemical,magnetic)
CV_result(model_result)


## -----------------------------------------------------------------------------
#data generate
chickwts_generate <- function(){
  attach(chickwts)
  x <- sort(as.vector(weight[feed == 'soybean']))
  y <- sort(as.vector(weight[feed == 'linseed']))
  detach(chickwts)
  data <- list(x=x,y=y)
  return(data)
}

## -----------------------------------------------------------------------------
#statistical inference
# Cramér-von Mises test function
cramer_von_mises <- function(x, y) {
  n <- length(x)
  m <- length(y)
  combined <- c(x, y)
  
  # Calculate empirical CDFs
  ecdf_x <- ecdf(x)
  ecdf_y <- ecdf(y)
  
  # Compute the Cramér-von Mises statistic
  W <- sum((ecdf_x(combined) - ecdf_y(combined))^2) * (1 / (n * m))
  
  return(W)
}

# Permutation test function
permutation_test <- function(x, y, n_permutations) {
  observed_stat <- cramer_von_mises(x, y)
  combined <- c(x, y)
  count <- 0
  
  for (i in 1:n_permutations) {
    permuted <- sample(combined)
    new_x <- permuted[1:length(x)]
    new_y <- permuted[(length(x) + 1):length(permuted)]
    permuted_stat <- cramer_von_mises(new_x, new_y)
    
    if (permuted_stat >= observed_stat) {
      count <- count + 1
    }
  }
  
  p_value <- count / n_permutations
  return(p_value)
}



## -----------------------------------------------------------------------------
#Result Reporting
permutation_result <- function(p_value){
  cat("P-value:", p_value, "\n")
}

## -----------------------------------------------------------------------------
set.seed(123)
R <- 1000
x <- chickwts_generate()$x
y <- chickwts_generate()$y
# 执行置换测试
p_value <- permutation_test(x, y, R)
permutation_result(p_value)

## -----------------------------------------------------------------------------
#data generating
random_generate <- function(n){
  x <- runif(n)
  y <- sin(x)
  return(list(x=x,y=y))
}

## -----------------------------------------------------------------------------
##Statistical inference
# Spearman rank correlation permutation test function
permutation_spearman_test <- function(x, y, n_permutations) {
  observed_cor <- cor(x, y, method = "spearman")
  count <- 0
  
  combined <- data.frame(x, y)
  
  for (i in 1:n_permutations) {
    # Permute the y-values
    permuted_y <- sample(y)
    
    # Calculate Spearman correlation on permuted data
    permuted_cor <- cor(x, permuted_y, method = "spearman")
    
    if (abs(permuted_cor) >= abs(observed_cor)) {
      count <- count + 1
    }
  }
  
  p_value <- count / n_permutations
  return(p_value)
}




## -----------------------------------------------------------------------------
#Result Reporting
Spearman_reporting <- function(x,y,permutation_p_value){
  # 计算 Spearman 相关性和 p 值
spearman_test_result <- cor.test(x, y, method = "spearman")
# 输出结果
cat("Spearman correlation coefficient:", spearman_test_result$estimate, "\n")
cat("P-value from cor.test:", spearman_test_result$p.value, "\n")
cat("P-value from permutation test:", permutation_p_value, "\n")
}

## -----------------------------------------------------------------------------
set.seed(123)
R <- 1000
n <- 500
x <- random_generate(n)$x
y <- random_generate(n)$y

# 执行置换测试
permutation_p_value <- permutation_spearman_test(x, y,R)
Spearman_reporting(x,y,permutation_p_value)

## -----------------------------------------------------------------------------
#upload package
library(ggplot2)

## -----------------------------------------------------------------------------
## data generating
# Metropolis-Hastings sampler function
metropolis_hastings <- function(n, proposal_sd,X0) {
  # Initialize the chain
  x <- numeric(n)
  x[1] <- X0  # Start from some given value
  
  for (i in 2:n) {
    # Current value
    current_x <- x[i - 1]
    
    # Propose a new value from a normal distribution
    proposed_x <- rnorm(1, mean = current_x, sd = proposal_sd)
    
    # Calculate acceptance probability
    acceptance_prob <- dcauchy(proposed_x) / dcauchy(current_x)
    
    # Accept or reject the proposed value
    if (runif(1) < acceptance_prob) {
      x[i] <- proposed_x
    } else {
      x[i] <- current_x
    }
  }
  
  return(x)
}

## -----------------------------------------------------------------------------
## statistical inference
deciles_calculating <- function(samples){
# Calculate deciles of the generated samples
sample_deciles <- quantile(samples, probs = seq(0, 1, 0.1))

# Calculate deciles of the standard Cauchy distribution
cauchy_deciles <- qcauchy(seq(0, 1, 0.1))

# Create a comparison data frame
comparison <- data.frame(
  Decile = seq(0, 1, 0.1),
  Sample_Deciles = sample_deciles,
  Cauchy_Deciles = cauchy_deciles
)

# Print the comparison data frame
return(comparison)
}



## -----------------------------------------------------------------------------
## result reporting
Comparison_plot <- function(comparison){
# Plot the comparison
ggplot(comparison, aes(x = Decile)) +
  geom_line(aes(y = Sample_Deciles, color = "Sample Deciles"), linewidth = 1) +
  geom_line(aes(y = Cauchy_Deciles, color = "Cauchy Deciles"), linewidth = 1) +
  labs(title = "Deciles Comparison", y = "Value", x = "Decile") +
  scale_color_manual(values = c("Sample Deciles" = "blue", "Cauchy Deciles" = "red")) +
  theme_minimal()
}

## -----------------------------------------------------------------------------
set.seed(123)
n_samples <- 10000
discard_samples <- 1000
X0 = -5
proposal_sd <- 1  # Standard deviation of the proposal distribution
samples <- metropolis_hastings(n_samples, proposal_sd,X0)
# Discard the first 1000 samples
samples <- samples[-(1:discard_samples)]

comparison <- deciles_calculating(samples)

Comparison_plot(comparison)


## -----------------------------------------------------------------------------
#data generating
# Gibbs sampler function
gibbs_sampler <- function(n_samples, n, a, b) {
  # Initialize vectors to store samples
  x_samples <- numeric(n_samples)
  y_samples <- numeric(n_samples)

  # Initial values
  x_samples[1] <- sample(0:n,1)
  y_samples[1] <- runif(1)  # Starting with a reasonable value for y

  for (i in 2:n_samples) {
    # Sample x given y
    y_current <- y_samples[i - 1]
    x_samples[i] <- rbinom(1, n, y_current)

    # Sample y given x
    x_current <- x_samples[i]
    alpha <- x_current + a
    beta <- n - x_current + b
    y_samples[i] <- rbeta(1, alpha, beta)
  }

  return(data.frame(X = x_samples, Y = y_samples))
}

## -----------------------------------------------------------------------------
## statistical inference

#I think no statistical inference need for this problem

## -----------------------------------------------------------------------------
## result reporting
gibbs_plot <- function(samples){
# Plot the results
ggplot(samples, aes(x = X, y = Y)) +
  geom_point(alpha = 0.2) +
  labs(title = "Gibbs Sampler: Bivariate Density", x = "X", y = "Y") +
  theme_minimal()
}

## -----------------------------------------------------------------------------
set.seed(123)
# Parameters
n <- 10  # Total trials
a <- 2   # Shape parameter for Beta distribution
b <- 3   # Shape parameter for Beta distribution
n_samples <- 10000  # Number of samples
discard_samples <- 1000


# Generate samples
samples <- gibbs_sampler(n_samples, n, a, b)
samples <- samples[-(1:discard_samples),]

gibbs_plot(samples)

## -----------------------------------------------------------------------------
#data generating

#Metropolis-Hasting function with multiple chains
metropolis_hastings_multiple_chains <- function(n_chains,n_samples,proposal_sd,x0){
  # Initialize a matrix to store samples for each chain
  chains <- matrix(0, nrow = n_samples, ncol = n_chains)
  # Run multiple chains
  
  for(chain in 1:n_chains){
    # Store samples in the chains matrix
    chains[,chain] <- metropolis_hastings(n_samples,proposal_sd,x0[chain])
  }
  return(chains)
}




# Gibbs sampler function with multiple chains
gibbs_sampler_multiple_chains <- function(n_chains, n_samples, n, a, b) {
  # Initialize a matrix to store samples for each chain
  chains_x <- matrix(0, nrow = n_samples, ncol = n_chains)
  chains_y <- matrix(0, nrow = n_samples, ncol = n_chains)

  # Run multiple chains
  
  for(chain in 1:n_chains){
    # Store samples in the chains matrix
    chains_x[, chain] <- gibbs_sampler(n_samples, n, a, b)$X
    chains_y[, chain] <- gibbs_sampler(n_samples, n, a, b)$Y
  }
  return(list(X = chains_x, Y = chains_y))
}


## -----------------------------------------------------------------------------
#statistical inference
# Function to compute Gelman-Rubin diagnostic
Gelman.Rubin <- function(psi) {
        # psi[i,j] is the statistic psi(X[i,1:j])
        # for chain in i-th row of X
        psi <- as.matrix(psi)
        n <- ncol(psi)
        k <- nrow(psi)

        psi.means <- rowMeans(psi)     #row means
        B <- n * var(psi.means)        #between variance est.
        psi.w <- apply(psi, 1, "var")  #within variances
        W <- mean(psi.w)               #within est.
        v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
        r.hat <- v.hat / W             #G-R statistic
        return(r.hat)
}


## -----------------------------------------------------------------------------
#Result Reporting
GB_monitor_plot <- function(chains,n_chains,b,n_samples){
  #compute diagnostic statistics
    psi <- t(apply(chains, 1, cumsum))
    for (i in 1:nrow(psi))
        psi[i,] <- psi[i,] / (1:ncol(psi))

    #plot psi for the four chains
#    par(mfrow=c(2,2))
    for (i in 1:n_chains)
      if(i==1){
        plot((b+1):n_samples,psi[i, (b+1):n_samples], type="l",
            xlab='Index', ylab=bquote(phi))
      }else{
        lines(psi[i, (b+1):n_samples], col=i)
    }
    par(mfrow=c(1,1)) #restore default
    #plot the sequence of R-hat statistics
    rhat <- rep(0, n_samples)
    for (j in 1:n_samples){
        rhat[j] <- Gelman.Rubin(psi[,1:j])}
    plot((b+1):n_samples,rhat[(b+1):n_samples], type="l", xlab="", ylab="R")
    abline(h=1.1, lty=2)
}

## -----------------------------------------------------------------------------
set.seed(123)
# Parameters
n <- 10  # Total trials
a <- 2   # Shape parameter for Beta distribution
b <- 3   # Shape parameter for Beta distribution
n_chains <- 4  # Number of chains
n_samples <- 15000  # Number of samples
proposal_sd <- 1
b <- 1000       #burn-in length
#choose overdispersed initial values
x0 <- c(-5,-10, 5, 10)

#Method1
# Run the Metropolis-Hastings algorithm until convergence
convergence_threshold <- 1.2
R_hat <- Inf
samples <- NULL

while (R_hat >= convergence_threshold) {
  chains1 <- metropolis_hastings_multiple_chains(n_chains,n_samples,proposal_sd,x0)
  chains1 <- t(chains1)
  R_hat<- Gelman.Rubin(chains1)
  
}
print(R_hat)
GB_monitor_plot(chains1,n_chains,b,n_samples)



## -----------------------------------------------------------------------------
set.seed(1234)
# Run the Gibbs sampler until convergence
convergence_threshold <- 1.2
R_hat <- Inf
samples <- NULL

while (R_hat >= convergence_threshold) {
  chains2 <- gibbs_sampler_multiple_chains(n_chains, n_samples, n, a, b)
  X <- t(chains2$X)
  Y <- t(chains2$Y)
  # Compute R_hat for both X and Y
  R_hat_x <- Gelman.Rubin(X)
  R_hat_y <- Gelman.Rubin(Y)
  R_hat <- max(R_hat_x, R_hat_y)
  
}

GB_monitor_plot(X,n_chains,b,n_samples)
GB_monitor_plot(Y,n_chains,b,n_samples)



## -----------------------------------------------------------------------------
kth_term <- function(k, a, d) {
  # compute Euclid norm ||a||
  norm_a <- sqrt(sum(a^2))
  
  # 计算每一部分
  part1 <- (-1)^k / (factorial(k) * 2^k)
  part2 <- norm_a^(2 * k + 2) / ((2 * k + 1) * (2 * k + 2))
  part3 <- gamma((d + 1) / 2) * gamma(k + 1.5) / gamma(k + d / 2 + 1)
  
  return(part1 * part2 * part3)
}

## -----------------------------------------------------------------------------
series_sum <- function(a, d, tol = 1e-10, max_terms = 1000) {
  total_sum <- 0
  k <- 0
  
  repeat {
    term <- kth_term(k, a, d)
    total_sum <- total_sum + term
    
    # if the value is less than threshold,exit!!!
    if (abs(term) < tol || k >= max_terms) {
      break
    }
    
    k <- k + 1
  }
  
  return(total_sum)
}

## -----------------------------------------------------------------------------
a <- c(1, 2)
d <- 2
result <- series_sum(a, d)
result

## -----------------------------------------------------------------------------
# Define S_k-1 and S_k
S_k_minus_1 <- function(a, k) {
  1 - pt(sqrt(a^2 * (k - 1) / (k - a^2)), df = k - 1)
}

S_k <- function(a, k) {
  1 - pt(sqrt(a^2 * k / (k + 1 - a^2)), df = k)
}

# Define function to find root using Bisection method
find_intersection_bisection <- function(k, tol = 1e-6, max_iter = 1000) {
  # Define the function to find zero
  f <- function(a) S_k_minus_1(a, k) - S_k(a, k)
  
   # Initial interval (0, sqrt(k))
  epsilon <- 1e-3
  lower <- epsilon
  upper <- sqrt(k)-epsilon
  
  if (f(lower) * f(upper) > 0) {
    stop("No sign change detected in the initial interval; adjust the interval.")
  }
  
  
  # Bisection method loop
  for (i in 1:max_iter) {
    mid <- (lower + upper) / 2
    f_mid <- f(mid)
    
    if (abs(f_mid) < tol) {
      return(mid)
    } else if (f(lower) * f_mid < 0) {
      upper <- mid
    } else {
      lower <- mid
    }
  }
  
  warning("Bisection method did not converge.")
  return(NA)
}

# Calculate intersection points for specified values of k
k_values <- c(4:25, 100,500,1000)

root <- numeric(length(k_values))
for(i in 1:length(k_values)){
  root[i] <- find_intersection_bisection(k_values[i])
}

names(root) <- paste0("k=", k_values)
root




## ----error=TRUE---------------------------------------------------------------
# Define c_k and c_{k-1} as functions of a and k
c_k <- function(a, k) {
  sqrt(a^2 * k / (k + 1 - a^2))
}

c_k_minus_1 <- function(a, k) {
  sqrt(a^2 * (k - 1) / (k - a^2))
}

# Define the left side integral
integral_left <- function(a, k) {
  c_km1 <- c_k_minus_1(a, k)
  integrand <- function(u) (1 + u^2 / (k - 1))^(-k / 2)
  return(integrate(integrand, 0, c_km1)$value)
}

# Define the right side integral
integral_right <- function(a, k) {
  c_k_val <- c_k(a, k)
  integrand <- function(u) (1 + u^2 / k)^(-(k + 1) / 2)
  return(integrate(integrand, 0, c_k_val)$value)
}

# Define the full equation difference function
equation_difference <- function(a, k) {
  gamma_ratio_left <- 2 * gamma(k / 2) / (sqrt(pi * (k - 1)) * gamma((k - 1) / 2))
  gamma_ratio_right <- 2 * gamma((k + 1) / 2) / (sqrt(pi * k) * gamma(k / 2))
  
  left_side <- gamma_ratio_left * integral_left(a, k)
  right_side <- gamma_ratio_right * integral_right(a, k)
  
  return(left_side - right_side)
}


# Define function to find root using Bisection method
find_a <- function(k, tol = 1e-6, max_iter = 1000) {
  # Define the function to find zero
  f <- function(a) equation_difference(a,k)
  
  # Initial interval (0, sqrt(k))
  epsilon <- 1e-2
  lower <- epsilon * k
  upper <- sqrt(k)-epsilon * k
  
  if (f(lower) * f(upper) > 0) {
    stop("No sign change detected in the initial interval; adjust the interval.")
  }
  # Bisection method loop
  for (i in 1:max_iter) {
    mid <- (lower + upper) / 2
    f_mid <- f(mid)
    
    if (abs(f_mid) < tol) {
      return(mid)
    } else if (f(lower) * f_mid < 0) {
      upper <- mid
    } else {
      lower <- mid
    }
  }
  
  warning("Bisection method did not converge.")
  return(NA)
}



# Compare solutions for specified k values
k_values <- c(4:25, 100, 500, 1000)
solutions <- numeric(length(k_values))
for(i in 1:length(k_values)){
  solutions[i] <- find_a(k_values[i])
  print(c(k_values[i],solutions[i]))
  
}

names(solutions) <- paste0("k=", k_values)
solutions

## -----------------------------------------------------------------------------
compare <- matrix(c(root,solutions),nrow = 2,byrow = T)
colnames(compare) <- k_values
rownames(compare) <- c('Exercise 11.4','Exercise 11.5')
compare

## -----------------------------------------------------------------------------
find_root_of_equation <- function(k) {

  # General integral function
  integral_expr <- function(u, n) {
    (1 + u^2 / (n - 1))^(-n / 2)
  }
  
  # Calculate constant c_k
  calculate_c <- function(n, a) {
    sqrt(a^2 * n / (n + 1 - a^2))
  }
  
  # Compute the left or right side of the equation
  compute_expression <- function(n, a) {
    c <- calculate_c(n - 1, a)
    integral_value <- integrate(function(u) integral_expr(u, n), lower = 0, upper = c)$value
    2 / sqrt(pi * (n - 1)) * exp(lgamma(n / 2) - lgamma((n - 1) / 2)) * integral_value
  }
  
  # Define the function for the difference between left and right side
  difference_function <- function(a) {
    left_side <- compute_expression(k, a)
    right_side <- compute_expression(k + 1, a)
    left_side - right_side
  }
  
  # Define a small epsilon for interval boundaries
  epsilon <- 1e-2
  # Check if there is a root within the interval
  if ((difference_function(epsilon) < 0 && difference_function(sqrt(k) - epsilon) > 0) ||
      (difference_function(epsilon) > 0 && difference_function(sqrt(k) - epsilon) < 0)) {
    root <- uniroot(difference_function, interval = c(epsilon, sqrt(k) - epsilon))$root
  } else {
    root <- NA
  }
  
  return(root)
}

# Apply the function to a range of k values
solutions <- sapply(c(4:25, 100, 500, 1000), function(k) {
  find_root_of_equation(k)
})
solutions

## -----------------------------------------------------------------------------
compare <- matrix(c(root,solutions),nrow = 2,byrow = T)
colnames(compare) <- k_values
rownames(compare) <- c('Exercise 11.4','Exercise 11.5')
compare

## -----------------------------------------------------------------------------
# Observed data
Y <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
tau <- 1
n <- length(Y)

# Initial estimate of lambda
lambda <- mean(Y)
tolerance <- 1e-6
diff <- 1
max_iter <- 1000
iter <- 0

# E-M Algorithm
while (diff > tolerance && iter < max_iter) {
  lambda_old <- lambda
  # Separate uncensored and censored observations
  uncensored <- Y[Y < tau]
  censored_count <- sum(Y == tau)
  
  # E-Step: Expected total for censored values
  expected_censored_sum <- censored_count * (tau + lambda_old)
  
  # M-Step: Update lambda
  lambda <- (sum(uncensored) + expected_censored_sum) / n
  
  # Check for convergence
  diff <- abs(lambda - lambda_old)
  iter <- iter + 1
}

cat("Estimated lambda using EM algorithm:", lambda, "\n")

## -----------------------------------------------------------------------------
# Calculate MLE directly for comparison (using uncensored data only)
lambda_mle <- mean(uncensored) / (1 - censored_count / n)

cat("Observed data MLE for lambda:", lambda_mle, "\n")

## ----warning=FALSE------------------------------------------------------------
# 加载 lpSolve 包
library(lpSolve)

## -----------------------------------------------------------------------------

# 定义目标函数的系数 (Minimize 4x + 2y + 9z)
objective <- c(4, 2, 9)

# 定义约束矩阵 (左侧系数)
constraints <- matrix(c(
  2,  1,  1,  # 2x + y + z <= 2
  1, -1,  3   # x - y + 3z <= 3
), nrow = 2, byrow = TRUE)

# 定义约束的右侧值
rhs <- c(2, 3)

# 定义约束类型
constraint_directions <- c("<=", "<=")

# 求解线性规划
solution <- lp("min", objective, constraints, constraint_directions, rhs, 
               all.int = FALSE) # 确保变量是非负连续的

# 输出结果
if (solution$status == 0) {
  cat("Optimal Solution Found:\n")
  cat("Objective Value (Minimum):", solution$objval, "\n")
  cat("Values of x, y, z:", solution$solution, "\n")
} else {
  cat("No Optimal Solution Found.\n")
}


## -----------------------------------------------------------------------------
# 加载数据集
data(mtcars)

# 定义公式列表
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

## -----------------------------------------------------------------------------
# 方法1: 使用 for 循环
lm_results_for <- list()  # 用于存储结果的列表
for (i in seq_along(formulas)) {
  lm_results_for[[i]] <- lm(formulas[[i]], data = mtcars)
}
names(lm_results_for) <- paste0("Model_", seq_along(formulas))  # 给结果命名
print("Results using for loop:")
print(lm_results_for)

## -----------------------------------------------------------------------------
# 方法2: 使用 lapply 函数
lm_results_lapply <- lapply(formulas, function(formula) lm(formula, data = mtcars))
names(lm_results_lapply) <- paste0("Model_", seq_along(formulas))  # 给结果命名
print("Results using lapply:")
print(lm_results_lapply)

## -----------------------------------------------------------------------------
# 引导样本生成
set.seed(123)  # 确保结果可重复
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), replace = TRUE)
  mtcars[rows, ]
})


## -----------------------------------------------------------------------------
# 初始化存储结果的列表
lm_results_for <- list()

# 使用 for 循环拟合模型
for (i in seq_along(bootstraps)) {
  lm_results_for[[i]] <- lm(mpg ~ disp, data = bootstraps[[i]])
}

# 检查结果
print("Results using for loop:")
print(lm_results_for)


## -----------------------------------------------------------------------------
# 定义命名函数用于拟合模型
fit_model <- function(data) {
  lm(mpg ~ disp, data = data)
}

# 使用 lapply 调用命名函数
lm_results_lapply <- lapply(bootstraps, fit_model)

# 检查结果
print("Results using lapply:")
print(lm_results_lapply)


## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared


## -----------------------------------------------------------------------------
# 假设模型列表来自之前的 for 循环 lm_results_for
r2_for <- numeric(length(lm_results_for))  # 初始化存储 R^2 的向量

for (i in seq_along(lm_results_for)) {
  r2_for[i] <- rsq(lm_results_for[[i]])
}

# 输出 R^2 结果
print("R^2 values using for loop:")
print(r2_for)


## -----------------------------------------------------------------------------
# 假设模型列表来自之前的 lapply lm_results_lapply
r2_lapply <- sapply(lm_results_lapply, rsq)  # sapply 将结果简化为向量

# 输出 R^2 结果
print("R^2 values using lapply:")
print(r2_lapply)


## -----------------------------------------------------------------------------
# 假设模型列表来自引导样本的 for 循环 lm_results_for
bootstrap_r2_for <- numeric(length(lm_results_for))

for (i in seq_along(lm_results_for)) {
  bootstrap_r2_for[i] <- rsq(lm_results_for[[i]])
}

# 输出 R^2 结果
print("Bootstrap R^2 values using for loop:")
print(bootstrap_r2_for)


## -----------------------------------------------------------------------------
# 假设模型列表来自引导样本的 lapply lm_results_lapply
bootstrap_r2_lapply <- sapply(lm_results_lapply, rsq)

# 输出 R^2 结果
print("Bootstrap R^2 values using lapply:")
print(bootstrap_r2_lapply)


## -----------------------------------------------------------------------------
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)

## -----------------------------------------------------------------------------
# 使用 sapply 和匿名函数提取 p 值
p_values <- sapply(trials, function(trial) trial$p.value)

# 查看结果
print("Extracted p-values:")
print(p_values)


## -----------------------------------------------------------------------------
# 检查第一个试验的 p 值
cat("P-value from the first trial:", p_values[1], "\n")

# 检查显著性水平（例如 0.05）的比例
significant <- mean(p_values < 0.05)
cat("Proportion of significant p-values (p < 0.05):", significant, "\n")


## -----------------------------------------------------------------------------
# 绘制 p 值的分布
hist(p_values, breaks = 20, main = "Distribution of p-values", xlab = "p-value", col = "skyblue")


## -----------------------------------------------------------------------------
parallel_lapply <- function(..., FUN, FUN.VALUE) {
  # 使用 Map 将多个输入组合成列表
  inputs <- Map(list, ...)
  
  # 使用 vapply 并行应用函数，确保输出为向量或矩阵
  vapply(inputs, function(args) do.call(FUN, args), FUN.VALUE = FUN.VALUE)
}


## -----------------------------------------------------------------------------
# 定义一个加法函数
add <- function(x, y) x + y

# 使用 parallel_lapply 实现并行加法
result <- parallel_lapply(1:5, 6:10, FUN = add, FUN.VALUE = numeric(1))
print(result)  # 输出: [1]  7  9 11 13 15


## -----------------------------------------------------------------------------
# 定义一个组合函数
combine <- function(x, y) c(x, y)

# 使用 parallel_lapply 并行组合
result <- parallel_lapply(1:3, 4:6, FUN = combine, FUN.VALUE = numeric(2))
print(result)


## -----------------------------------------------------------------------------
# 加速版卡方检验函数
fast_chisq_test <- function(observed, expected) {
  # 确保输入是数值向量，且长度一致
  if (!is.numeric(observed) || !is.numeric(expected)) {
    stop("输入必须是数值向量。")
  }
  if (length(observed) != length(expected)) {
    stop("观测值和期望值的长度必须一致。")
  }
  
  # 确保没有缺失值
  if (any(is.na(observed)) || any(is.na(expected))) {
    stop("输入向量中不能包含缺失值。")
  }
  
  # 计算卡方统计量
  chi_square_statistic <- sum((observed - expected)^2 / expected)
  
  return(chi_square_statistic)
}

## -----------------------------------------------------------------------------
# 定义观测值和期望值
observed <- c(10, 20, 30, 40)
expected <- c(12, 18, 35, 35)

# 使用加速版函数计算卡方统计量
result <- fast_chisq_test(observed, expected)
print(result)  

## -----------------------------------------------------------------------------
# 使用 chisq.test() 验证
chisq_result <- chisq.test(x = observed, p = expected / sum(expected), rescale.p = TRUE)
print(chisq_result$statistic)  # 输出的卡方统计量应该与 fast_chisq_test 一致

## -----------------------------------------------------------------------------
# 针对两个整数向量的快速 table 函数
fast_table <- function(x, y) {
  # 确保输入是整数向量并且长度一致
  if (!is.integer(x) || !is.integer(y)) {
    stop("两个输入必须是整数向量。")
  }
  if (length(x) != length(y)) {
    stop("两个向量的长度必须一致。")
  }
  
  # 获取 x 和 y 的唯一值
  x_levels <- unique(x)
  y_levels <- unique(y)
  
  # 初始化列联表
  contingency_table <- matrix(0, nrow = length(x_levels), ncol = length(y_levels))
  rownames(contingency_table) <- x_levels
  colnames(contingency_table) <- y_levels
  
  # 填充列联表
  for (i in seq_along(x)) {
    row_index <- match(x[i], x_levels)
    col_index <- match(y[i], y_levels)
    contingency_table[row_index, col_index] <- contingency_table[row_index, col_index] + 1
  }
  
  return(contingency_table)
}

## -----------------------------------------------------------------------------
# 示例输入
x <- as.integer(c(1, 2, 1, 2, 3, 3, 1, 2))
y <- as.integer(c(1, 1, 2, 2, 3, 1, 2, 3))

# 计算列联表
result_table <- fast_table(x, y)
print(result_table)


## -----------------------------------------------------------------------------
# 基于 fast_table 的快速卡方检验
fast_chisq_test_with_table <- function(x, y) {
  # 计算列联表
  contingency_table <- fast_table(x, y)
  
  # 计算行和、列和及总数
  row_totals <- rowSums(contingency_table)
  col_totals <- colSums(contingency_table)
  grand_total <- sum(contingency_table)
  
  # 计算期望频数
  expected <- outer(row_totals, col_totals, FUN = "*") / grand_total
  
  # 计算卡方统计量
  chi_square_statistic <- sum((contingency_table - expected)^2 / expected)
  
  return(chi_square_statistic)
}

## -----------------------------------------------------------------------------
# 示例输入
x <- as.integer(c(1, 2, 1, 2, 3, 3, 1, 2))
y <- as.integer(c(1, 1, 2, 2, 3, 1, 2, 3))

# 计算卡方统计量
chisq_result <- fast_chisq_test_with_table(x, y)
print(chisq_result)

## -----------------------------------------------------------------------------
library(Rcpp)

## ----warning=FALSE------------------------------------------------------------
dir_cpp <- '../src/'

# Can create source file in Rstudio
sourceCpp(paste0(dir_cpp,"gibbs_sampler.cpp"))


# 设置参数
n <- 10
a <- 2
b <- 3
num_samples <- 5000
burn_in <- 500

# 调用Gibbs采样函数
Rcpp_results <- gibbs_sampler(n, a, b, num_samples, burn_in)

Rcpp_results <- data.frame(X = Rcpp_results[,1],Y = Rcpp_results[,2])
# 检查生成的样本
head(Rcpp_results)

# 可视化样本
plot(Rcpp_results[,1], Rcpp_results[,2], main = "Gibbs Sampler Samples", xlab = "x", ylab = "y", pch = 20)

## -----------------------------------------------------------------------------
# Gibbs sampler function
r_gibbs_sampler <- function(n, a, b ,num_samples, burn_in) {
  # Initialize vectors to store samples
  x_samples <- numeric(num_samples)
  y_samples <- numeric(num_samples)

  # Initial values
  x_samples[1] <- sample(0:n,1)
  y_samples[1] <- runif(1)  # Starting with a reasonable value for y

  for (i in 2:num_samples) {
    # Sample x given y
    y_current <- y_samples[i - 1]
    x_samples[i] <- rbinom(1, n, y_current)

    # Sample y given x
    x_current <- x_samples[i]
    alpha <- x_current + a
    beta <- n - x_current + b
    y_samples[i] <- rbeta(1, alpha, beta)
  }

  return(data.frame(X = x_samples[-(1:burn_in)], Y = y_samples[-(1:burn_in)]))
}

## -----------------------------------------------------------------------------
R_results <- r_gibbs_sampler(n, a, b, num_samples, burn_in)

# 检查生成的样本
head(R_results)

# 可视化样本
plot(R_results[,1], R_results[,2], main = "Gibbs Sampler Samples", xlab = "x", ylab = "y", pch = 20)

## -----------------------------------------------------------------------------
# Q-Q 图比较
par(mfrow = c(1, 2)) # 设置图形布局为 1 行 2 列
qqplot(R_results[, 1], Rcpp_results[, 1], 
       main = "QQ Plot of R vs Rcpp (X)",,xlab = 'R',ylab = 'Rcpp')
abline(0, 1, col = "red")
qqplot(R_results[, 2], Rcpp_results[, 2], 
       main = "QQ Plot of R vs Rcpp (Y)",,xlab = 'R',ylab = 'Rcpp')
abline(0, 1, col = "blue")

## -----------------------------------------------------------------------------
# 性能比较
library(microbenchmark)
ts <- microbenchmark(gibbsR=r_gibbs_sampler(n, a, b, num_samples, burn_in),gibbsRcpp=gibbs_sampler(n, a, b, num_samples, burn_in))
summary(ts)[,c(1,3,5,6)]

