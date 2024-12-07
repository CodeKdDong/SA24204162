#include <Rcpp.h>
#include <random>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix gibbs_sampler(int n, double a, double b, int num_samples, int burn_in) {
  // 初始化样本存储矩阵
  NumericMatrix samples(num_samples, 2);

  // 初始化x和y的值
  int x = n / 2;  // 随意选择一个初始值
  double y = 0.5; // 初始值在(0, 1)范围内

  // 设置随机数种子
  std::random_device rd;
  std::mt19937 gen(rd());

  // 采样主循环
  for (int i = 0; i < num_samples + burn_in; ++i) {
    // Step 1: 从条件分布 x | y ~ Binomial(n, y) 中抽样
    x = R::rbinom(n, y);

    // Step 2: 从条件分布 y | x ~ Beta(x + a, n - x + b) 中抽样
    y = R::rbeta(x + a, n - x + b);

    // 如果超出了burn-in阶段，将样本存储下来
    if (i >= burn_in) {
      samples(i - burn_in, 0) = x;
      samples(i - burn_in, 1) = y;
    }
  }

  return samples;
}
