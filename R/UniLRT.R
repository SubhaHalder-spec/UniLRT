#' Likelihood ratio test for testing equality means against umbrella-ordered alternatives (peak = ceiling(k/2) = h) in one-way ANOVA
#' @export
#' @param sample_data list
#' @param significance_level numeric
#' @return Critical value numeric
#' @return Test statistic value numeric
#' @return Result Character
#' @details Testing of H_0:mu_1 = mu_2 = ... = mu_k vs H_1:mu_1 <=.....<= mu_(h-1)<= mu_h >= mu_(h+1)>=....>= mu_k (at least one strict inequality), where mu_i represents the population means of the i-th treatment. The input consists of two variables: sample_data and significance_level. The output consists of the critical value, the UniLRT test statistic value, and the result, which indicates whether to reject or not reject the null hypothesis.
#' @import stats
#' @importFrom Iso ufit
#' @author Subha Halder
UniLRT <- function(sample_data, significance_level){
  set.seed(456)
  LRT_H0_new <- function(sample_data_list) {
    means <- sapply(sample_data_list, mean)
    sample_sizes <- sapply(sample_data_list, length)
    S <- unlist(sample_data_list)
    mu1 <- mean(S)
    var1 <- sapply(1:length(sample_data_list), function(i) (sum((sample_data_list[[i]] - means[i])^2)) / sample_sizes[i])
    u1 <- sample_sizes / var1

    repeat {
      new_mu1 <- (sum(u1 * means)) / sum(u1)
      new_var1 <- sapply(1:length(sample_data_list), function(i) (sum((sample_data_list[[i]] - new_mu1)^2)) / sample_sizes[i])
      new_u1 <- sample_sizes / new_var1

      if (max(abs(new_mu1 - mu1)) <= 0.0000001) {
        break  # Exit the loop if the difference is less than epsilon
      }

      u1 <- new_u1
      mu1 <- new_mu1
      var1 <- new_var1
    }

    return(var1)
  }

  LRT_H1_new <- function(sample_data_list) {
    n <- sapply(sample_data_list, length)
    mu0 <- sapply(sample_data_list, mean)
    var0 <- sapply(1:length(sample_data_list), function(i) (sum((sample_data_list[[i]] - mu0[i])^2)) / n[i])
    w0 <- n / var0
    repeat {
      new_mu0 <- ufit(y = mu0, imode = ceiling((length(sample_data_list))/2), w = w0)[[2]]
      new_var0 <- sapply(1:length(sample_data_list), function(i) (sum((sample_data_list[[i]] - new_mu0[i])^2)) / n[i])
      new_w0 <- n / new_var0

      if (max(abs(new_mu0 - mu0)) <= 0.0000001) {
        break  # Exit the loop if the difference is less than epsilon
      }

      w0 <- new_w0
      mu0 <- new_mu0
      var0 <- new_var0
    }

    return(var0)
  }

  sample_data <- lapply(sample_data, function(x) x[!is.na(x)])
  num_samples = 20000
  num_datasets <- length(sample_data)
  n <- sapply(sample_data, length)
  lambda_values_star <- numeric(num_samples)
  for (i in 1:num_samples) {
    bootstrap_samples <- lapply(sample_data, function(x) rnorm(n = length(x), mean = 0, sd = sqrt(var(x))))
    V_R_star <- LRT_H1_new(bootstrap_samples) / LRT_H0_new(bootstrap_samples)
    weights <- sapply(1:num_datasets, function(i) V_R_star[i]^(length(bootstrap_samples[[i]]) / 2))
    lambda_values_star[i] <- prod(weights)
  }
  sort_lambda_star <- sort(lambda_values_star)
  quantile_value <- quantile(sort_lambda_star, probs = significance_level)
  V_R <- LRT_H1_new(sample_data) / LRT_H0_new(sample_data)
  weights <- sapply(1:num_datasets, function(i) V_R[i]^(length(sample_data[[i]]) / 2))
  lambda <- prod(weights)
  if (lambda < quantile_value) {
    result <- "Reject null hypothesis"
  } else {
    result <- "Do not reject null hypothesis"
  }
  return(paste("Critical value:", quantile_value, "; UniLRT Test statistic:", lambda, "; Result:", result))
}
