

# Function to compute pooled mean and standard deviation
pooled_stats <- function(n, means, sds) {
  # Input validation
  if(length(n) != length(means) || length(means) != length(sds)) {
    stop("All input vectors (n, means, sds) must have the same length.")
  }
  
  # Step 1: Compute pooled (weighted) mean
  overall_mean <- sum(n * means) / sum(n)
  
  # Step 2: Compute pooled variance
  within_var <- sum((n - 1) * (sds^2))                      # within-group variance
  between_var <- sum(n * (means - overall_mean)^2)         # between-group variance
  total_var <- (within_var + between_var) / (sum(n) - 1)
  
  # Step 3: Compute pooled standard deviation
  overall_sd <- sqrt(total_var)
  
  # Return as a named list
  return(list(
    mean = overall_mean,
    sd = overall_sd,
    n_total = sum(n)
  ))
}
