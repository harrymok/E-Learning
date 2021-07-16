CONFIG_DEFAULT <- within(list(), { # list of configuration parameters, default
  n <- 400  # training sample size
  p <- 10   # covariate dimension regardless the intercept
  s <- 5    # effective covariate dimension regardless the intercept
  K <- 3    # number of treatment arms
  SNR <- 1  # SNR = interaction_effect^2 / (treatment_free_effect^2 + variance)
  ### specification setup
  heteroscedasticity <- FALSE
})

if(is.null(get0("left_full_join_list"))) {
  # set list elements according to another list
  left_full_join_list <- function(List_left = list(), List_right = list()) {
    # return a list including all elements from List_left, and elements from List_right without a match in List_Left 
    names_set <- setdiff(names(List_right), names(List_left))
    List_left[names_set] <- List_right[names_set]
    return(List_left)
  }
} 

# overwrite CONFIG with defaults if unspecified
CONFIG <- left_full_join_list(get0("CONFIG"), CONFIG_DEFAULT)

### data generation for nonlinear experiments
unit_vec <- function(p) {
  # generate uniformly from the p-dimensional unit sphere
  v <- rnorm(p)
  return(v / sqrt(sum(v^{2})))
}
config_data_generation <- function(config = CONFIG) {
  return(within(config, {
    s <- pmax(s, K)
    p <- pmax(p, s)
    range_X <- replicate(p, c(-3, 3), simplify = FALSE)
    names(range_X) <- paste0("X", 1:p)
    levels_A <- 1:K
    if(is.null(config[["coefficients"]])) {
      coefficients <- replicate(K, unit_vec(s + 1))
      coefficients <- rbind(coefficients, matrix(0, nrow = p - s, ncol = K))
      ### multiply by sqrt(SNR)
      coefficients <- sqrt(SNR) * coefficients
      ### subtract treatment-free effect: row-sum-to-zero
      coefficients <-  (1 - 1/K)^{-1/2} * (coefficients - rowMeans(coefficients)) 
      rownames(coefficients) <- c("Intercept", paste0("X", 1:p))
      colnames(coefficients) <- 1:K
    }
  }))
}
data_generation <- function(config = CONFIG) {
  return(with(config_data_generation(config), within(list(), {
    # covariates
    X <- matrix(rnorm(n*p), nrow = n, ncol = p)
    
    # model specifications
    Exp_X <- exp(sqrt(2) * X[,1:K])
    ### treatment-free effects
    treatment_free_effect <- sqrt(K) * rowMeans(X[,1:K]^2)
    
    ### propensity scores
    propensity_score <- exp(X[,1:K] / 2)
    propensity_score <- propensity_score / rowSums(propensity_score)
    
    ### variances
    if(heteroscedasticity) {
      ### heteroscedasticity
      variance <- Exp_X^2
    } else {
      ### homoscedasticity
      variance <- matrix(1, nrow = n, ncol = K)
    }
    
    # treatment assignments
    A <- apply(propensity_score, 1, function(prob) sample(K, 1, prob = prob))
    
    # interaction effects, optimal treatments and interaction effects at optimal
    interaction_effect <- cbind(1, X^2) %*% coefficients
    optimal_treatment <- apply(interaction_effect, 1, which.max)
    interaction_effect_optimal <- extract(interaction_effect, optimal_treatment)
    
    # propensity scores, variances and interaction effects at treated
    colnames(interaction_effect) <- colnames(propensity_score) <- colnames(variance) <- levels_A
    propensity_score_treated <- extract(propensity_score, A)
    variance_treated <- extract(variance, A)
    interaction_effect_treated <- extract(interaction_effect, A)
    
    # format X, A, optimal_treatment
    X <- data.frame(X)
    A <- factor(A, levels = levels_A)
    optimal_treatment <- factor(optimal_treatment, levels = levels_A)
    
    # outcomes
    mean_treated <- treatment_free_effect + interaction_effect_treated
    mean_optimal <- treatment_free_effect + interaction_effect_optimal
    value_optimal <- mean(mean_optimal)
    error <- rnorm(n)
    Y <- mean_treated + sqrt(variance_treated) * error
  })))
}
