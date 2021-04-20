rm(list = ls())

library(dplyr)
library(doParallel)

# setup
CONFIG <- within(list(), {
  # working directory
  dir <- ".."
  # file names to source from
  file_source <- "eLearn.R"
  file_data_generation <- "data_generation.R"
  # file name to save to
  file_save <- "experiment.RData"
  # data generation setups
  n <- 400
  n_test <- 10000
  p <- 10
  s <- 5
  K <- 3
  SNR <- 1
  ### specification setup
  treatment_free_effect <- TRUE
  propensity_score <- TRUE
  heteroscedasticity <- FALSE
  # training setups
  methods_sigma2 <- c("MARS", "regression_forest", "COSSO")
  methods <- c(
    "l1PLS", 
    "dLearn", 
    "rdLearn", 
    paste0("eLearn_", methods_sigma2),
    "eLearn_oracle"
  )
  ### fitting setup
  n_fold_tf <- 10
  n_fold_prop <- 10
  n_fold_sigma2 <- 1
  method_sigma2 <- "regression_forest"
  ### tuning setup
  lambda <- NULL
  lambda_max <- 1
  lambda_min_ratio <- 0.01
  lambda_seq_length <- 10
  n_fold_lambda <- 10
  # interface
  verbose <- FALSE
})
CONFIG <- c(CONFIG, with(CONFIG, {
  is_methods <- replicate(length(methods), list(TRUE))
  names(is_methods) <- paste0("is_", methods)
  is_methods
}))

# import setup from command line --args and override CONFIG
args <- commandArgs(trailingOnly = TRUE); args
for(input in names(CONFIG)) {
  id <- grep(paste0("^", input, "="), args)
  if(length(id)) { # input from args
    CONFIG[[input]] <- eval(parse(text = gsub(paste0("^", input, "="), "", args[id])))
  }
  print(paste0(input, " = ", CONFIG[[input]]))
}

# source functions
with(CONFIG, {
  source(paste0(dir, "/", file_source))
  source(paste0(dir, "/", file_data_generation))
})

# training data generation
config <- config_data_generation(CONFIG)
data <- data_generation(config) %>% within({
  # estimate treatment-free effect and propensity score
  prop <- propensity_score(list(X = X, A = A), config)
  tf <- treatment_free_effect(list(X = X, A = A, Y = Y, prop = prop), config)
  # optimal working variance for E-Learning
  sigma2_optimal <- (treatment_free_effect - tf)^2 + variance
})

# training
####################################################################################################
### l1-PLS
if(config[["is_l1PLS"]]) {
  time_l1PLS <- system.time({
    model_l1PLS <- l1PLS(data, config)
  })
  print(paste("elapsed time", time_l1PLS[["elapsed"]]))
}
####################################################################################################
### D-Learning
if(config[["is_dLearn"]]) {
  time_dLearn <- system.time({
    model_dLearn <- dLearn(data, config[["lambda"]], config)
  })
  print(paste("elapsed time", time_dLearn[["elapsed"]]))
}
####################################################################################################
### RD-Learning
if(config[["is_rdLearn"]]) {
  time_rdLearn <- system.time({
    model_rdLearn <- rdLearn(within(data, Y <- Y - tf), config[["lambda"]], config)
  })
  print(paste("elapsed time", time_rdLearn[["elapsed"]]))
}
####################################################################################################
### E-Learning
##### E-learning with oracle working variance
if(config[["is_eLearn_oracle"]]) {
  time_eLearn_oracle <- system.time({
    model_eLearn_oracle <- eLearn(within(data, { Y <- Y - tf; sigma2 <- sigma2_optimal }), 
                                  config[["lambda"]], config)
  })
  print(paste("elapsed time", time_eLearn_oracle[["elapsed"]]))
}
####################################################################################################
##### E-learning with estimated working variance
for(method in config[["methods_sigma2"]]) {
  cat(paste0(paste(rep("#", 100), collapse = "")), "\n")
  if(config[[paste0("is_eLearn_", method)]]) {
    cat(paste0("##### E-Learning with ", method, "\n"))
    time_eLearn <- system.time({
      model_eLearn <- eLearn(within(data, { Y <- Y - tf; sigma2 <- NULL }), config[["lambda"]], 
                             within(config, method_sigma2 <- method))
    })
    assign(paste0("model_eLearn_", method), model_eLearn)
    assign(paste0("time_eLearn_", method), time_eLearn)
    print(paste("elapsed time", time_eLearn[["elapsed"]]))
  }
}
####################################################################################################

# testing prediction
### testing data generation
data_test <- within(data_generation(within(config, n <- n_test)), { prop <- propensity_score; sigma2 <- variance })
result <- foreach(method = config[["methods"]], .combine = rbind) %do% {
  if(config[[paste0("is_", method)]]) {
    model <- get0(paste0("model_", method))
    time <- get0(paste0("time_", method))
    predict_function <- get(paste0("predict.", gsub("([[:alnum:]]*)([_][[:print:]]*)", "\\1", method)))
    prediction <- predict_function(model, data_test)
    
    RMSE <- sqrt(mean((model[["coefficients"]] - config[["coefficients"]])^2))
    
    statistics <- with(data_test, {
      A <- prediction[["optimal_treatment"]]
      interaction_effect_predicted <- extract(interaction_effect, as.numeric(A))
      value <- mean(treatment_free_effect + interaction_effect_predicted)
      regret <- value_optimal - value
      misclass <- mean(A != optimal_treatment)
      
      return(c(value = value, regret = regret, misclass = misclass))
    })
    
    return(data.frame(
      method = method,
      time = time[["elapsed"]],
      t(statistics),
      RMSE = RMSE
    ))
  }
}
####################################################################################################
results <- within(list(), {
  result <- result
  config <- config
})
####################################################################################################
# testing results
results[["result"]]
####################################################################################################
with(CONFIG, save(results, file = file_save))
