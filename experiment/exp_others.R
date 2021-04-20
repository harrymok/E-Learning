rm(list = ls())

library(dplyr)
library(doParallel)

# setup
CONFIG <- within(list(), {
  # working directory
  dir <- ".."
  # file names to source from
  file_source <- "others.R"
  file_data_generation <- "data_generation.R"
  # file name to save to
  file_save <- "exp_others.RData"
  # data generation setups
  n <- 400
  n_test <- 10000
  p <- 10
  s <- 5
  K <- 2
  SNR <- 1
  ### specification setup
  treatment_free_effect <- TRUE
  propensity_score <- TRUE
  heteroscedasticity <- FALSE
  # methods
  methods <- c("gest", "dwols", "aLearn", "subgroup", "owl", "rwl", "earl", "policyTree")
  # tuning setup
  lambda <- NULL
  lambda_max <- 1
  lambda_min_ratio <- 1e-4
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
####################################################################################################
# training data generation
config <- config_data_generation(CONFIG)
data <- data_generation(config)
####################################################################################################
# training
train <- foreach(method = config[["methods"]], .combine = rbind, .errorhandling = "remove") %do% {
  if(config[[paste0("is_", method)]]) {
    print(paste0("training ", method))
    time <- system.time({
      fit <- get0(method)
      model <- fit(data, config)
    })
    list(method = method, model = model, time = time)
  }
}
####################################################################################################
# testing
data_test <- data_generation(within(config, n <- n_test))
result <- foreach(id = 1:nrow(train), .combine = rbind) %do% {
  method <- train[[id, "method"]]
  model <- train[[id, "model"]]
  time <- train[[id, "time"]]
  predict_function <- get(paste0("predict.", method))
  prediction <- predict_function(model, data_test)
  
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
    t(statistics)
  ))
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
