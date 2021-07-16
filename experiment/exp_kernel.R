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
  file_save <- "exp.RData"
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
})

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
# data_generation
config <- config_data_generation(CONFIG)
data <- data_generation(config) %>% within({
  # estimate treatment-free effect and propensity score
  prop <- propensity_score(list(X = X, A = A), config)
  prop_A <- extract(prop, A)
  design_A <- format_A(A)[["data"]][["design"]]
})
data_test <- data_generation(within(config, n <- 10000))
test_statistics <- with(data_test, {
  function(A) {
    interaction_effect_predicted <- extract(interaction_effect, as.numeric(A))
    value <- mean(treatment_free_effect + interaction_effect_predicted)
    regret <- value_optimal - value
    misclass <- mean(A != optimal_treatment)
    
    return(c(value = value, regret = regret, misclass = misclass))
  }
})
####################################################################################################
# training
### Kernel Q-Learning
time.kqLearn <- system.time({
  model.kqLearn <- krdLearn(with(config[c("n", "K")], within(data, {
    prop <- matrix(1/K, n, K)
    prop_A <- rep(1/K, n)
  })))
})[["elapsed"]]
print(paste("elapsed time", time.kqLearn))
####################################################################################################
### Kernel RD-Learning
time.krdLearn <- system.time({
  model.krdLearn <- krdLearn(data)
})[["elapsed"]]
print(paste("elapsed time", time.krdLearn))
### Kernel D-Learning
time.kdLearn <- system.time({
  model.kdLearn <- keLearn(data, within(config, method <- "kdLearn"))
})[["elapsed"]]
print(paste("elapsed time", time.kdLearn))
####################################################################################################
### Kernel DR-Learning
time.kdrLearn <- system.time({
  model.kdrLearn <- keLearn(data, within(config, method <- "kdrLearn"))
})[["elapsed"]]
print(paste("elapsed time", time.kdrLearn))
####################################################################################################
### Kernel E-Learning
time.keLearn <- system.time({
  model.keLearn <- keLearn(data, within(config, method <- "keLearn"))
})[["elapsed"]]
print(paste("elapsed time", time.keLearn))
####################################################################################################
# testing
result.kqLearn <- predict.krdLearn(model.kqLearn, data_test)
result.krdLearn <- predict.krdLearn(model.krdLearn, data_test)
result.kdLearn <- predict.keLearn(model.kdLearn, data_test)
result.kdrLearn <- predict.keLearn(model.kdrLearn, data_test)
result.keLearn <- predict.keLearn(model.keLearn, data_test)

methods <- c("kqLearn", "krdLearn", "kdLearn", "kdrLearn", "keLearn")
result <- foreach(method = methods, .combine = rbind) %do% {
  data.frame(
    method = method, 
    time = get(paste0("time.", method)),
    t(test_statistics(get(paste0("result.", method))[["optimal_treatment"]]))
  )
}
####################################################################################################
# testing results
results <- list(result = result, config = config)
results[["result"]]
####################################################################################################
with(CONFIG, save(results, file = file_save))
