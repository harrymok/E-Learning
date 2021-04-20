# model fitting and optimization
library(glmnet)
library(earth)
library(grf)

library(apg)  # accerlated proximal gradient descent
# devtools::install_github("jpvert/apg")

# utilities
library(dplyr)
library(tidyr)
library(ggplot2)
library(doParallel)

# import functions from packages
library(rlang)  # as_function

CONFIG_DEFAULT <- within(list(), { # list of configuration parameters, default
  # working directory
  dir <- ".."
  # file names to source from
  file_basic <- "basic.R"
  file_data_generation <- "data_generation.R"
  file_cosso_interaction <- "cosso_interaction.R"
  # fitting setup
  n_fold_tf <- 10
  n_fold_prop <- 10
  n_fold_sigma2 <- 10
  method_sigma2 <- "regression_forest"
  args_sigma2 <- within(list(), {
    MARS <- within(list(), {
      pmethod <- "cv"
      nfold <- 5
    })
    regression_forest <- within(list(), {
      num.trees <- 200
      tune.parameters <- c(
        "sample.fraction", 
        "mtry", 
        "min.node.size", 
        # "honesty.fraction",
        "honesty.prune.leaves", 
        "alpha", 
        "imbalance.penalty"
      )
      compute.oob.predictions <- FALSE
      num.threads <- 1
    })
    COSSO <- list()
  })
  standardize <- TRUE
  penalty_grouped <- TRUE
  penalty_adaptive <- FALSE
  # tuning setup
  lambda <- NULL
  lambda_max <- 1
  lambda_min_ratio <- 0.01
  lambda_seq_length <- 10
  n_fold_lambda <- 10
  # interface
  verbose <- FALSE
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

# source functions
with(CONFIG, {
  source(paste0(dir, "/", file_basic))
  source(paste0(dir, "/", file_data_generation))
  source(paste0(dir, "/", file_cosso_interaction))
})

# templates
### in-sample/cross fitting and return prediction
fit_template <- function(fit.data, fit.predict)
  # fit.data = data preparation function
  ### input = data, config
  ### output = a list of: data_train, data_predict, config
  # fit.predict = model fitting and prediction function
  ### input = data_train, data_predict, config
  ### output = an n-by-K-matrix prediction
  function(data, config = CONFIG) {
    # config = a list of n_fold_fit, verbose and config for fit.data
    # output = an n-by-K-matrix prediction
    
    # format input
    assign_list(fit.data(data, config), c("data_train", "data_predict", "config"))
    assign_list(config, c("n", "K", "levels_A", "n_fold_fit", "verbose"))
    
    # in-sample/cross fitting and prediction
    prediction <- matrix(nrow = n, ncol = K, dimnames = list(NULL, levels_A))
    if(n_fold_fit >= 2) {
      if(verbose) {
        print(paste0(n_fold_fit, "-fold cross fitting"))
      }
      ### cross fitting
      ids_cross <- cut(sample.int(n), breaks = n_fold_fit, labels = 1:n_fold_fit)
      for(id_fold in 1:n_fold_fit) {
        if(verbose) {
          print(paste0("cross fitting: the ", id_fold, "-th fold"))
        }
        # sample splitting
        data_out <- lapply(data_train, as_function(~ subset(., id_fold != ids_cross)))
        data_in <- lapply(data_predict, as_function(~ subset(., id_fold == ids_cross)))
        # prediction
        prediction[id_fold == ids_cross, ] <- fit.predict(data_out, data_in, config)
      }
    } else { # n_fold_fit == 1
      if(verbose) {
        print("in-sample fitting")
      }
      ### in-sample fitting
      prediction[] <- fit.predict(data_train, data_predict, config)
    }
    return(prediction)
  }
### prediction
predict_template <- function(predict.data, predict.predict)
  # predict.data = data preparation function
  ### input = data, config
  ### output = a list of: data_predict, config
  # predict.predict = prediction function
  ### input = object, data_predict
  ### output = an n-by-K-matrix prediction
  function(object, newdata = NULL, ...) {
    # object = a list of: object elements, config[, prediction, data]
    # output = prediction, optimal_treatment, prediction_optimal[, prediction_treated]
    # ...: unused argument
    
    if(is.null(newdata)) {
      # in-sample prediction
      if(!is.null(object[["prediction"]])) {
        return(object[["prediction"]])
      } else {
        newdata <- object[["data"]]
      }
    }
    
    # format input
    assign_list(predict.data(newdata, object[["config"]]), c("data", "data_predict", "config"))
    assign_list(config, c("n", "K", "levels_A"))
    
    # prediction
    prediction <- matrix(nrow = n, ncol = K, dimnames = list(NULL, levels_A))
    prediction[] <- predict.predict(object, data_predict)
    # optimal treatment
    A_opt_id <- apply(prediction, 1, which.max)
    optimal_treatment <- factor(levels_A[A_opt_id], levels = levels_A)
    # prediction at optimal treatment
    prediction_optimal <- extract(prediction, A_opt_id)
    # prediction at treated
    prediction_treated <- extract(prediction, data[["A_id"]])
    
    return(extract_list(
      as.list(environment()),
      c("prediction", "optimal_treatment", "prediction_optimal", "prediction_treated")
    ))
  }
### value function
value_function <- function(value.predict) {
  # value.predict = prediction function
  ### input = object, data_predict
  ### output = an n-by-K-matrix prediction
  
  value.prescribe <- function(object, data_predict) {
    levels_A <- levels(data_predict[["A"]])
    prediction <- value.predict(object, data_predict)
    A_opt_id <- apply(prediction, 1, which.max)
    return(factor(levels_A[A_opt_id], levels = levels_A))
  }
  
  function(object, data_predict) {
    A_prescribe <- value.prescribe(object, data_predict)
    return(with(data_predict, mean((A == A_prescribe) / prop_A * Y)))
  }
}
### tuning lambda
tune_lambda <- function(tune.fit, tune.value)
  # tune.fit = model fitting function
  ### input = data_train, lambda, config
  ### output = object
  # tune.value = prediction value function
  ### input = object, data_predict
  ### output = prediction value
  function(data_train, data_predict, lambda_seq = NULL, config) {
    if(is.null(lambda_seq)) {
      lambda_seq <- lambda_seq_generation(tune.fit)(data_train, config)
    }
    
    assign_list(extract_list(config, c("n", "n_fold_lambda", "verbose")))
    config <- within(config, verbose <- FALSE)
    
    value <- rep(0, length(lambda_seq))
    names(value) <- ifelse(lambda_seq >= 1e-3, round(lambda_seq, digits = 3), paste0("exp(", round(log(lambda_seq), 3), ")"))
    ids_cross <- cut(sample.int(n), breaks = n_fold_lambda, labels = 1:n_fold_lambda)
    for(id_fold in 1:n_fold_lambda) {
      if(verbose) {
        print(paste0("tuning lambda by ", n_fold_lambda, "-fold cross validation: the ", id_fold, "-th fold"))
      }
      # sample splitting
      data_out <- lapply(data_train, as_function(~ subset(., id_fold != ids_cross)))
      data_in <- lapply(data_predict, as_function(~ subset(., id_fold == ids_cross)))
      # prediction value
      value <- value + sapply(lambda_seq, function(lambda) {
        model <- tune.fit(data_train = data_out, lambda = lambda, config = config)
        return(tune.value(model, data_in))
      })
    }
    value <- value/n_fold_lambda
    id_max <- which.max(value)
    lambda <- lambda_seq[id_max]
    if(verbose) {
      print("cross-validation values")
      print(value)
      print(paste0("tuned lambda = ", lambda))
    }
    return(lambda)
  }
### generating lambda_seq
lambda_seq_generation <- function(lambda_seq.fit)
  # lambda_seq.fit = model fitting function
  ### input = data_train, lambda, config
  ### output = a list containing coefficients
  function(data_train, config) {
    assign_list(extract_list(config, c("lambda_max", "lambda_min_ratio", "lambda_seq_length", "verbose")))
    config <- within(config, verbose <- FALSE)
    
    # searching for the smallest lambda_max such that lambda_seq.fit(data = data, lambda = lambda_max, config = config)[["coefficients"]] are all 0's except for Intercept
    flag <- FALSE
    repeat { 
      coef <- lambda_seq.fit(data_train = data_train, lambda = lambda_max, config = config)[["coefficients"]]
      if(all(abs(coef[-1,]) < 1e-5)) { # search downward
        if(flag) { break }
        lambda_max <- lambda_max / 2
      } else { # search upward
        flag <- TRUE
        lambda_max <- lambda_max * 2
      }
    }
    if(verbose) {
      print(paste0("tuning lambda_seq: lambda_max = ", lambda_max))
    }
    return(lambda_max * exp(seq(0, log(lambda_min_ratio), length.out = lambda_seq_length)))
  }

# estimate treatment-free effect
treatment_free_effect <- with(list(), {
  # define internal functions
  ### define data preparation function
  treatment_free_effect.data <- function(data, config) {
    assign_list(format_data(extract_list(data, c("X", "A", "Y", "prop")), config), 
                c("data", "config"))
    if(config[["standardize"]]) {
      # standardize design_X as in glmnet
      if(is.null(config[["center_X"]]) | is.null(config[["scale_X"]])) {
        # compute the centers and scaling constants for non-intercept columns
        config <- with(data, within(config, {
          center_X <- apply(design_X[, -1], 2, mean)
          scale_X <- apply(design_X[, -1] - rep(1, nrow(design_X)) %*% t(center_X), 
                           2, function(x) sqrt(sum(x^2) / length(x)))
        }))
      }
      data <- with(config, within(data, {
        design_X[, -1] <- scale(design_X[, -1], center_X, scale_X)
      }))
    }
    data <- with(config, within(data, {
      treatment_free_design <- design_X[, -1]  # excluding the intercept
      interaction_design <- (1 - 1/K) * (design_A %x_row% design_X)
    }))
    data_train <- with(data, within(list(), {
      x <- cbind(treatment_free_design, interaction_design)
      y <- Y
      weights <- prop_A^{-1}
    }))
    data_predict <- with(data, {
      interaction_design[] <- 0
      within(list(), {
        newx <- cbind(treatment_free_design, interaction_design)
      })
    })
    return(extract_list(as.list(environment()), c("data_train", "data_predict", "config")))
  }
  ### define model fitting and prediction function
  treatment_free_effect.predict <- function(data_train, data_predict, config) {
    model <- with(data_train, cv.glmnet(x, y, weights, standardize = FALSE))
    return(with(data_predict, predict(model, newx, s = "lambda.min")[, 1]))
  }
  # treatment free effect function
  function(data, config = CONFIG, ...) {
    # ...: unused argument
    
    if(config[["verbose"]]) {
      print("fitting treatment-free effect function")
    }
    config <- within(config, n_fold_fit <- n_fold_tf)
    fit_template(treatment_free_effect.data, treatment_free_effect.predict)(data, config)[, 1]
  }
})
####################################################################################################
##### UNIT TEST for "treatment_free_effect"
# # K = 2
# config <- config_data_generation(within(CONFIG, {
#   K <- 2
#   verbose <- TRUE
# }))
# data <- within(data_generation(config = config), {
#   prop <- propensity_score
# })
# tf <- treatment_free_effect(data, within(config, n_fold_tf <- 1))
# tf_crossfit <- treatment_free_effect(data, config)
# qqplot(data[["treatment_free_effect"]], tf); abline(a = 0, b = 1)
# qqplot(data[["treatment_free_effect"]], tf_crossfit); abline(a = 0, b = 1)
# qqplot(tf, tf_crossfit); abline(a = 0, b = 1)
####################################################################################################
# # K = 3
# config <- config_data_generation(within(CONFIG, {
#   K <- 3
#   verbose <- TRUE
# }))
# data <- within(data_generation(config = config), {
#   prop <- propensity_score
# })
# tf <- treatment_free_effect(data, within(config, n_fold_tf <- 1))
# tf_crossfit <- treatment_free_effect(data, config)
# qqplot(data[["treatment_free_effect"]], tf); abline(a = 0, b = 1)
# qqplot(data[["treatment_free_effect"]], tf_crossfit); abline(a = 0, b = 1)
# qqplot(tf, tf_crossfit); abline(a = 0, b = 1)
####################################################################################################

# estimate propensity score
propensity_score <- with(list(), {
  # define internal functions
  ### define data preparation function
  propensity_score.data <- function(data, config) {
    assign_list(format_data(extract_list(data, c("X", "A")), config), 
                c("data", "config"))
    if(config[["standardize"]]) {
      # standardize design_X as in glmnet
      if(is.null(config[["center_X"]]) | is.null(config[["scale_X"]])) {
        # compute the centers and scaling constants for non-intercept columns
        config <- with(data, within(config, {
          center_X <- apply(design_X[, -1], 2, mean)
          scale_X <- apply(design_X[, -1] - rep(1, nrow(design_X)) %*% t(center_X), 
                           2, function(x) sqrt(sum(x^2) / length(x)))
        }))
      }
      data <- with(config, within(data, {
        design_X[, -1] <- scale(design_X[, -1], center_X, scale_X)
      }))
    }
    data_train <- with(data, within(list(), {
      x <- design_X[, -1]  # excluding the intercept
      y <- A
    }))
    data_predict <- with(data, within(list(), {
      newx <- design_X[, -1]  # excluding the intercept
    }))
    return(extract_list(as.list(environment()), c("data_train", "data_predict", "config")))
  }
  ### define model fitting and prediction function
  propensity_score.predict <- function(data_train, data_predict, config) {
    if(config[["K"]] == 2) {
      model <- with(data_train, cv.glmnet(x, y, family = "binomial", standardize = FALSE))
      pred_1 <- with(data_predict, predict(model, newx, s = "lambda.min", type = "response")[, 1])
      return(cbind(1 - pred_1, pred_1))
    } else { # K > 2
      model <- with(data_train, cv.glmnet(x, y, family = "multinomial", standardize = FALSE, type.multinomial = "grouped"))
      return(with(data_predict, predict(model, newx, s = "lambda.min", type = "response")[, , 1]))
    }
  }
  # propensity score function
  function(data, config = CONFIG, ...) {
    # ...: unused argument
    
    if(config[["verbose"]]) {
      print("fitting propensity score function")
    }
    config <- within(config, n_fold_fit <- n_fold_prop)
    fit_template(propensity_score.data, propensity_score.predict)(data, config)
  }
})
####################################################################################################
##### UNIT TEST for "propensity_score"
# ### K = 2
# config <- config_data_generation(within(CONFIG, {
#   K <- 2
#   verbose <- TRUE
# }))
# data <- data_generation(config = config)
# prop <- propensity_score(data, config = within(config, n_fold_prop <- 1))
# prop_crossfit <- propensity_score(data, config)
# qqplot(data[["propensity_score"]][, 1], prop[, 1]); abline(a = 0, b = 1)
# qqplot(data[["propensity_score"]][, 1], prop_crossfit[, 1]); abline(a = 0, b = 1)
# qqplot(prop[, 1], prop_crossfit[, 1]); abline(a = 0, b = 1)
####################################################################################################
# ### K = 3
# config <- config_data_generation(within(CONFIG, {
#   K <- 3
#   verbose <- TRUE
# }))
# data <- data_generation(config = config)
# prop <- propensity_score(data, within(config, n_fold_prop <- 1))
# prop_crossfit <- propensity_score(data, config)
# qqplot(data[["propensity_score"]], prop); abline(a = 0, b = 1)
# qqplot(data[["propensity_score"]], prop_crossfit); abline(a = 0, b = 1)
# qqplot(prop, prop_crossfit); abline(a = 0, b = 1)
####################################################################################################
# ### estimate treatmemt-free effect with estimated propensity scores
# config <- config_data_generation(within(CONFIG, {
#   K <- 3
#   verbose <- TRUE
# }))
# data <- data_generation(config = config)
# tf <- treatment_free_effect(within(data, {
#   prop <- propensity_score(data, config)
# }), config)
# qqplot(data[["treatment_free_effect"]], tf); abline(a = 0, b = 1)
# ### misspecified propensity score
# config[["propensity_score"]] <- FALSE
# data <- within(data_generation(config = config), {
#   propensity_score_estimate <- propensity_score(data, config)
# })
# with(data, qqplot(propensity_score, propensity_score_estimate)); abline(a = 0, b = 1)
# ### estimate treatment-free effect with estimated propensity scores
# tf <- treatment_free_effect(within(data, {
#   prop <- propensity_score_estimate
# }), config)
# tf_correct_prop <- treatment_free_effect(within(data, {
#   prop <- propensity_score
# }), config)
# qqplot(data[["treatment_free_effect"]], tf); abline(a = 0, b = 1)
# qqplot(data[["treatment_free_effect"]], tf_correct_prop); abline(a = 0, b = 1)
# qqplot(tf_correct_prop, tf); abline(a = 0, b = 1)
####################################################################################################

# l1-PLS
### define internal functions
l1PLS_internal <- within(list(), {
  ### define data preparation function
  l1PLS.data <- function(data, config) {
    assign_list(format_data(extract_list(data, c("X", "A", "Y")), config),
                c("data", "config"))
    data_train <- with(data, within(list(), {
      x <- model.matrix(~ A * ., data.frame(X, A))[, -1]  # [X, A, X * A]
      y <- Y
    }))
    data_predict <- with(c(data, config), {
      X <- X[rep(1:n, K), ]
      A <- factor(rep(levels_A, each = n), levels = levels_A)
      within(list(), {
        newx <- model.matrix(~ A * ., data.frame(X, A))[, -1]
      })
    })
    return(extract_list(as.list(environment()), c("data", "data_train", "data_predict", "config")))
  }
  ### define prediction function
  l1PLS.predict <- function(object, data_predict) {
    return(with(data_predict, predict(object[["model"]], newx, s = "lambda.min")[, 1]))
  }
})
### model fitting function
l1PLS <- with(l1PLS_internal, {
  # Q-Learning with LASSO: E(Y|X,A) ~ 1 + X + A + X * A
  function(data, config = CONFIG, ...) {
    # ...: unused argument
    
    # format input
    assign_list(l1PLS.data(data, config), c("data", "data_train", "config"))
    
    # fit l1-PLS
    model <- with(data_train, cv.glmnet(x, y))
    
    output <- extract_list(as.list(environment()), c("model", "data", "config"))
    # coefficients
    coefficients <- coefficients.l1PLS(output)
    # in-sample prediction
    prediction <- predict.l1PLS(output, data)
    return(c(output, list(coefficients = coefficients, prediction = prediction)))
  }
})
### prediction function
predict.l1PLS <- with(l1PLS_internal, {
  predict_template(l1PLS.data, l1PLS.predict)
})
### extracting coefficients
coefficients.l1PLS <- function(object, ...) {
  # object = a list of: object elements, config
  # output = an n-by-K-matrix prediction
  # ...: unused argument
  
  assign_list(object, c("model", "config"))
  
  data_predict <- with(config, {
    design_X <- diag(p)[, -1]  # excluding the intercept
    design_X <- design_X[rep(1:p, K), ]
    A <- factor(rep(levels_A, each = p), levels = levels_A)
    within(list(), {
      newx <- model.matrix(~ A * ., data.frame(design_X, A))[, -1]
    })
  })
  prediction <- with(config, matrix(nrow = p, ncol = K, dimnames = list(NULL, levels_A)))
  prediction[] <- with(data_predict, predict(model, newx, s = "lambda.min")[, 1])
  
  with(config, {
    design_X <- diag(p); design_X[, 1] <- 1  # intercept column
    coefficients <- solve(design_X, prediction)
    coefficients <- coefficients - rowMeans(coefficients)
    rownames(coefficients) <- names_X
    return(coefficients)
  })
}
####################################################################################################
##### UNIT TEST for "l1PLS"
# # K = 2
# data <- data_binary()
# model_l1PLS <- l1PLS(data)
# round(model_l1PLS[["coefficients"]], digits = 2)
# 
# # K = 3
# data <- data_multiple()
# model_l1PLS <- l1PLS(data)
# round(model_l1PLS[["coefficients"]], digits = 2)
####################################################################################################

# q-COSSO
### define internal functions
qCOSSO_internal <- within(list(), {
  ### define data preparation function
  qCOSSO.data <- function(data, config) {
    assign_list(format_data(extract_list(data, c("X", "A", "Y")), config),
                c("data", "config"))
    data_train <- with(data, within(list(), {
      X1 <- X
      X2 <- data.frame(A = A)
      Y <- Y
    }))
    data_predict <- with(c(data, config), within(list(), {
      X1 <- X[rep(1:n, K), ]
      X2 <- data.frame(A = factor(rep(levels_A, each = n), levels = levels_A))
    }))
    config <- within(config, {
      levels_X1 <- levels_X
      range_X1 <- range_X
      levels_X2 <- list(A = levels_A)
    })
    return(extract_list(as.list(environment()), c("data", "data_train", "data_predict", "config")))
  }
  qCOSSO.predict <- function(object, data_predict) {
    return(predict.cosso_interaction(object[["model"]], data_predict))
  }
})
### model fitting function
qCOSSO <- with(qCOSSO_internal, {
  # Q-Learning with COSSO: SS-ANOVA model for Q-Learning: E(Y|X,A) ~ 1 + X + A + X * A
  function(data, config = CONFIG, ...) {
    # ...: unused argument
    
    # format input
    assign_list(qCOSSO.data(data, config), c("data", "data_train", "config"))
    
    # fit q-COSSO
    model <- within(cosso_interaction(data_train, config), rm("data"))
    
    output <- extract_list(as.list(environment()), c("model", "data", "config"))
    # in-sample prediction
    prediction <- predict.qCOSSO(output, data)
    return(c(output, list(prediction = prediction)))
  }
})
### prediction function
predict.qCOSSO <- with(qCOSSO_internal, {
  predict_template(qCOSSO.data, qCOSSO.predict)
})
####################################################################################################
##### UNIT TEST for "qCOSSO"
# config <- within(list(), {
#   n <- 400
#   p <- 100
#   range_X <- replicate(p, c(-3, 3), simplify = FALSE)
#   names(range_X) <- paste0("X", 1:p)
# })
### K = 2
# data <- data_binary(config)
# time_qCOSSO <- system.time({
#   model_qCOSSO <- qCOSSO(data, config)
# })
# ### in-sample prediction
# table(truth = data[["optimal_treatment"]],
#       predict = with(model_qCOSSO, prediction[["optimal_treatment"]]))
# qqplot(x = data[["interaction_effect"]],
#        y = with(model_qCOSSO, prediction[["prediction"]][, "1"]),
#        xlab = "truth", ylab = "predict"); abline(0, 1)
# ### testing prediction
# data <- data_binary(within(config, n <- 10000))
# pred_qCOSSO <- predict.qCOSSO(model_qCOSSO, data)
# table(truth = data[["optimal_treatment"]],
#       predict = pred_qCOSSO[["optimal_treatment"]])
# qqplot(x = data[["interaction_effect"]],
#        y = pred_qCOSSO[["prediction"]][, "1"],
#        xlab = "truth", ylab = "predict"); abline(0, 1)
####################################################################################################
### K = 3
# config <- within(list(), {
#   n <- 400
#   p <- 10
#   K <- 3
#   range_X <- replicate(p, c(-3, 3), simplify = FALSE)
#   names(range_X) <- paste0("X", 1:p)
# })
# data <- data_multiple(config)
# time_qCOSSO <- system.time({
#   model_qCOSSO <- qCOSSO(data, config)
# })
# ### in-sample prediction
# table(truth = data[["optimal_treatment"]],
#       predict = with(model_qCOSSO, prediction[["optimal_treatment"]]))
# qqplot(x = data[["interaction_effect"]],
#        y = with(model_qCOSSO, prediction[["prediction"]]),
#        xlab = "truth", ylab = "predict"); abline(0, 1)
# ### testing prediction
# data <- data_multiple(within(config, { n <- 10000 }))
# pred_qCOSSO <- predict.qCOSSO(model_qCOSSO, data)
# table(truth = data[["optimal_treatment"]],
#       predict = pred_qCOSSO[["optimal_treatment"]])
# qqplot(x = data[["interaction_effect"]],
#        y = pred_qCOSSO[["prediction"]],
#        xlab = "truth", ylab = "predict"); abline(0, 1)
####################################################################################################
### experiment
# config <- config_data_generation(CONFIG)
# data <- data_generation(config = config)
# time_qCOSSO <- system.time({
#   model_qCOSSO <- qCOSSO(data, config)
# })
# ### in-sample prediction
# table(truth = data[["optimal_treatment"]],
#       predict = with(model_qCOSSO, prediction[["optimal_treatment"]]))
# qqplot(x = data[["interaction_effect"]],
#        y = with(model_qCOSSO, prediction[["prediction"]]),
#        xlab = "truth", ylab = "predict"); abline(0, 1)
# ### testing prediction
# data <- data_generation(within(config, { n <- 10000 }))
# pred_qCOSSO <- predict.qCOSSO(model_qCOSSO, data)
# table(truth = data[["optimal_treatment"]],
#       predict = pred_qCOSSO[["optimal_treatment"]])
# qqplot(x = data[["interaction_effect"]],
#        y = pred_qCOSSO[["prediction"]],
#        xlab = "truth", ylab = "predict"); abline(0, 1)
####################################################################################################

# D-Learning
dLearn_internal <- within(list(), {
  dLearn.data <- function(data, config) {
    assign_list(format_data(extract_list(data, c("X", "A", "Y", "prop")), config),
                c("data", "config"))
    if(config[["standardize"]]) {
      # standardize design_X as in glmnet
      if(is.null(config[["center_X"]]) | is.null(config[["scale_X"]])) {
        # compute the centers and scaling constants for non-intercept columns
        config <- with(data, within(config, {
          center_X <- apply(design_X[, -1], 2, mean)
          scale_X <- apply(design_X[, -1] - rep(1, nrow(design_X)) %*% t(center_X), 
                           2, function(x) sqrt(sum(x^2) / length(x)))
        }))
      }
      data <- with(config, within(data, {
        design_X[, -1] <- scale(design_X[, -1], center_X, scale_X)
      }))
    }
    data_train <- extract_list(data, c("design_X", "design_A", "prop_A", "Y"))
    data_predict <- extract_list(data, c("design_X", "A", "prop_A", "Y"))
    return(extract_list(as.list(environment()), c("data", "data_train", "data_predict", "config")))
  }
  dLearn.predict <- function(object, data_predict) {
    return(data_predict[["design_X"]] %*% object[["coefficients"]])
  }
  dLearn.value <- value_function(dLearn.predict)
  dLearn.fit <- function(data_train, lambda, config) {
    assign_list(config, c("n", "p", "K", "names_X", "groups_X", "levels_A", "Omega", "penalty_grouped", "penalty_adaptive"))
    
    coefficients <- matrix(nrow = p, ncol = K, dimnames = list(names_X, levels_A))
    if(lambda == 0) { # without regularization
      model_lm <- with(data_train, lm(K*Y*design_A ~ design_X - 1, weights = prop_A^{-1}))
      coef_code <- model_lm[["coefficients"]]
    } else { # with row-wise group-LASSO regularization
      # create groups for prox.grouplasso
      if(penalty_grouped) {
        # groups <- split(1:(p*(K-1)), groups_X)
        groups <- split(1:(p*(K-1)), 1:p)
      } else {
        groups <- as.list(seq(p*(K-1)))
      }
      # compute group weights if penalty_adaptive
      if(penalty_adaptive) {
        coef <- dLearn.fit(data_train, 1e-5, within(config, penalty_adaptive <- FALSE))[["coefficients"]]
        coef_code <- coef %*% Omega
        group_norms <- sapply(groups, function(ids) sqrt(sum(coef_code[ids]^2)))
        group_weights <- lambda / group_norms
      } else {
        group_weights <- lambda
      }
      optim_apg <- apg(dLearn.estimating_function(data_train, config), prox.grouplasso, p*(K-1),
                       list(groups = groups, groupweights = group_weights, QUIET = TRUE))
      coef_code <- matrix(optim_apg[["x"]], nrow = p, ncol = K-1)
    }
    coefficients[] <- (1 - 1/K) * coef_code %*% t(Omega)
    
    return(list(coefficients = coefficients))
  }
  dLearn.estimating_function <- function(data_train, config) {
    # gradient function of (2*n*K)^{-1} * || prop_A^{-1/2} * (K*Y*design_A - design_X %*% coef_matrix) ||_F^2
    
    assign_list(data_train, c("design_X", "design_A", "prop_A", "Y"))
    assign_list(config, c("n", "p", "K"))
    
    Corr <- t(prop_A^{-1} * design_X) %*% (Y * design_A)
    Gram <- t(prop_A^{-1} * design_X) %*% (1/K * design_X)
    
    function(coef_vec, opts = NULL) {
      # input: vectorized coefficients
      # output: vectorized gradient
      coef_matrix <- matrix(coef_vec, nrow = p, ncol = K-1)
      return(n^{-1} * as.numeric(Gram %*% coef_matrix - Corr))
    }
  }
})
### model fitting function
dLearn <- with(dLearn_internal, {
  # D-Learning (Qi et al., 2020)
  function(data, lambda = config[["lambda"]], config = CONFIG, ...) {
    # ...: unused argument
    
    # format input
    assign_list(dLearn.data(data, config), c("data", "data_train", "data_predict", "config"))
    
    # tune lambda
    if(is.null(lambda)) {
      lambda <- tune_lambda(dLearn.fit, dLearn.value)(data_train, data_predict, , config)
    } else if(length(lambda) > 1) {
      lambda <- tune_lambda(dLearn.fit, dLearn.value)(data_train, data_predict, lambda, config)
    }
    
    # fit D-Learning
    model <- dLearn.fit(data_train, lambda, config)
    
    output <- c(model, list(lambda = lambda, data = data, config = config))
    # in-sample prediction
    prediction <- predict.dLearn(output, data)
    return(c(output, list(prediction = prediction)))
  }
})
### prediction function
predict.dLearn <- with(dLearn_internal, {
  predict_template(dLearn.data, dLearn.predict)
})
####################################################################################################
#### UNIT TEST for "dLearn"
# config <- config_data_generation(within(CONFIG, {
#   verbose <- TRUE
#   coefficients <- matrix(0, nrow = p + 1, ncol = K, dimnames = list(c("Intercept", paste0("X", 1:p)), 1:K))
#   coefficients[1:(s+1), 1:3] <- t(replicate(s + 1, c(-1, 0, 1)))
# }))
# data <- within(data_generation(config), prop <- propensity_score)
# ### specified a lambda
# model <- dLearn(data, 0, config); model[["coefficients"]]
# ### specified a tuning sequence of lambdas
# model <- dLearn(data, exp(seq(0, log(0.01), length.out = 10)), config); model[["coefficients"]]
# ### unspecified tuning sequence
# model <- dLearn(data, , config); model[["coefficients"]]
# ### in-sample prediction
# pred <- predict.dLearn(model)
# table(truth = data[["optimal_treatment"]], predict = pred[["optimal_treatment"]])
# qqplot(x = data[["interaction_effect"]], y = pred[["prediction"]], xlab = "truth", ylab = "predict"); abline(0, 1)
# ### testing prediction
# data <- data_generation(within(config, n <- 10000))
# pred <- predict.dLearn(model, data)
# table(truth = data[["optimal_treatment"]], predict = pred[["optimal_treatment"]])
# qqplot(x = data[["interaction_effect"]], y = pred[["prediction"]], xlab = "truth", ylab = "predict"); abline(0, 1)
####################################################################################################

# RD-Learning
rdLearn_internal <- within(list(), {
  rdLearn.data <- function(data, config) {
    assign_list(format_data(extract_list(data, c("X", "A", "Y", "prop")), config),
                c("data", "config"))
    if(config[["standardize"]]) {
      # standardize design_X as in glmnet
      if(is.null(config[["center_X"]]) | is.null(config[["scale_X"]])) {
        # compute the centers and scaling constants for non-intercept columns
        config <- with(data, within(config, {
          center_X <- apply(design_X[, -1], 2, mean)
          scale_X <- apply(design_X[, -1] - rep(1, nrow(design_X)) %*% t(center_X), 
                           2, function(x) sqrt(sum(x^2) / length(x)))
        }))
      }
      data <- with(config, within(data, {
        design_X[, -1] <- scale(design_X[, -1], center_X, scale_X)
      }))
    }
    data_train <- extract_list(with(config, within(data, {
      design_interaction <- (1 - 1/K) * (design_A %x_row% design_X)
    })), c("design_interaction", "prop_A", "Y"))
    data_predict <- extract_list(data, c("design_X", "A", "prop_A", "Y"))
    return(extract_list(as.list(environment()), c("data", "data_train", "data_predict", "config")))
  }
  rdLearn.predict <- function(object, data_predict) {
    return(data_predict[["design_X"]] %*% object[["coefficients"]])
  }
  rdLearn.value <- value_function(rdLearn.predict)
  rdLearn.fit <- function(data_train, lambda, config) {
    assign_list(config, c("n", "p", "K", "names_X", "groups_X", "levels_A", "Omega", "penalty_grouped", "penalty_adaptive"))
    
    coefficients <- matrix(nrow = p, ncol = K, dimnames = list(names_X, levels_A))
    if(lambda == 0) { # without regularization
      model_lm <- with(data_train, lm(Y ~ design_interaction - 1, weights = prop_A^{-1}))
      coef_code <- matrix(model_lm[["coefficients"]], nrow = p, ncol = K-1)
    } else { # with row-wise group-LASSO regularization
      # create groups for prox.grouplasso: grouping by rows except for the interecpt
      if(penalty_grouped) {
        # groups <- split(1:(p*(K-1)), groups_X)
        groups <- split(1:(p*(K-1)), 1:p)
      } else {
        groups <- as.list(seq(p*(K-1)))
      }
      # compute group weights if penalty_adaptive
      if(penalty_adaptive) {
        coef <- rdLearn.fit(data_train, 1e-5, within(config, penalty_adaptive <- FALSE))[["coefficients"]]
        coef_code <- coef %*% Omega
        group_norms <- sapply(groups, function(ids) sqrt(sum(coef_code[ids]^2)))
        group_weights <- lambda / group_norms
      } else {
        group_weights <- lambda
      }
      optim_apg <- apg(rdLearn.estimating_function(data_train, config), prox.grouplasso, p*(K-1),
                       list(groups = groups, groupweights = group_weights, QUIET = TRUE))
      coef_code <- matrix(optim_apg[["x"]], nrow = p, ncol = K-1)
    }
    coefficients[] <- (1 - 1/K) * coef_code %*% t(Omega)
    
    return(list(coefficients = coefficients))
  }
  rdLearn.estimating_function <- function(data_train, config) {
    # gradient function of (2*n)^{-1} * || prop_A^{-1/2} * (Y - design_interaction %*% coef) ||_2^2
    
    assign_list(data_train, c("design_interaction", "prop_A", "Y"))
    assign_list(config, "n")
    
    Corr <- t(prop_A^{-1} * design_interaction) %*% Y
    Gram <- t(prop_A^{-1} * design_interaction) %*% design_interaction
    
    function(coef_vec, opts = NULL) {
      # input: vectorized coefficients
      # output: vectorized gradient
      return(n^{-1} * as.numeric(Gram %*% coef_vec - Corr))
    }
  }
})
### model fitting function
rdLearn <- with(rdLearn_internal, {
  # RD-Learning (Meng and Qiao, 2020)
  function(data, lambda = config[["lambda"]], config = CONFIG, ...) {
    # ...: unused argument
    
    # format input
    assign_list(rdLearn.data(data, config), c("data", "data_train", "data_predict", "config"))
    
    # tune lambda
    if(is.null(lambda)) {
      lambda <- tune_lambda(rdLearn.fit, rdLearn.value)(data_train, data_predict, , config)
    } else if(length(lambda) > 1) {
      lambda <- tune_lambda(rdLearn.fit, rdLearn.value)(data_train, data_predict, lambda, config)
    }
    
    # fit RD-Learning
    model <- rdLearn.fit(data_train, lambda, config)
    
    output <- c(model, list(lambda = lambda, data = data, config = config))
    # in-sample prediction
    prediction <- predict.rdLearn(output, data)
    return(c(output, list(prediction = prediction)))
  }
})
### prediction function
predict.rdLearn <- with(rdLearn_internal, {
  predict_template(rdLearn.data, rdLearn.predict)
})
####################################################################################################
##### UNIT TEST for "rdLearn"
# config <- config_data_generation(within(CONFIG, {
#   verbose <- TRUE
#   coefficients <- matrix(0, nrow = p + 1, ncol = K, dimnames = list(c("Intercept", paste0("X", 1:p)), 1:K))
#   coefficients[1:(s+1), 1:3] <- t(replicate(s + 1, c(-1, 0, 1)))
# }))
# data <- within(data_generation(config), prop <- propensity_score)
# ### specified a lambda
# model <- rdLearn(data, 0, config); model[["coefficients"]]
# ### specified a tuning sequence of lambdas
# model <- rdLearn(data, exp(seq(0, log(0.01), length.out = 10)), config); model[["coefficients"]]
# ### unspecified tuning sequence
# model <- rdLearn(data, , config); model[["coefficients"]]
# ### in-sample prediction
# pred <- predict.rdLearn(model)
# table(truth = data[["optimal_treatment"]], predict = pred[["optimal_treatment"]])
# qqplot(x = data[["interaction_effect"]], y = pred[["prediction"]], xlab = "truth", ylab = "predict"); abline(0, 1)
# ### testing prediction
# data <- data_generation(within(config, n <- 10000))
# pred <- predict.rdLearn(model, data)
# table(truth = data[["optimal_treatment"]], predict = pred[["optimal_treatment"]])
# qqplot(x = data[["interaction_effect"]], y = pred[["prediction"]], xlab = "truth", ylab = "predict"); abline(0, 1)
####################################################################################################

# E-Learning
eLearn_internal <- within(list(), {
  eLearn.data <- function(data, config) {
    assign_list(format_data(extract_list(data, c("X", "A", "Y", "prop", "sigma2")), config),
                c("data", "config"))
    if(config[["standardize"]]) {
      # standardize design_X as in glmnet
      if(is.null(config[["center_X"]]) | is.null(config[["scale_X"]])) {
        # compute the centers and scaling constants for non-intercept columns
        config <- with(data, within(config, {
          center_X <- apply(design_X[, -1], 2, mean)
          scale_X <- apply(design_X[, -1] - rep(1, nrow(design_X)) %*% t(center_X), 
                           2, function(x) sqrt(sum(x^2) / length(x)))
        }))
      }
      data <- with(config, within(data, {
        design_X[, -1] <- scale(design_X[, -1], center_X, scale_X)
      }))
    }
    data_train <- extract_list(with(config, within(data, {
      if(is.null(sigma2)) {
        sigma2 <- matrix(1, nrow = n, ncol = K, dimnames = list(NULL, levels_A))
      }
      if(any(is.infinite(sigma2))) {
        stop("infinite values in sigma2")
      }
      if(any(is.na(sigma2))) {
        stop("NAs in sigma2")
      }
      design_interaction <- (1 - 1/K) * (design_A %x_row% design_X)
      V_A <- t(Omega) %*% Omega
      instrument <- t(sapply(1:n, function(id) {
        gradient_X <- diag(K-1) %x% design_X[id, ]
        instrument_A <- design_A[id, ] / prop_A[id]
        V <- t(Omega) %*% diag(sigma2[id, ]/prop[id, ]) %*% Omega
        V_inv <- with(svd(V), {
          d_inv <- ifelse(d >= 1e-5, d^{-1}, 0)
          D_inv <- diag(d_inv, nrow = K-1, ncol = K-1)
          return(u %*% D_inv %*% t(v))
        })
        return(gradient_X %*% V_A %*% V_inv %*% instrument_A)
      }))
    })), c("design_interaction", "instrument", "Y"))
    data_predict <- extract_list(data, c("design_X", "A", "prop_A", "Y"))
    return(extract_list(as.list(environment()), c("data", "data_train", "data_predict", "config")))
  }
  eLearn.predict <- function(object, data_predict) {
    return(data_predict[["design_X"]] %*% object[["coefficients"]])
  }
  eLearn.value <- value_function(eLearn.predict)
  eLearn.fit <- function(data_train, lambda, config) {
    assign_list(config, c("n", "p", "K", "names_X", "groups_X", "levels_A", "Omega", "penalty_grouped", "penalty_adaptive"))
    
    coefficients <- matrix(nrow = p, ncol = K, dimnames = list(names_X, levels_A))
    # create groups for prox.grouplasso: grouping by rows except for the interecpt
    if(penalty_grouped) {
      # groups <- split(1:(p*(K-1)), groups_X)
      groups <- split(1:(p*(K-1)), 1:p)
    } else {
      groups <- as.list(seq(p*(K-1)))
    }
    # compute group weights if penalty_adaptive
    if(penalty_adaptive) {
      coef <- eLearn.fit(data_train, 1e-5, within(config, penalty_adaptive <- FALSE))[["coefficients"]]
      coef_code <- coef %*% Omega
      group_norms <- sapply(groups, function(ids) sqrt(sum(coef_code[ids]^2)))
      group_weights <- lambda / group_norms
    } else {
      group_weights <- lambda
    }
    optim_apg <- apg(eLearn.estimating_function(data_train, config), prox.grouplasso, p*(K-1),
                     list(groups = groups, groupweights = group_weights, QUIET = TRUE))
    coef_code <- matrix(optim_apg[["x"]], nrow = p, ncol = K-1)
    coefficients[] <- (1 - 1/K) * coef_code %*% t(Omega)
    
    return(list(coefficients = coefficients))
  }
  eLearn.estimating_function <- function(data_train, config) {
    # solving n^{-1} * t(instrument) * (Y - design_interaction %*% coef) = 0
    
    assign_list(data_train, c("design_interaction", "instrument", "Y"))
    assign_list(config, "n")
    
    Corr <- t(instrument) %*% Y
    Gram <- t(instrument) %*% design_interaction
    
    function(coef_vec, opts = list()) {
      # input: vectorized coefficients
      # output: vectorized gradient
      return(n^{-1} * as.numeric(Gram %*% coef_vec - Corr))
    }
  }
  eLearn.variance_function <- with(list(), {
    # define internal functions
    ### compute residual
    eLearn.residual <- function(data, config) {
      data <- extract_list(data, c("X", "A", "Y", "prop"))
      if(config[["verbose"]]) {
        print("computing E-Learning residuals")
        config[["verbose"]] <- FALSE
      }
      # set sigma2 = 1
      data[["sigma2"]] <- with(format_data(data, config)[["config"]], {
        sigma2 <- matrix(1, nrow = n, ncol = K, dimnames = list(NULL, levels_A))
      })
      # fit E-Learning and obtain in-sample prediction at treated
      prediction <- eLearn(data, config[["lambda"]], config)[[c("prediction", "prediction_treated")]]
      return(with(data, Y - prediction))
    }
    ### define data preparation function
    variance_function.data <- function(data, config) {
      data[["residual"]] <- eLearn.residual(data, config)
      data <- extract_list(data, c("X", "A", "residual"))
      
      return(switch(
        config[["method_sigma2"]],
        MARS = variance_function.data.MARS(data, config),
        regression_forest = variance_function.data.regression_forest(data, config),
        COSSO = variance_function.data.COSSO(data, config)
      ))
    }
    variance_function.data.MARS <- function(data, config) {
      assign_list(format_data(extract_list(data, c("X", "A", "residual")), config),
                  c("data", "config"))
      data_train <- with(data, within(list(), {
        x <- model.matrix(~ A * ., data.frame(X, A))[, -1]  # [X, A, X * A]
        y <- residual^2
      }))
      data_predict <- with(c(data, config), {
        X <- X[rep(1:n, K), ]
        A <- factor(rep(levels_A, each = n), levels = levels_A)
        within(list(), {
          newdata <- model.matrix(~ A * ., data.frame(X, A))[, -1]
        })
      })
      return(extract_list(as.list(environment()), c("data_train", "data_predict", "config")))
    }
    variance_function.data.regression_forest <- function(data, config) {
      assign_list(format_data(extract_list(data, c("X", "A", "residual")), config),
                  c("data", "config"))
      data_train <- with(data, within(list(), {
        X <- model.matrix(~ A * ., data.frame(X, A))[, -1]  # [X, A, X * A]
        Y <- residual^2
      }))
      data_predict <- with(c(data, config), {
        X <- X[rep(1:n, K), ]
        A <- factor(rep(levels_A, each = n), levels = levels_A)
        within(list(), {
          newdata <- model.matrix(~ A * ., data.frame(X, A))[, -1]
        })
      })
      return(extract_list(as.list(environment()), c("data_train", "data_predict", "config")))
    }
    variance_function.data.COSSO <- function(data, config) {
      assign_list(format_data(extract_list(data, c("X", "A", "residual")), config),
                  c("data", "config"))
      data_train <- with(data, within(list(), {
        X1 <- X
        X2 <- data.frame(A = A)
        Y <- residual^2
      }))
      data_predict <- with(c(data, config), within(list(), {
        X1 <- X[rep(1:n, K), ]
        X2 <- data.frame(A = factor(rep(levels_A, each = n), levels = levels_A))
      }))
      config <- within(config, {
        levels_X1 <- levels_X
        range_X1 <- range_X
        levels_X2 <- list(A = levels_A)
      })
      return(extract_list(as.list(environment()), c("data_train", "data_predict", "config")))
    }
    ### define prediction functions
    variance_function.predict <- function(data_train, data_predict, config) {
      sigma2 <- switch(
        config[["method_sigma2"]],
        MARS = variance_function.predict.MARS(data_train, data_predict, config),
        regression_forest = variance_function.predict.regression_forest(data_train, data_predict, config),
        COSSO = variance_function.predict.COSSO(data_train, data_predict, config)
      )
      sigma2 <- pmax(0, sigma2)
      return(sigma2)
    }
    variance_function.predict.MARS <- function(data_train, data_predict, config) {
      model <- do.call(earth, c(data_train, config[[c("args_sigma2", "MARS")]]))
      return(with(data_predict, predict(model, newdata)))
    }
    variance_function.predict.regression_forest <- function(data_train, data_predict, config) {
      model <- do.call(regression_forest, c(data_train, config[[c("args_sigma2", "regression_forest")]]))
      return(with(data_predict, predict(model, newdata))[["predictions"]])
    }
    variance_function.predict.COSSO <- function(data_train, data_predict, config) {
      model <- cosso_interaction(data_train, config)
      return(predict.cosso_interaction(model, data_predict))
    }
    # variance function
    function(data, config) {
      is_method_sigma2 <- FALSE
      if(!is.null(config[["method_sigma2"]])) 
        if(config[["method_sigma2"]] %in% c("MARS", "regression_forest", "COSSO")) {
          is_method_sigma2 <- TRUE
        }
      if(!is_method_sigma2) {
        warning("incorrect method_sigma2; set method_sigma2 = regression_forest")
        config[["method_sigma2"]] <- "regression_forest"
      }
      
      if(config[["verbose"]]) {
        print(paste0("fitting variance function by ", config[["method_sigma2"]]))
      }
      config <- within(config, n_fold_fit <- n_fold_sigma2)
      fit_template(variance_function.data, variance_function.predict)(data, config)
    }
  })
})
### model fitting function
eLearn <- with(eLearn_internal, {
  # E-Learning
  function(data, lambda = config[["lambda"]], config = CONFIG, ...) {
    # ...: unused argument
    
    if(is.null(data[["sigma2"]])) {
      # estimate sigma2
      data[["sigma2"]] <- eLearn.variance_function(data, config)
    }
    
    # format input
    assign_list(eLearn.data(data, config), c("data", "data_train", "data_predict", "config"))
    
    # tune lambda
    if(is.null(lambda)) {
      lambda <- tune_lambda(eLearn.fit, eLearn.value)(data_train, data_predict, , config)
    } else if(length(lambda) > 1) {
      lambda <- tune_lambda(eLearn.fit, eLearn.value)(data_train, data_predict, lambda, config)
    }
    
    # fit E-Learning
    model <- eLearn.fit(data_train, lambda, config)
    
    output <- c(model, list(lambda = lambda, data = data, config = config))
    # in-sample prediction
    prediction <- predict.eLearn(output, data)
    return(c(output, list(prediction = prediction)))
  }
})
### prediction function
predict.eLearn <- with(eLearn_internal, {
  predict_template(eLearn.data, eLearn.predict)
})
####################################################################################################
##### UNIT TEST for "eLearn"
# config <- config_data_generation(within(CONFIG, {
#   method_sigma2 <- "COSSO"
#   verbose <- TRUE
#   coefficients <- matrix(0, nrow = p + 1, ncol = K, dimnames = list(c("Intercept", paste0("X", 1:p)), 1:K))
#   coefficients[1:(s+1), 1:3] <- t(replicate(s + 1, c(-1, 0, 1)))
# }))
# data <- within(data_generation(config), prop <- propensity_score)
# ### specified a lambda
# model <- eLearn(data, 0, config); model[["coefficients"]]
# ### specified a tuning sequence of lambdas
# model <- eLearn(data, exp(seq(0, log(0.01), length.out = 10)), config); model[["coefficients"]]
# ### unspecified tuning sequence
# model <- eLearn(data, , config); model[["coefficients"]]
# ### in-sample prediction
# pred <- predict.eLearn(model)
# table(truth = data[["optimal_treatment"]], predict = pred[["optimal_treatment"]])
# qqplot(x = data[["interaction_effect"]], y = pred[["prediction"]], xlab = "truth", ylab = "predict"); abline(0, 1)
# ### testing prediction
# data <- within(data_generation(within(config, n <- 10000)), prop <- propensity_score)
# pred <- predict.eLearn(model, data)
# table(truth = data[["optimal_treatment"]], predict = pred[["optimal_treatment"]])
# qqplot(x = data[["interaction_effect"]], y = pred[["prediction"]], xlab = "truth", ylab = "predict"); abline(0, 1)
####################################################################################################

# # E-RD-Ensemble
# config_erdLearn <- function(config = CONFIG) {
#   config <- within(config, {
#     names_assign <- c("alpha_seq_by", "n_fold_alpha", "verbose")
#   })
#   retain_name <- with(config, unique(c(
#     names_assign, "alpha", "names_assign",
#     names(config_rdLearn(config)),
#     names(config_eLearn(config))
#   )))
#   config <- config[retain_name]
#   return(config)
# }
# predict.erdLearn <- predict_template
# erdLearn <- function(data, alpha = CONFIG[["alpha"]], config = CONFIG, ...) {
#   # data: a list of: X, A, Y[, prop, sigma2, levels_X, range_X]
#   # Assumption: Y has been subtracted from a treatment-free effect estimate
#   # Assumption: if prop is unspecified, then propensity scores are computed by frequencies in A
#   # alpha: ensemble weight
#   # ...: unused argument
# 
#   config <- config_erdLearn(config)
#   assign_list(config)
# 
#   if(verbose) {
#     print("E-RD-Ensemble configurations:")
#     print(unlist(within(config, rm(names_assign))))
#   }
# 
#   # format input data
#   data <- format_data(data[c("X", "A", "Y", "prop", "sigma2")])
#   data_train <- data[c("X", "A", "Y", "prop", "sigma2")]
#   data_others <- data[c("levels_X", "range_X")]
#   data_pred <- data[c("X", "A", "Y", "prop_A")]
#   assign_list(data, c("n", "K", "levels_A"))
# 
#   # tune alpha
#   if(is.null(alpha)) {
#     alpha_seq <- seq(0, 1, alpha_seq_by)
#     if(verbose) {
#       print(paste0("tuning alpha by ", n_fold_alpha, "-fold cross validation"))
#     }
#     predictor_rdLearn <- predictor_eLearn <- matrix(0, nrow = n, ncol = K)
#     ids_cross <- cut(sample.int(n), breaks = n_fold_alpha, labels = 1:n_fold_alpha)
#     for(id_fold in 1:n_fold_alpha) {
#       if(verbose) {
#         print(paste0("tuning: the ", id_fold, "-th fold"))
#       }
# 
#       data_out <- lapply(data_train, as_function(~ subset(., id_fold != ids_cross)))
#       data_out <- c(data_out, data_others)
#       data_in  <- lapply(data_pred,  as_function(~ subset(., id_fold == ids_cross)))
# 
#       model_rdLearn <- rdLearn(data_out, , within(config, { verbose <- FALSE }))
#       model_eLearn  <- eLearn(data_out, ,  within(config, { verbose <- FALSE }))
#       pred_rdLearn  <- predict.rdLearn(model_rdLearn, data_in)
#       pred_eLearn   <- predict.eLearn(model_eLearn, data_in)
# 
#       predictor_rdLearn[id_fold == ids_cross, ] <- pred_rdLearn[["interaction_effect"]]
#       predictor_eLearn[id_fold == ids_cross, ]  <- pred_eLearn[["interaction_effect"]]
#     }
#     value_cv <- sapply(alpha_seq, function(alpha) {
#       predictor_ensemble <- alpha * predictor_rdLearn + (1 - alpha) * predictor_eLearn
#       A_pred_id <- apply(predictor_ensemble, 1, which.max)
#       A_pred <- factor(levels_A[A_pred_id], levels = levels_A)
#       return(with(data_in, mean((A == A_pred) / prop_A * Y)))
#     })
#     id_max <- which.max(value_cv)
#     alpha <- alpha_seq[id_max]
#     if(verbose) {
#       names(value_cv) <- alpha_seq
#       print("cross-validation values")
#       print(value_cv)
#       print(paste0("tuned alpha = ", alpha))
#     }
#   }
# 
#   ### training
#   output <- within(list(), {
#     model <- list(
#       rdLearn = rdLearn(data, , config),
#       eLearn = eLearn(data, , config)
#     )
#     data <- with(model, eLearn[["data"]])
#     alpha <- alpha
#     coefficients <- with(model, {
#       alpha * rdLearn[["coefficients"]] + (1 - alpha) * eLearn[["coefficients"]]
#     })
#   })
#   return(within(output, {
#     # in-sample prediction
#     prediction <- predict.erdLearn(output, data)
#   }))
# }
####################################################################################################
##### UNIT TEST for "erdLearn"
# config <- within(CONFIG, {
#   n_fold_lambda <- 2
#   n_fold_alpha <- 2
#   verbose <- TRUE
# })
# data <- data_multiple()
# ### rdLearn
# time_rdLearn <- system.time({
#   coef_rdLearn <- rdLearn(data, , config)[["coefficients"]]
# })
#
# ### eLearn
# time_eLearn <- system.time({
#   coef_eLearn <- eLearn(data, , config)[["coefficients"]]
# })
#
# ### erdLearn with alpha = 0.5
# time_erdLearn.5 <- system.time({
#   coef_erdLearn.5 <- erdLearn(data, 0.5, config)[["coefficients"]]
# })
#
### erdLearn with tuned alpha
# time_erdLearn <- system.time({
#   coef_erdLearn <- erdLearn(data, , config)[["coefficients"]]
# })
####################################################################################################
