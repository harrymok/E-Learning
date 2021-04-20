CONFIG_DEFAULT <- within(list(), { # list of configuration parameters, default
  # working directory
  dir <- ".."
  # file names to source from
  file_basic <- "basic.R"
  # tuning setup
  lambda <- NULL
  lambda_max <- 1
  lambda_min_ratio <- 1e-4
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
})

####################################################################################################
### generic internal functions
internal_functions <- within(list(), {
  # prediction from model coefficients
  predict_from_coefficients <- function(object, newdata = NULL, ...) {
    # object = a list of: coefficients, config[, prediction, data]
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
    assign_list(format_data(newdata, object[["config"]]), "data")
    
    # prediction
    prediction <- with(c(object, data), design_X %*% coefficients)
    # optimal treatment
    A_opt_id <- apply(prediction, 1, which.max)
    optimal_treatment <- with(config, factor(levels_A[A_opt_id], levels = levels_A))
    # prediction at optimal treatment
    prediction_optimal <- extract(prediction, A_opt_id)
    # prediction at treated
    prediction_treated <- extract(prediction, data[["A_id"]])
    
    return(extract_list(
      as.list(environment()),
      c("prediction", "optimal_treatment", "prediction_optimal", "prediction_treated")
    ))
  }
})
####################################################################################################
### internal functions for methods in DTRreg
DTRreg_internal <- within(internal_functions, {
  # fit model
  DTRreg <- function(data, config = CONFIG, ...) {
    # ...: unused argument
    
    # format input
    assign_list(format_data(extract_list(data, c("X", "A", "Y")), config),
                c("data", "config"))
    data_train <- with(c(data, config), {
      data.frame(design_X[, -1], A = ifelse(A == levels_A[2], 1, 0), Y)
    })
    
    # fit model
    formula_X <- with(config, as.formula(paste("~", paste(names_X[-1], collapse = " + "))))
    ### blip model
    blip.mod <- list(formula_X)
    ### treatment model
    treat.mod <- list(update(formula_X, A ~ .))
    ### treatment-free model
    tf.mod <- list(formula_X)
    ### fit DTRreg
    model <- DTRreg::DTRreg(
      Y, blip.mod, treat.mod, tf.mod, data = data_train,
      method = config[["method"]], verbose = config[["verbose"]]
    )
    psi <- model[["psi"]][[1]]
    coef <- cbind(0, psi)
    coef <- coef - rowMeans(coef)
    dimnames(coef) <- with(config, list(names_X, levels_A))
    output <- list(coefficients = coef, config = config)
    
    # in-sample prediction
    prediction <- with(data, predict.DTRreg(output, data))
    
    return(c(output, list(prediction = prediction)))
  }
  # prediction
  predict.DTRreg <- predict_from_coefficients
})
####################################################################################################
### G-Estimation (Robins, 2004)
gest <- with(DTRreg_internal, {
  function(data, config = CONFIG, ...) {
    DTRreg(data, within(config, method <- "gest"))
  }
})
predict.gest <- with(DTRreg_internal, predict.DTRreg)
####################################################################################################
### dWOLS (Wallace and Moodie, 2015)
dwols <- with(DTRreg_internal, {
  function(data, config = CONFIG, ...) {
    DTRreg(data, within(config, method <- "dwols"))
  }
})
predict.dwols <- with(DTRreg_internal, predict.DTRreg)
####################################################################################################
### internal functions for methods in personalized
personalized_internal <- within(internal_functions, {
  # estimate propensity score
  prop.func <- function(x, trt) {
    propens.model <- glmnet::cv.glmnet(y = trt, x = x, family = "binomial")
    prop <- predict(propens.model, newx = x, s = "lambda.min", type = "response")[ ,1]
    return(prop)
  }
  # estimate treatment-free effect
  augment.func <- within(list(), {
    ### method = "a_learning"
    a_learning <- function(x, y) {
      data <- data.frame(x, y)
      
      lmod <- glmnet::cv.glmnet(y = y, x = x)
      pred  <- predict(lmod, x, s = "lambda.min")[, 1]
      
      return(pred)
    }
    ### method = "weighting"
    weighting <- function(x, y, trt) {
      data <- data.frame(x, y, trt)
      xm <- model.matrix(y~trt*x-1, data = data)
      
      lmod <- glmnet::cv.glmnet(y = y, x = xm)
      ## get predictions when trt = 1
      data$trt <- 1
      xm <- model.matrix(y~trt*x-1, data = data)
      preds_1  <- predict(lmod, xm, s = "lambda.min")[, 1]
      
      ## get predictions when trt = -1
      data$trt <- -1
      xm <- model.matrix(y~trt*x-1, data = data)
      preds_n1  <- predict(lmod, xm, s = "lambda.min")[, 1]
      
      ## return predictions averaged over trt
      return(0.5 * (preds_1 + preds_n1))
    }
  })
  # fit model
  fit.subgroup <- function(data, config = CONFIG, ...) {
    # ...: unused argument
    
    # format input
    assign_list(format_data(extract_list(data, c("X", "A", "Y")), config),
                c("data", "config"))
    data <- with(config, within(data, {
      A_01 <- ifelse(A == levels_A[2], 1, 0)
    }))
    
    # fit model
    method <- config[["method"]]
    model <- with(data, personalized::fit.subgroup(
      x = design_X[, -1], y = Y, trt = A_01,
      propensity.func = prop.func, augment.func = augment.func[[method]],
      method = method, loss = "sq_loss_lasso", nfolds = config[["n_fold_lambda"]]
    ))
    coefficients <- model[["coefficients"]][-1, 1]
    coef <- cbind(- coefficients, coefficients)
    dimnames(coef) <- with(config, list(names_X, levels_A))
    output <- list(coefficients = coef, config = config)
    
    # in-sample prediction
    prediction <- with(data, predict.fit.subgroup(output, data))
    
    return(c(output, list(prediction = prediction)))
  }
  # prediction
  predict.fit.subgroup <- predict_from_coefficients
})
####################################################################################################
### A-Learning (Lu et al., 2013; Shi et al., 2018)
aLearn <- with(personalized_internal, {
  function(data, config = CONFIG, ...) {
    fit.subgroup(data, within(config, method <- "a_learning"))
  }
})
predict.aLearn <- with(personalized_internal, predict.fit.subgroup)
####################################################################################################
### subgroup identification (Tian et al., 2014; Chen et al., 2017)
subgroup <- with(personalized_internal, {
  function(data, config = CONFIG, ...) {
    fit.subgroup(data, within(config, method <- "weighting"))
  }
})
predict.subgroup <- with(personalized_internal, predict.fit.subgroup)
####################################################################################################
### internal functions for methods in DynTxRegime
cv.glmnet.formula <- get("cv.glmnet.formula", environment(glmnetUtils::cv.glmnet))
DynTxRegime_internal <- within(internal_functions, {
  # generate lambda_seq for tuning
  lambda_seq_generation <- function(lambda_seq.fit)
    # lambda_seq.fit = model fitting function
    ### input = data, config
    ### output = a list containing coefficients
    function(data_train, config) {
      assign_list(extract_list(config, c("lambda_max", "lambda_min_ratio", "lambda_seq_length", "verbose")))
      config <- within(config, verbose <- FALSE)
      
      # searching for the smallest lambda_max such that lambda_seq.fit(data = data, lambda = lambda_max, config = config)[["coefficients"]] are all 0's except for Intercept
      flag <- FALSE
      repeat {
        config[["lambda"]] <- lambda_max
        coef <- lambda_seq.fit(data = data, config = config)[["coefficients"]]
        if(all(abs(coef[-1,]) < 1e-2)) { # search downward
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
  # prediction
  predict.DynTxRegime <- predict_from_coefficients
})
####################################################################################################
### OWL (Zhao et al., 2012)
owl <- with(DynTxRegime_internal, {
  function(data, config = CONFIG, ...) {
    # ...: unused argument
    
    # format input
    assign_list(format_data(extract_list(data, c("X", "A", "Y")), config),
                c("data", "config"))
    
    # generate lambda sequence for tuning
    if(is.null(config[["lambda"]])) {
      config[["lambda"]] <- lambda_seq_generation(owl)(data, config)
    }
    
    # fit model
    formula_X <- with(config, as.formula(paste("~", paste(names_X[-1], collapse = " + "))))
    ### propensity model
    moPropen <- modelObj::buildModelObj(
      model = formula_X,
      solver.method = "cv.glmnet.formula",
      solver.args = list(use.model.frame = TRUE, family = "binomial"),
      predict.method = "predict",
      predict.args = list(s = "lambda.min", type = "response")
    )
    ### fit owl
    model <- with(data, DynTxRegime::owl(
      moPropen = moPropen,
      data = data.frame(design_X[, -1], A), reward = Y, txName = "A",
      regime = formula_X, surrogate = "hinge", kernel = "linear", 
      lambdas = 1/config[["lambda"]], cvFolds = config[["n_fold_lambda"]], verbose = config[["verbose"]]
    ))
    regimeCoef <- DynTxRegime::regimeCoef(model)
    coef <- cbind(- regimeCoef, regimeCoef)
    dimnames(coef) <- with(config, list(names_X, levels_A))
    output <- list(coefficients = coef, config = config)
    
    # in-sample prediction
    prediction <- with(data, predict.earl(output, data))
    
    return(c(output, list(prediction = prediction)))
  }
})
predict.owl <- with(DynTxRegime_internal, predict.DynTxRegime)
####################################################################################################
### RWL (Zhou et al., 2017)
rwl <- with(DynTxRegime_internal, {
  function(data, config = CONFIG, ...) {
    # ...: unused argument
    
    # generate lambda sequence for tuning
    if(is.null(config[["lambda"]])) {
      config[["lambda"]] <- lambda_seq_generation(rwl)(data, config)
    }
    
    # format input
    assign_list(format_data(extract_list(data, c("X", "A", "Y")), config),
                c("data", "config"))
    
    # fit model
    formula_X <- with(config, as.formula(paste("~", paste(names_X[-1], collapse = " + "))))
    ### propensity model
    moPropen <- modelObj::buildModelObj(
      model = formula_X,
      solver.method = "cv.glmnet.formula",
      solver.args = list(use.model.frame = TRUE, family = "binomial"),
      predict.method = "predict",
      predict.args = list(s = "lambda.min", type = "response")
    )
    ### outcome models
    moMain <- modelObj::buildModelObj(
      model = formula_X,
      solver.method = "cv.glmnet.formula",
      solver.args = list(use.model.frame = TRUE),
      predict.method = "predict",
      predict.args = list(s = "lambda.min")
      
    )
    ### fit rwl
    model <- with(data, DynTxRegime::rwl(
      moPropen = moPropen, moMain = moMain,
      data = data.frame(design_X[, -1], A), reward = Y, txName = "A",
      regime = formula_X, kernel = "linear", 
      lambdas = config[["lambda"]], cvFolds = config[["n_fold_lambda"]], verbose = config[["verbose"]]
    ))
    regimeCoef <- DynTxRegime::regimeCoef(model)
    coef <- cbind(- regimeCoef, regimeCoef)
    dimnames(coef) <- with(config, list(names_X, levels_A))
    output <- list(coefficients = coef, config = config)
    
    # in-sample prediction
    prediction <- with(data, predict.earl(output, data))
    
    return(c(output, list(prediction = prediction)))
  }
})
predict.rwl <- with(DynTxRegime_internal, predict.DynTxRegime)
####################################################################################################
### EARL (Zhao et al., 2019)
earl <- with(DynTxRegime_internal, {
  function(data, config = CONFIG, ...) {
    # ...: unused argument
    
    # format input
    assign_list(format_data(extract_list(data, c("X", "A", "Y")), config),
                c("data", "config"))
    
    # generate lambda sequence for tuning
    if(is.null(config[["lambda"]])) {
      config[["lambda"]] <- lambda_seq_generation(earl)(data, config)
    }
    
    # fit model
    formula_X <- with(config, as.formula(paste("~", paste(names_X[-1], collapse = " + "))))
    ### propensity model
    moPropen <- modelObj::buildModelObj(
      model = formula_X,
      solver.method = "cv.glmnet.formula",
      solver.args = list(use.model.frame = TRUE, family = "binomial"),
      predict.method = "predict",
      predict.args = list(s = "lambda.min", type = "response")
    )
    ### outcome models
    moMain <- modelObj::buildModelObj(
      model = formula_X,
      solver.method = "cv.glmnet.formula",
      solver.args = list(use.model.frame = TRUE),
      predict.method = "predict",
      predict.args = list(s = "lambda.min")
      
    )
    moCont <- modelObj::buildModelObj(
      model = formula_X,
      solver.method = "cv.glmnet.formula",
      solver.args = list(use.model.frame = TRUE),
      predict.method = "predict",
      predict.args = list(s = "lambda.min")
    )
    ### fit earl
    model <- with(data, DynTxRegime::earl(
      moPropen = moPropen, moMain = moMain, moCont = moCont,
      data = data.frame(design_X[, -1], A), response = Y, txName = "A",
      regime = formula_X, surrogate = "hinge", kernel = "linear", 
      lambdas = config[["lambda"]], cvFolds = config[["n_fold_lambda"]], verbose = config[["verbose"]]
    ))
    regimeCoef <- DynTxRegime::regimeCoef(model)
    coef <- cbind(- regimeCoef, regimeCoef)
    dimnames(coef) <- with(config, list(names_X, levels_A))
    output <- list(coefficients = coef, config = config)
    
    # in-sample prediction
    prediction <- with(data, predict.earl(output, data))
    
    return(c(output, list(prediction = prediction)))
  }
})
predict.earl <- with(DynTxRegime_internal, predict.DynTxRegime)
####################################################################################################
### doubly robust policy learning (Athey and Wager, 2021; Zhou et al., 2018)
policyTree <- function(data, config = CONFIG, ...) {
  # ...: unused argument
  
  # format input
  assign_list(format_data(extract_list(data, c("X", "A", "Y")), config),
              c("data", "config"))
  
  # fit model
  multi.forest <- with(data, policytree::multi_causal_forest(X = design_X[, -1], Y = Y, W = A))
  DR.scores <- policytree::double_robust_scores(multi.forest)
  model <- with(data, policytree::policy_tree(design_X[, -1], DR.scores, depth = 2))
  output <- list(model = model, config = config)
  
  # in-sample prediction
  prediction <- with(data, predict.policyTree(output, data))
  
  return(c(output, list(prediction = prediction)))
}
predict.policyTree <- function(object, newdata = NULL, ...) {
  # object = output of policyTree
  # output = optimal_treatment
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
  assign_list(format_data(newdata, object[["config"]]), "data")
  
  # prediction
  predict.policy_tree <- get("predict.policy_tree", environment(policytree::policy_tree))
  optimal_treatment <- with(data, predict(object[["model"]], design_X[, -1]))
  
  return(list(optimal_treatment = optimal_treatment))
}