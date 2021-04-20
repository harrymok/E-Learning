# cosso_interaction: 
### fit COSSO of the SS-ANOVA model: Y ~ 1 + X1 + X2 + X1 * X2, family = "Gaussian"

library(cosso)  # depend on cosso.Gaussian, bigGram; bigGram, genK.cat, rescale are overwritten
library(abind)  # "abind" function
library(doParallel)  # "foreach" function

CONFIG_DEFAULT <- within(list(), { # list of configuration parameters, default
  # working directory
  dir <- "."
  # file names to source from
  file_basic <- "basic.R"
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

# source basic functions
with(CONFIG, {
  source(paste0(dir, "/", file_basic))
})

# overwrite functions from cosso
### overwrite cosso::twostep.Gaussian
unlockBinding("twostep.Gaussian", environment(cosso.Gaussian))
assign("twostep.Gaussian", envir = environment(cosso.Gaussian), value = {
  function (Gramat1, Gramat2, y, wt, lambda, mm)  {
    n <- length(y)
    nbasis <- dim(Gramat1)[2]
    d <- length(wt)
    cb0 <- sspline(Gramat1, Gramat2, y, 1/(wt^2), lambda)
    c0 <- cb0[-1]
    b0 <- cb0[1]
    G1 <- matrix(0, n, d)
    G2 <- matrix(0, nbasis, d)
    for (j in 1:d) {
      G1[, j] = Gramat1[, , j] %*% c0 * (wt[j]^(-2))
      G2[, j] = Gramat2[, , j] %*% c0 * (wt[j]^(-2))
    }
    dvec <- 2 * t(G1) %*% (y - b0) - n * lambda * t(G2) %*% c0
    Dmat <- 2 * t(G1) %*% G1
    constant <- svd(Dmat)[["d"]][1]
    Dmat <- Dmat + constant * 1e-5 * diag(d)
    Amat <- rbind(diag(1, d), diag(-1, d))
    bvec <- c(rep(0, d), rep(-mm, d))
    theta1 <- solve.QP(Dmat/constant, dvec/constant, t(Amat), bvec)[[1]]
    theta1[theta1 < 1e-08] <- 0
    cb1 <- sspline(Gramat1, Gramat2, y, theta1/(wt^2), lambda)
    result <- list(coefs = cb1[-1], intercept = cb1[1], theta = theta1)
    return(result)
  }
})
lockBinding("twostep.Gaussian", environment(cosso.Gaussian))

### overwrite cosso::genK.cat: generate the kernel matrix for an L-level factor
genK.cat <- function(x1, x2, L = length(unique(c(x1, x2)))) {
  if(L < length(unique(c(x1, x2)))) {
    stop(paste0("x1, x2 have more than L = ", L, " levels"))
  }
  
  n1 <- length(x1)
  n2 <- length(x2)
  x1 <- rep(x1, times = n2)
  x2 <- rep(x2, each = n1)
  
  return(matrix((x1 == x2) - 1/L, nrow = n1, ncol = n2))
}

### overwrite cosso::bigGram: generate coordinate-wise kernel matrices
bigGram <- function (x1, x2 = x1) {
  if(is.vector(x1)) {
    x1 <- matrix(x1, dimnames = list(names(x1)))
  } else if(!is.matrix(x1) & !is.data.frame(x1)) {
    stop("x1 should be a vector, a matrix or a data.frame")
  }
  if(is.vector(x2)) {
    x2 <- matrix(x2, dimnames = list(names(x2)))
  } else if(!is.matrix(x2) & !is.data.frame(x2)) {
    stop("x2 should be a vector, a matrix or a data.frame")
  }
  
  # format rows
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  if(is.null(rownames(x1))) {
    rownames(x1) <- 1:n1
  }
  if(is.null(rownames(x2))) {
    rownames(x2) <- 1:n2
  }
  
  # format columns
  if(ncol(x1) != ncol(x2)) {
    stop("columns in x1 and x2 mismatch")
  } else {
    d <- ncol(x2)
    if(is.null(colnames(x2))) {
      colnames(x2) <- paste0("X", 1:d)
    }
  }
  
  Gram <- array(dim = c(n1, n2, d), dimnames = list(rownames(x1), rownames(x2), colnames(x2)))
  levels <- replicate(d, NULL)
  range <- replicate(d, NULL)
  names(levels) <- names(range) <- colnames(x2)
  for (j in 1:d) {
    x1_j <- x1[, j]
    x2_j <- x2[, j]
    
    if(is.factor(x2_j)) {
      levels[[j]] <- levels_j <- levels(x2_j)
      x1_j <- factor(x1_j, levels = levels_j)
      if(any(is.na(x1_j))) {
        stop(paste0("x1[, ", j, "] has unknown levels"))
      }
    } else if(length(levels_j <- sort(unique(c(x1_j, x2_j)))) <= 6) {
      levels[[j]] <- levels_j
      x1_j <- factor(x1_j, levels = levels_j)
      x2_j <- factor(x2_j, levels = levels_j)
    } else {
      levels_j <- NULL
      if(min(c(x1_j, x2_j)) < 0 | max(c(x1_j, x2_j)) > 1) {
        message("rescale x1, x2 to [0, 1]")
        range_j <- range(c(x1_j, x2_j))
        names(range_j) <- c("lower", "upper")
        range[[j]] <- range_j
        
        x1_j <- (x1_j - range_j["lower"])/diff(range_j)
        x2_j <- (x2_j - range_j["lower"])/diff(range_j)
      }
    }
    
    if(is.null(levels_j)) { # x1_j, x2_j are numeric
      Gram[, , j] <- genK(x1_j, x2_j)
    } else { # x1_j, x2_j are categorical
      Gram[, , j] <- genK.cat(x1_j, x2_j, length(levels_j))
    }
  }
  attributes(Gram) <- c(attributes(Gram), list(levels = levels, range = range))
  return(Gram)
}
####################################################################################################
##### UNIT TEST for "bigGram"
# X <- with(list(n = 10, p = 4, K = 3, p_cat = 2), {
#   X_cat <- replicate(p_cat, sample(K, n, replace = TRUE))
#   p_num <- p - p_cat
#   X_num <- matrix(rnorm(n*p_num), nrow = n, ncol = p_num)
#   return(cbind(X_cat, X_num))
# })
# bigGram(X, X)

# library(dplyr)
# X <- data.frame(X) %>%
#   mutate(across("X1", ~ factor(., levels = 1:3)))
# bigGram(X, X)
# bigGram(matrix(1, nrow = 4, ncol = 4), X)
# bigGram(matrix(c(rep(1, 4), 5:8, rnorm(8)), nrow = 4, ncol = 4), X)
# bigGram(matrix(1:16, nrow = 4, ncol = 4), X)
####################################################################################################

### overwrite cosso::rescale: rescale a data.frame by column with a list of specified ranges
rescale <- function(X, range) {
  for(var in names(range)) {
    if(is.null(X[[var]])) {
      stop("variable names in X and range mismatch")
    }
    # check legitimacy of range[[var]]
    range_var <- range[[var]]
    is_range <- FALSE
    if(is.numeric(range_var) & length(range_var) >= 2) if(range_var[1] <= range_var[2]) {
      # range_x is legitimate  
      is_range <- TRUE
      range_var <- range_var[1:2]
    }
    if(!is_range) {
      stop(paste0("range[[", var, "]] is illegal"))
    }
    X[[var]] <- pmax(0, pmin(1, (X[[var]] - range_var[1])/diff(range_var)))
  }
  return(X)
}

# generate kernel matrices for the SS-ANOVA model: ~ 1 + X1 + X2 + X1 * X2
kernel_matrix_interaction <- function(X1, X2, id_basis = rownames(X2)) {
  # id_basis: sample representers in X1, X2, specified as a character subset of rownames(X1)
  
  # format rows
  n <- nrow(X1)
  if(nrow(X2) != n) {
    stop("sample sizes in X1 and X2 mismatch")
  }
  if(!all(id_basis %in% rownames(X1))) {
    stop("id_basis is not a subset of rownames(X1)")
  }
  if(is.null(rownames(X1))) {
    rownames(X1) <- id_basis <- as.character(1:nrow(X1))
  }
  n_basis <- length(id_basis)
  rownames(X2) <- rownames(X1)
  
  # format columns
  d1 <- ncol(X1)
  d2 <- ncol(X2)
  if(is.null(colnames(X1))) {
    colnames(X1) <- paste0("X1.", 1:d1)
  }
  if(is.null(colnames(X2))) {
    colnames(X2) <- paste0("X2.", 1:d2)
  }
  
  # format output Gram matrix
  d <- d1 + d2 + d1 * d2
  names_X1 <- colnames(X1)
  names_X2 <- colnames(X2)
  names_X12 <- paste(names_X1[rep(1:d1, times = d2)], names_X2[rep(1:d2, each = d1)], sep = ":")
  names_Gram <- c(names_X1, names_X2, names_X12)
  Gram <- array(dim = c(n, n_basis, d), dimnames = list(rownames(X1), id_basis, names_Gram))
  
  K1 <- bigGram(X1, X1[id_basis, , drop = FALSE])
  K2 <- bigGram(X2, X2[id_basis, , drop = FALSE])
  K12 <- K1[, , rep(1:d1, times = d2), drop = FALSE] * K2[, , rep(1:d2, each = d1), drop = FALSE]
  Gram[] <- abind(K1, K2, K12, along = 3)
  
  attr_K1 <- attributes(K1)[c("levels", "range")]
  names(attr_K1) <- paste0(names(attr_K1), "_X1")
  attr_K2 <- attributes(K2)[c("levels", "range")]
  names(attr_K2) <- paste0(names(attr_K2), "_X2")
  attributes(Gram) <- c(attributes(Gram), attr_K1, attr_K2)
  
  return(Gram)
}
####################################################################################################
##### UNIT TEST for "kernel_matrix_interaction"
# data <- with(list(n = 10, K = 3, p_cat = 2, p_num = 2), {
#   X_cat <- replicate(p_cat, sample(K, n, replace = TRUE))
#   X_num <- matrix(runif(n*p_num), nrow = n, ncol = p_num)
#   return(list(X1 = data.frame(X_cat), X2 = data.frame(X_num)))
# })
# with(data, kernel_matrix_interaction(X1, X2))
# with(data, kernel_matrix_interaction(X1, X2, 2*(1:5)))
# 
# # large n, p
# data <- with(list(n = 400, K = 3, p_cat = 1, p_num = 100), {
#   X_cat <- replicate(p_cat, sample(K, n, replace = TRUE))
#   X_num <- matrix(runif(n*p_num), nrow = n, ncol = p_num)
#   return(list(X1 = data.frame(X_cat), X2 = data.frame(X_num)))
# })
# time1 <- system.time({
#   K1 <- with(data, kernel_matrix_interaction(X1, X2))
# })
# time2 <- system.time({
#   id_basis <- with(data, {
#     n_basis <- max(40, ceiling(12 * n^(2/9)))
#     return(sort(sample(n, n_basis)))
#   })
#   
#   K2 <- with(data, kernel_matrix_interaction(X1, X2, id_basis))
# })
# time1; time2; format(object.size(K1), units = "Mb"); format(object.size(K2), units = "Mb")
####################################################################################################

# COSSO of the SS-ANOVA model: Y ~ 1 + X1 + X2 + X1 * X2, family = "Gaussian"
### define internal functions
cosso_interaction_internal <- within(list(), {
  ### define data preparation function
  cosso_interaction.data <- function(data = NULL, data_representer = NULL, config = CONFIG) {
    # if is.null(data_representer), return data_train; else, return data_predict
    
    assign_list(extract_list(data, c("X1", "X2", "Y")))
    assign_list(extract_list(config, c("levels_X1", "range_X1", "levels_X2", "range_X2")))
    
    # format X1, X2
    format_X1 <- format_X(X1, rename_list(config, c("range" = "range_X1", "levels" = "levels_X1")))
    format_X2 <- format_X(X2, rename_list(config, c("range" = "range_X2", "levels" = "levels_X2")))
    data <- list(X1 = format_X1[[c("data", "X")]], X2 = format_X2[[c("data", "X")]], Y = Y)
    config <- c(
      with(format_X1, {
        c(with(data, list(d1 = ncol(X), names_X1 = colnames(X))),
          rename_list(config, c("levels_X1" = "levels", "range_X1" = "range")))
      }),
      with(format_X2, {
        c(with(data, list(d2 = ncol(X), names_X2 = colnames(X))),
          rename_list(config, c("levels_X2" = "levels", "range_X2" = "range")))
      })
    )
    
    # rescale numerical variables in X1, X2 to [0, 1]
    X1 <- with(format_X1, rescale(data[["X"]], config[["range"]]))
    X2 <- with(format_X2, rescale(data[["X"]], config[["range"]]))
    n <- nrow(X1)
    if(nrow(X2) != n) {
      stop("sample sizes in X1 and X2 mismatch")
    } else {
      rownames(X2) <- rownames(X1)
    }
    
    if(is.null(data_representer)) {
      is_train <- TRUE
      # sample representers
      if(n < 40) {
        warning("sample size < 40")
        n_basis <- n
      } else {
        n_basis <- max(40, ceiling(12 * n^(2/9)))
      }
      id_basis <- rownames(X1)[sort(sample(n, n_basis))]
      data_representer <- list(
        X1_basis = X1[id_basis, , drop = FALSE],
        X2_basis = X2[id_basis, , drop = FALSE]
      )
    } else { # !is.null(data_representer)
      is_train <- FALSE
      n_basis <- nrow(data_representer[["X1_basis"]])
      id_basis <- rownames(data_representer[["X1_basis"]])
    }
    config <- c(config, extract_list(as.list(environment()), c("n", "n_basis", "id_basis")))
    
    if(is_train) {
      # generate training kernel matrices
      kernel_matrix <- kernel_matrix_interaction(X1, X2, id_basis)
      data_train <- list(Kmat = kernel_matrix, y = Y)
      config <- c(config, list(d = dim(kernel_matrix)[3]))
    } else { # !is_train
      data_predict <- with(data_representer, within(list(), {
        K1 <- bigGram(X1, X1_basis)
        K2 <- bigGram(X2, X2_basis)
      }))
    }
    return(extract_list(
      as.list(environment()), 
      c("data", "data_train", "data_representer", "data_predict", "config")
    ))
  }
})
### model fitting function
cosso_interaction <- with(cosso_interaction_internal, {
  function(data, config = CONFIG, ...) {
    # ...: unused argument
    
    # format input
    assign_list(cosso_interaction.data(data, , config), 
                c("data", "data_train", "data_representer", "config"))
    
    # fit COSSO
    training_time <- system.time({
      model <- with(c(data_train, config), within(list(), {
        Kmat <- Kmat
        y <- y
        wt <- rep(1, d)
        assign_list(cosso.Gaussian(Kmat, y, wt, n_basis, id_basis))
      }))
      # extract coefficients = list(intercept, representer coefficients, kernel weights)
      options(warn = -1)
      coefs <- predict.cosso(model, type = "coefficients"); class(coefs) <- "list"
      options(warn = 0)
    })
    
    # format coefficients
    coefficients <- with(c(coefs, config), within(list(), {
      # intercept
      intercept <- intercept
      # representer coefficients
      names(coef) <- id_basis
      # kernel weights
      theta1 <- theta[1:d1]
      names(theta1) <- names_X1
      theta2 <- theta[d1+(1:d2)]
      names(theta2) <- names_X2
      theta12 <- matrix(theta[-(1:(d1+d2))], nrow = d1, ncol = d2,
                        dimnames = list(names_X1, names_X2))
    }))
    
    output <- extract_list(
      as.list(environment()), 
      c("data_representer", "coefficients", "training_time", "data", "config")
    )
    # in-sample prediction
    prediction <- predict.cosso_interaction(output, data)
    return(c(output, list(prediction = prediction)))
  }
})
### prediction function
predict.cosso_interaction <- with(cosso_interaction_internal, {
  function(object, newdata = NULL, ...) {
    # object = a list of: data_representer, coefficients, config[, prediction, data]
    # return = an n-by-K-matrix prediction
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
    assign_list(with(object, cosso_interaction.data(newdata, data_representer, config)),
                c("data_predict", "config"))
    assign_list(object[["coefficients"]])
    assign_list(data_predict)
    assign_list(config, c("n", "n_basis"))
    
    return(intercept + foreach(id_basis = 1:n_basis, .combine = "+") %do% {
      coef[id_basis] * sapply(1:n, function(id) {
        k1 <- K1[id, id_basis, ]
        k2 <- K2[id, id_basis, ]
        theta1 %*% k1 + theta2 %*% k2 + k1 %*% theta12 %*% k2
      })
    })
  }
})
####################################################################################################
##### UNIT TEST for "cosso_interaction"
# config <- within(list(), {
#   n <- 400
#   p1_cat <- 1
#   p1_num <- 5
#   K1 <- 3
#   p1 <- p1_cat + p1_num
#   p2_cat <- 1
#   p2_num <- 5
#   K2 <- 3
#   p2 <- p2_cat + p2_num
# 
#   levels_X1 <- replicate(p1_cat, 1:K1, simplify = FALSE)
#   names(levels_X1) <- paste0("X", 1:p1_cat)
#   levels_X2 <- replicate(p2_cat, 1:K2, simplify = FALSE)
#   names(levels_X2) <- paste0("X", 1:p2_cat)
#   range_X1 <- replicate(p1_num, c(-3, 3), simplify = FALSE)
#   names(range_X1) <- paste0("X", p1_cat+(1:p1_num))
#   range_X2 <- replicate(p2_num, c(-3, 3), simplify = FALSE)
#   names(range_X2) <- paste0("X", p2_cat+(1:p2_num))
# })
# data_interaction <- function(config) {
#   data <- with(config, within(list(), {
#     X1_cat <- matrix(sample(K1, n*p1_cat, replace = TRUE), nrow = n, ncol = p1_cat)
#     X1_num <- matrix(rnorm(n*p1_num), nrow = n, ncol = p1_num)
#     X1 <- cbind(X1_cat, X1_num)
#     X2_cat <- matrix(sample(K2, n*p2_cat, replace = TRUE), nrow = n, ncol = p2_cat)
#     X2_num <- matrix(rnorm(n*p2_num), nrow = n, ncol = p2_num)
#     X2 <- cbind(X2_cat, X2_num)
# 
#     mean_Y <- rowMeans(X1) + row_means(X2) +
#       extract(X1_num, X2_cat[, 1]) + extract(X2_num, X1_cat[, 1])
#     Y <- mean_Y + rnorm(n)
#   }))
# 
#   return(c(data, config = list(config)))
# }
# library(dplyr)
# # training
# data <- data_interaction(config)
# time_cosso <- system.time({
#   model_cosso <- with(data, cosso(cbind(X1, X2), Y))
# })
# time_cosso_interaction <- system.time({
#   model_cosso_interaction <- cosso_interaction(data, config)
# })
# # in-sample prediction
# within(data, prediction <- predict.cosso(model_cosso, cbind(X1, X2))) %>%
#   with(qqplot(mean_Y, prediction)); abline(0, 1)
# with(c(data, model_cosso_interaction), qqplot(mean_Y, prediction)); abline(0, 1)
# # testing prediction
# data <- data_interaction(within(config, { n <- 10000 }))
# within(data, prediction <- predict.cosso(model_cosso, cbind(X1, X2))) %>%
#   with(qqplot(mean_Y, prediction)); abline(0, 1)
# within(data, prediction <- predict.cosso_interaction(model_cosso_interaction, data)) %>%
#   with(qqplot(mean_Y, prediction)); abline(0, 1)
####################################################################################################
