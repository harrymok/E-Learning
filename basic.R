# quick computation of sapply(1:nrow(matrix), function(i) matrix[i,id[i]])
extract <- function(matrix, id) {
  n <- nrow(matrix)
  m <- ncol(matrix)
  if(any(!id %in% 1:m)) {
    stop("column id overflows")
  } else {
    id <- as.numeric(id)
  }
  return(t(matrix)[id + m*(1:n - 1)])
}
####################################################################################################
##### UNIT TEST for "extract"
# X <- matrix(1:9, nrow = 3, ncol = 3)
# id <- c(2, 3, 1)
# extract(X, id)
# sapply(1:nrow(X), function(i) X[i,id[i]])
####################################################################################################

# generate angle-based codes
code_generation <- function(K) {
  rbind((K-1)^{-1/2} * rep(1,K-1), - (1+sqrt(K))/(K-1)^{3/2} + (K/(K-1))^{1/2} * diag(1,K-1))
}
####################################################################################################
##### UNIT TEST for "code_generation"
# code_generation(2)  # c(1, -1)
# K <- 4
# code_generation(K) %*% t(code_generation(K))  # diagonal = 1, off-diagonal = -1/(K-1)
# t(code_generation(K)) %*% code_generation(K)  # (K/(K-1)) * identity
####################################################################################################

# row-wise dot product
`%*_row%` <- function(A, B) {
  rowSums(A * B)
}
# row-wise Kronecker product
`%x_row%` <- function(A, B) {
  as.matrix(t(Matrix::KhatriRao(t(A), t(B))))
}
####################################################################################################
##### UNIT TEST for `%x_row%`
# diag(3) %x_row% matrix(1:9, nrow = 3, ncol = 3)  # diagonal blocks of 1:3, 4:6, 7:9
####################################################################################################

# expit
expit <- function(u) {
  exp(u) / (1 + exp(u))
}
####################################################################################################
##### UNIT TEST for expit
# plot(expit, xlim = c(-5, 5))
####################################################################################################

# extract list elements according to a character vector of element names
extract_list <- function(List = list(), names_extract = character()) {
  names_keep <- intersect(names_extract, names(List))
  names_NULL <- setdiff(names_extract, names_keep)
  List_keep <- List[names_keep]
  List_NULL <- replicate(length(names_NULL), NULL)
  names(List_NULL) <- names_NULL
  return(c(List_keep, List_NULL))
}
####################################################################################################
##### UNIT TEST for "extract_list"
# extract_list(list(x=1), c("x", "y"))
####################################################################################################

# set list elements according to another list
left_full_join_list <- function(List_left = list(), List_right = list()) {
  # return a list including all elements from List_left, and elements from List_right without a match in List_Left 
  names_set <- setdiff(names(List_right), names(List_left))
  List_left[names_set] <- List_right[names_set]
  return(List_left)
}
####################################################################################################
##### UNIT TEST for "extract_list"
# left_full_join_list(list(x = 1, y = 1), list(x = 2, z = 2))
# left_full_join_list(NULL, list(x = 1))
####################################################################################################

# rename list elements according to a named character vector: new_name = original_name
rename_list <- function(List = list(), renames = character()) {
  for(id_renames in 1:length(renames)) {
    name <- renames[[id_renames]]
    id_List <- grep(paste0("^", name, "$"), names(List))
    names(List)[id_List] <- names(renames)[id_renames]
  }
  return(List)
}
####################################################################################################
##### UNIT TEST for "rename_list"
# rename_list(list(x = 1, y = 2), c("z" = "x", "zz" = "xx"))
####################################################################################################

# assign list elements to the parent frame
assign_list <- function(List = list(), names_assign = character()) {
  # names_assign: a character vector of names to assign; 
  ### if NULL, overwritten by List[["names_assign"]]
  ### if List[["names_assign"]] == NULL, assign all elements in List
  if(length(names_assign) == 0) {
    if(is.null(List[["names_assign"]])) {
      names_assign <- names(List)
    } else {
      names_assign <- List[["names_assign"]]
    }
  }
  for(name in names_assign) {
    assign(name, List[[name]], pos = parent.frame())
  }
}
####################################################################################################
##### UNIT TEST for "assign_list"
# rm(a)  # "a" does not exist in the current environment
# assign_list(list(a = 1)); a  # "a" is assigned as 1
# 
# rm(a)  # "a" does not exist in the current environment
# assign_list(list(a = 1, names_assign = "a")); a  # "a" is assigned as 1
# 
# rm(a, b)  # "a", "b" do not exist in the current environment
# assign_list(list(a = 1), names_assign = "b")
# a  # not assigned
# b  # NULL
# 
# rm(a)  # "a" does not exist in the current environment
# fun <- function(x) { assign_list(list(a = 1, names_assign = "a")); return(x + a) }
# fun(1)  # "a" is assigned as 1 within "fun"
# a  # "a" does not exist in the current environment
####################################################################################################

# format inputs
### covariates
format_X <- function(X, config = CONFIG) {
  config <- extract_list(config, c("levels", "range"))
  assign_list(config)
  
  # check eligibility
  if(is.vector(X) | is.factor(X)) {
    X <- data.frame(X = X)
  } else if(!is.matrix(X) & !is.data.frame(X)) {
    stop("X should be a vector, a factor, a matrix or a data.frame")
  }
  
  # check intercept
  is_intercept <- FALSE
  if(is.numeric(X[, 1])) if(all(abs(X[, 1] - 1) <= 1e-10)) {
    # the first column is an intercept
    X <- X[, -1]
  }
  
  # convert to a data.frame
  if(is.matrix(X)) {
    X <- data.frame(X)
  } 
  
  # check names
  vars <- colnames(X)  # all variables
  vars_cat <- names(levels)  # categorical variables
  vars_num <- names(range)   # numerical variables
  if(length(intersect(vars_cat, vars_num)) > 0) {
    stop("levels and range cannot be simultaneously specified")
  }
  if(any(sort(vars) != sort(c(vars_cat, vars_num)))) { # note: if vars_cat, vars_num = NULL, then it returns FALSE
    stop("variables in X and (levels, range) mismatch")
  }
  
  # identify categorical and numerical variables in X
  ### categorical variable: converted to a factor and update levels_X
  ### numerical variable: update range_X
  for(var in vars) {
    x <- X[[var]]
    if(!is.null(levels_x <- levels[[var]])) {
      # x is a specified as a factor
      range[[var]] <- NULL
      
      # convert X[[var]] to a factor with the given levels
      X[[var]] <- x <- factor(x, levels = levels_x)
      if(any(is.null(x))) {
        stop(paste0(var, " has unknown levels"))
      }
    } else if(is.factor(x)) {
      # x is a factor
      levels[[var]] <- levels(x)
      range[[var]] <- NULL
    } else if(is.null(range[[var]]) & length(levels_x <- sort(unique(x))) <= 6) {
      # x has <= 6 unique values, treated as categorical
      levels[[var]] <- levels_x
      
      # convert X[[var]] to a factor
      X[[var]] <- factor(x, levels = levels_x)
    } else { 
      # x is numeric
      range_x <- range[[var]]
      ### check and overwrite range[[var]]
      is_range <- FALSE
      if(is.numeric(range_x) & length(range_x) >= 2) if(range_x[1] <= range_x[2]) {
        # range_x is legitimate
        is_range <- TRUE
        range_x <- range_x[1:2]
      }
      if(!is_range) { 
        message(paste0("overwrite the range of ", var))
        range_x <- range(x)
      }
      names(range_x) <- c("lower", "upper")
      range[[var]] <- range_x
    }
  }
  
  # construct the design matrix, including the intercept
  design <- model.matrix(~ ., data.frame(X))
  
  data <- list(X = X, design = design)
  config <- list(levels = levels, range = range, groups_design = attr(design, "assign"))
  return(list(data = data, config = config))
}
####################################################################################################
##### UNIT TEST for "format_X"
# format_X(1:3)
# format_X(factor(1:3))
# X <- matrix(1:14, ncol = 2, nrow = 7)
# format_X(X)
# format_X(cbind(1,X))
# format_X(X, levels = list(X1 = 1:7))
# format_X(X, levels = list(X1 = 1:7), range = list(X2 = NULL))
# colnames(X) <- paste0("x", 1:2)
# format_X(X)
# X <- data.frame(X)
# X[, 1] <- factor(X[, 1])
# format_X(X)
####################################################################################################

##### treatment
format_A <- function(A = NULL, config = CONFIG) {
  config <- extract_list(config, "levels")
  assign_list(config)
  
  if(!is.null(levels)) {
    levels <- levels
  } else if(is.null(A)) {
    stop("either A or levels should be specified")
  } else if(is.factor(A)) {
    levels <- levels(A)
  } else {
    levels <- sort(unique(A))
  }
  K <- length(levels)
  if(K == 1) {
    levels <- c(levels, NA)
    K <- 2
  }
  
  # generate code matrix Omega
  Omega <- code_generation(K)
  rownames(Omega) <- levels
  colnames(Omega) <- 1:(K-1)
  
  data <- within(list(), {
    if(!is.null(A)) {
      A <- factor(A, levels = levels)
      id <- as.numeric(A)
      design <- Omega[id, , drop = FALSE]
    }
  })
  config <- list(levels = levels, K = K, Omega = Omega)
  return(list(data = data, config = config))
}
####################################################################################################
##### UNIT TEST for "format_A"
# format_A(rep(c(2, 1, 3, 4), each = 3))
# format_A(rep(c("B", "A", "C", "D"), each = 3))
# format_A(config = list(levels = 1:4))
# format_A(1:3, config = list(levels = 1:4))
####################################################################################################

##### format data
format_data <- function(data, config = CONFIG) {
  # keep copies
  data_input <- data
  config_input <- config
  
  # format_input
  assign_list(extract_list(data, c("X", "A", "prop")))
  config <- extract_list(config, c("levels_X", "range_X", "levels_A"))
  
  # format covariates
  format_X <- format_X(X, rename_list(config, c("range" = "range_X", "levels" = "levels_X")))
  data_X <- rename_list(format_X[["data"]], c("design_X" = "design"))
  config_X <- with(format_X, {
    c(with(data, list(n = nrow(X), p = ncol(design), names_X = colnames(design))),
      rename_list(config, c("range_X" = "range", "levels_X" = "levels", "groups_X" = "groups_design")))
  })
  
  # format treatments
  format_A <- format_A(A, rename_list(config, c("levels" = "levels_A")))
  data_A <- rename_list(format_A[["data"]], c("A_id" = "id", "design_A" = "design"))
  config_A <- with(format_A, rename_list(config, c("levels_A" = "levels")))
  if(!is.null(data_A[["A_id"]]) & !is.null(prop)) {
    data_A <- within(data_A, {
      colnames(prop) <- config_A[["levels_A"]]
      prop_A <- extract(prop, A_id)
    })
  }
  
  data <- left_full_join_list(c(data_X, data_A), data_input)
  config <- left_full_join_list(c(config_X, config_A), config_input)
  return(list(data = data, config = config))
}
####################################################################################################
##### UNIT TEST for "format_data"
# data <- within(list(), {
#   X <- cbind(rep(1:3, 4),
#              matrix(1:9, ncol = 3, nrow = 12))
#   A <- rep(c(2, 1, 3, 4), each = 3)
#   prop <- matrix(1/4, 12, 4)
#   Y <- rnorm(12)
# })
# config <- within(list(), {
#   range_X <- replicate(3, c(-3, 3), simplify = FALSE)
#   names(range_X) <- c("X2", "X3", "X4")
#   levels_X <- list(X1 = 1:4)
#   levels_A <- 1:4
# })
# data <- format_data(data, config)
####################################################################################################