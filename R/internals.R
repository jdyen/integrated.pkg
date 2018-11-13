# internal helper functions for integrated R package

# count stages that an individual successfully transitions out of
count_stages_survived <- function(x, classes) {
  out <- rep(0, classes)
  xpos <- x[x > 0]
  xsurv <- xpos[xpos < max(xpos)]
  out[seq_len(classes) %in% xsurv] <- 1
  tmp_var <- runif(classes)
  out[tmp_var < seq(0.3, 1, length = classes)] <-
    ifelse(seq_len(classes) %in% xpos, 1, 0)[tmp_var < seq(0.3, 1, length = classes)]
  out
}

# count stages that an individual is observed in
count_stages_lived <- function(x, classes) {
  out <- rep(0, classes)
  x <- x[x > 0]
  out[seq_len(classes) %in% x] <- 1
  out
}

# calculate stage at last observation
calculate_final_stage <- function(x) {
  x[max(which(x != 0))]
}

# set an object class
as_class <- function (object, name, type = c("function", "list")) {
  
  type <- match.arg(type)
  stopifnot(inherits(object, type))
  class(object) <- c(name, class(object))
  
  object
  
}

# test if values are equal to zero with a non-zero tolerance
almost_equal <- function(x, y, tolerance = 1e-10) {
  
  diff <- abs(x - y)
  mag <- pmax(abs(x), abs(y))
  
  ifelse(mag > tolerance, (diff / mag) <= tolerance, diff <= tolerance)
  
}

# test if multiple values are all equal
all_equal <- function (..., tolerance = 1e-10) {
  
  test <- list(...)
  
  if (length(test) > 1) {
    out <- rep(NA, length(test))
    for (i in seq_along(test)) {
      out[i] <- almost_equal(test[[i]], test[[1]], tolerance = tolerance)
    }
  } else {
    out <- TRUE
  }
  
  all(out)
  
}

# check parameters for dims and class
check_params <- function(params, classes) {

  class_mismatch <- FALSE
  length_mismatch <- FALSE
  vec_mismatch <- FALSE
  dim_mismatch <- FALSE
  type <- sapply(params, class)
  lengths <- sapply(params, length)
  
  # all can be scalar
  if (!all(type %in% c('numeric', 'matrix')))
    class_mismatch <- TRUE
  
  # fec, dens, and capture params must be scalar
  if (any(lengths[c(1:2, 7:10)] > 1))
    length_mismatch <- TRUE
  
  # surv params can be classes x 1
  if (lengths[3] > 1) {
    if (lengths[3] != classes)
      vec_mismatch <- TRUE
  }
  if (lengths[4] > 1) {
    if (lengths[4] != classes)
      vec_mismatch <- TRUE
  }
  
  # mat params can also be classes x classes mat
  if (lengths[5] > 1) {
    if (type[5] != 'matrix') {
      class_mismatch <- TRUE
    } else {
      if (nrow(params[[5]]) != classes)
        dim_mismatch <- TRUE
    }
  }
  if (lengths[6] > 1) {
    if (type[6] != 'matrix') {
      class_mismatch <- TRUE
    } else {
      if (nrow(params[[6]]) != classes)
        dim_mismatch <- TRUE
    }
  }
  
  list(length_error = length_mismatch,
       class_error = class_mismatch,
       vec_error = vec_mismatch,
       dim_error = dim_mismatch)
  
}

get_node <- greta:::get_node

# expand dimensions of prior distributions to match data
expand_prior <- function(prior, new_dim) {
  
  # I want the good bits from the prior
  node <- get_node(prior)
  
  # which distribution am I dealing with?
  distribution_name <- node$distribution$distribution_name
  
  # which distributions are not OK?
  discrete <- c("bernoulli", "binomial", "negative_binomial", "poisson",
                "multinomial", "categorical", "dirichlet_multinomial")
  too_hard <- "wishart_distribution"
  
  # is this distribution ok?
  if (distribution_name %in% discrete)
    stop(paste0("unable to set discrete (", distribution_name, ") distribution as a prior"), call. = FALSE)
  if (distribution_name %in% too_hard)
    stop(paste0("unable to set ", distribution_name, " distribution as a prior"), call. = FALSE)
  
  # is the distribution truncated?
  lower <- node$lower
  upper <- node$upper
  
  # what are its parameters?
  params <- lapply(node$distribution$parameters, function(x) greta:::as.greta_array(x))
  
  # if it's uniform, we need to turn params back to numerics
  if (distribution_name == "uniform")
    params <- lapply(params, function(x) as.numeric(x))
  
  # what are the arguments required for this distribution?
  arglist <- names(formals(node$initialize))
  
  # we need to add a new dim argument to univariate parameters
  if (!(distribution_name %in% c("multivariate_normal", "dirichlet"))) {
    params$dim <- new_dim
  } else {
    if (length(new_dim) > 1)
      stop("unable to expand a multivariate distribution in two dimensions", call. = FALSE)
    params$n_realisations <- new_dim
  }
  
  # the f distribution's name is "d" for some reason
  if (distribution_name == "d")
    distribution_name <- "f"
  
  # do we also need to truncate the distribution?
  if ("truncation" %in% arglist)
    params$truncation <- c(lower, upper)
  
  # now we need a new distribution with these parameters
  new_prior <- do.call(distribution_name, params)
  
  # we need the new node so we can cut it off from the old prior distribution
  new_node <- get_node(new_prior)
  
  # remove original variable greta array
  for (i in seq_along(node$distribution$parameters))
    node$distribution$parameters[[i]]$remove_parent(new_node)

  # return a new prior with expanded dims
  new_prior
  
}
