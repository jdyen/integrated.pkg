# extract a node from a greta array
get_node <- greta:::get_node

# convert node to greta array
as.greta_array <- greta:::as.greta_array

# extract unique values without keeping an extra data node around
extract_unique <- function(params) {
  
  node <- get_node(params)
  
  unique(node$value())
  
}

# change dimensions of greta distributions
change_dims <- function(array, new_dim) {
  
  # I want the good bits from the input greta array
  node <- get_node(array)
  
  # is this an operation or distribution node?
  if ("distribution_node" %in% class(node$distribution)) {
    
    new_array <- change_dims_distrib(node, new_dim)
    
  } else {
  
    # this has to be an operation node
    if (!("operation_node" %in% class(node)))
      stop("can only change dimensions of probability distributions", call. = FALSE)
    
    new_array <- change_dims_op(node, new_dim)

  }
  
  # return a new array with new dimensions
  new_array
  
}

change_dims_distrib <- function(node, new_dim) {
  
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
  params <- lapply(node$distribution$parameters, function(x) as.greta_array(x))
  
  # is it a multivariate distribution
  is_multivariate <- distribution_name %in% c("multivariate_normal", "dirichlet")
  
  # what are the dims of the input?
  dims <- node$dim
  
  # input distribution can be multidimensional in some cases
  if (is_multivariate) {
    
    # is it multidimensional?
    if (dims[1] > 1) {
      
      # are all params equal?
      if (!identical(dim(params[[1]]), dims))
        stop(paste0("cannot set multivariate priors with dimensions ", paste(dims, collapse = " x "), "unless all parameters are equal"), call. = FALSE)
      
      # the covariance has to be a single matrix for this to work
      if (distribution_name == "multivariate_normal") {
        if (!identical(dim(params[[2]]), rep(dims[2], 2)))
          stop("multivariate_normal prior covariance must match dimension of distribution", call. = FALSE)
      }
      
    }
    
  } else {
    
    # can be multidimensional univariate if all params are identical
    if (!all(sapply(params, function(x) length(unique(x))) == 1))
      stop(paste0("cannot set priors with dimensions ", paste(dims, collapse = " x "), " unless all parameters are equal"), call. = FALSE)
    
    # need to collapse parameter vectors down if all equal
    idx <- which(sapply(params, length) > 1)
    if (length(idx)) {
      params[idx] <- lapply(params[idx], extract_unique)
    }
    
  }
  
  # if it's uniform, we need to turn params back to numerics
  if (distribution_name == "uniform")
    params <- lapply(params, function(x) as.numeric(x))
  
  # what are the arguments required for this distribution?
  arglist <- names(formals(node$initialize))
  
  # we need to add a new dim argument to univariate parameters
  if (!is_multivariate) {
    params$dim <- new_dim
  } else {
    # we need to add a new n_realisations argument to multivariate parameters; this must be scalar
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
  new_distrib <- do.call(distribution_name, params)
  
  # we need the new node so we can cut it off from the input distribution
  new_node <- get_node(new_distrib)
  
  # remove original variable greta array
  for (i in seq_along(node$distribution$parameters))
    node$distribution$parameters[[i]]$remove_parent(new_node)

  # return a new distribution with new dims
  new_distrib
  
}

# change dims of a greta operation node
change_dims_op <- function(node, new_dim) {

  # what is the operation?
  op <- node$operation_name
  
  # which operations are OK?
  valid_ops <- c("ilogit", "iprobit", "icloglog", "icauchit",
                 "log1pe", "imultilogit",
                 "add", "subtract", "multiply", "divide", "power",
                 "log", "exp", "log1p", "expm1",
                 "cos", "sin", "tan", "acos", "asin", "atan")
  
  # is this `op` OK?
  if (!(op %in% valid_ops)) 
    stop(paste0("unable to change dimensions for operation: ", op), call. = FALSE)
  
  # `add` doesn't exist as a greta operation
  if (op %in% c("add", "subtract", "multiply", "divide", "power")) {
    op <- switch(op,
                 "add" = "+",
                 "subtract" = "-",
                 "multiply" = "*",
                 "divide" = "/",
                 "power" = "^")
  }
  
  # we need the node for the underlying distribution
  distrib_node <- node$children
  
  # want to make sure the children are all distributions
  if (any(sapply(distrib_node, function(x) is.null(x$distribution))))
    stop(paste0("one of the input nodes to ", op, " is not a probability distribution"), call. = FALSE)
  if (!all(sapply(distrib_node, function(x) "distribution_node" %in% class(x$distribution))))
    stop(paste0("one of the input nodes to ", op, " is not a probability distribution"), call. = FALSE)
  
  # change dims of underlying distributions
  new_children <- lapply(distrib_node, change_dims_distrib, new_dim = new_dim)
  
  # pull out extra arguments to op
  op_args <- node$operation_args
  
  # recombine children nodes and apply operation
  if (length(op_args)) {
    new_array <- do.call(op, new_children, op_args)
  } else {
    new_array <- do.call(op, new_children)
  }
  
  # return a new array with new dims
  new_array

}

# extract lower and upper bounds of a greta distribution
extract_bounds <- function(x) {
  
  x_node <- get_node(x)
  
  c(lower = x_node$lower, upper = x_node$upper)
  
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
