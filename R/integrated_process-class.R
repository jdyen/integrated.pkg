#' @name integrated_process
#' @title create integrated process objects
#'
#' @description An \code{integrated_process} object contains the underlying
#'   process model for an integrated population analysis
#' 
#' @param classes something
#' @param density something
#' @param replicates number of distinct populations
#' @param params named list of parameters
#' @param ... additional arguments
#'
#' @details something
#'
#' @return An object of class \code{integrated_process}, which can be passed to
#'    \link[integrated]{integrated_data} and \link[integrated]{integrated_model}
#' 
#' @import greta
#' 
#' @examples
#' \dontrun {
#' 
#' library(integrated)
#' 
#' # prepare an example model
#' mpm_process <- leslie(5, density = "none", replicates = 3)
#' }


#' @export
#' @rdname integrated_process
#' 
leslie <- function(classes, density = "none", params = list()) {
  
  # set type
  type <- "leslie"
  
  # initialise model parameters
  ## WHAT'S A BETTER WAY TO SET LINKS?
  param_list <- list(density = uniform(0, 1),
                     density_link = "identity",
                     reproductive = max(classes),
                     survival = beta(1, 1),
                     survival_link = "identity",
                     fecundity = normal(0, 10, truncation = c(0, Inf)),
                     fecundity_link = "identity",
                     initials = normal(0, 10, truncation = c(0, Inf)),
                     random = normal(0, 10, truncation = c(0, Inf)))
  param_list[names(params)] <- params

  # do the link functions make sense?
  # e.g. check if gaussian with identity for probabs
  #   (important for surv/growth probabs)
  
  # update check_params function to test dists  
  test_params <- check_params(param_list, classes)
  # dims of each term should be 1 (expanded) or classes
  # all elements should be greta arrays with distributions

  # how many observations and predictor variables?
  nfec <- length(param_list$reproductive)

  # collate and return outputs  
  process <- list(type = type,
                  classes = classes,
                  density = density,
                  params = param_list)
  
  # return outputs
  as.integrated_process(process)

}

#' @export
#' @rdname integrated_process
#' 
lefkovitch <- function(classes, density = "none", params = list()) {
  
  # set type
  type <- "lefkovitch"
  
  # initialise model parameters
  param_list <- list(density = uniform(0, 1),
                     density_link = "identity",
                     reproductive = max(classes),
                     survival = beta(1, 1),
                     survival_link = "identity",
                     growth = beta(1, 1),
                     growth_link = "identity",
                     fecundity = normal(0, 10, truncation = c(0, Inf)),
                     fecundity_link = "identity",
                     initials = normal(0, 10, truncation = c(0, Inf)),
                     random = normal(0, 10, truncation = c(0, Inf)))
  param_list[names(params)] <- params
  
  # update check_params function to test dists  
  test_params <- check_params(param_list, classes)
  # dims of each term should be 1 (expanded) or classes
  # all elements should be greta arrays with distributions
  
  # how many observations and predictor variables?
  nfec <- length(param_list$reproductive)
  
  # collate and return outputs  
  process <- list(type = type,
                  classes = classes,
                  density = density,
                  params = param_list)
  
  # return outputs
  as.integrated_process(process)
  
}

#' @export
#' @rdname integrated_process
#' 
ipm <- function(classes, density = "none", params = list()) {
  
  # set type
  type <- "ipm"
  
  # will default params work for ipm setup?
  
  # should we warn that classes is now computational not process-based?
  
  # initialise model parameters
  param_list <- list(density = uniform(0, 1),
                     density_link = "identity",
                     reproductive = max(classes),
                     survival = beta(1, 1),
                     survival_link = "identity",
                     growth = beta(1, 1),
                     growth_link = "identity",
                     fecundity = normal(0, 10, truncation = c(0, Inf)),
                     fecundity_link = "identity",
                     initials = normal(0, 10, truncation = c(0, Inf)),
                     random = normal(0, 10, truncation = c(0, Inf)))
  param_list[names(params)] <- params
  
  # update check_params function to test dists  
  test_params <- check_params(param_list, classes)
  # dims of each term should be 1 (expanded) or classes
  # all elements should be greta arrays with distributions
  
  # how many observations and predictor variables?
  nfec <- length(param_list$reproductive)
  
  # collate and return outputs  
  process <- list(type = type,
                  classes = classes,
                  density = density,
                  params = param_list)
  
  # set class
  process <- as.integrated_process(process)
  
  # return outputs
  as.integrated_process(process)
  
}


integrated_process <- function(type, classes,
                               structure = 'stage',
                               density_dependence = 'none',
                               replicates = 1,
                               params = list()) {
  
  if (!(type %in% c('IPM', 'MPM'))) {
    stop('type must be one of IPM or MPM',
         call. = FALSE)
  }
  
  # account for repeated observations of multiple replicates
  if (length(replicates) > 1) {
    replicate_id <- replicates
    replicates <- length(unique(replicates))
  } else {
    if (replicates > 1) {
      replicate_id <- seq_len(replicates)
    } else {
      replicate_id <- NULL
    }
  }
  
  # fill params_list
  params_list <- list(fec_lower = 0,
                      fec_upper = 1000,
                      surv_param1 = 1,
                      surv_param2 = 1,
                      surv_mat_param1 = 1,
                      surv_mat_param2 = 1,
                      density_lower = 0.999,
                      density_upper = 1.001,
                      capture_lower = 0.999,
                      capture_upper = 1.0)
  params_list[names(params)] <- params
  
  # check dimensions and classes of parameters
  test_params <- check_params(params_list, classes)
  if (test_params$length_error)
    stop('fecundity, density, and capture probabilities should be scalar values', call. = FALSE)
  if (test_params$dim_error)
    stop(paste0('survival parameters should be scalar values or vectors with ', classes, ' elements'), call. = FALSE)
  if (test_params$vec_error)
    stop(paste0('survival matrix parameters should be scalar values or matrices with ', classes, ' rows and columns'), call. = FALSE)
  if (test_params$class_error)
    stop('all parameters should be scalar, vector, or matrix numeric objects', call. = FALSE)
  
  ### pass distrbns not params (include params in dist call)
  # prior_list <- list(fecundity = greta::uniform(min = 0, max = 1000, dim = 1),
  #                    survival = greta::beta(shape1 = 1, shape2 = 1, dim = 1)
  #                    class_survival = greta::beta(shape1 = 1, shape2 = 1, dim = 1),
  #                    density_dependence = greta::uniform(min = 0, max = 2, dim = 1),
  #                    capture_probability = greta::uniform(min = 0, max = 1, dim = 1))
  # # but need to expand dimensions of each
  
  # initialise empty parameters list
  parameters <- list()
  
  if (type == 'IPM') {
    
    params <- lapply(seq_len(replicates),
                     integrated_ipm(classes = classes,
                                    gp_tol = 1e-5))
    
    mu_initial <- lapply(seq_len(replicates),
                         function(x) greta::lognormal(meanlog = 0.0,
                                                      sdlog = 2.0,
                                                      dim = classes))
    
    parameters$survival <- lapply(seq_len(replicates), function(i) params[[i]]$survival) 
    parameters$survival_vec <- lapply(seq_len(replicates), function(i) params[[i]]$survival_vec) 
    parameters$fecundity <- lapply(seq_len(replicates), function(i) params[[i]]$fecundity) 
    parameters$density_parameter <- greta::uniform(min = params_list$density_lower,
                                                   max = params_list$density_upper, dim = 1)
    
    parameters$capture_probability <- lapply(seq_len(replicates),
                                             function(i) greta::uniform(min = params_list$capture_lower,
                                                                        max = params_list$capture_upper,
                                                                        dim = c(classes, 1)))
    structure <- 'IPM'
    
  }
  
  if (type == 'MPM') {
    
    # create transition matrix    
    params <- lapply(seq_len(replicates), function(x) get(structure)(classes = classes,
                                                                     params = params_list))
    
    mu_initial <- lapply(seq_len(replicates),
                         function(x) greta::lognormal(meanlog = 0.0,
                                                      sdlog = 4.0,
                                                      dim = classes))
    
    parameters$survival <- lapply(seq_len(replicates), function(i) params[[i]]$survival) 
    parameters$survival_vec <- lapply(seq_len(replicates), function(i) params[[i]]$survival_vec) 
    parameters$fecundity <- lapply(seq_len(replicates), function(i) params[[i]]$fecundity) 
    parameters$density_parameter <- lapply(seq_len(replicates),
                                           function(i) greta::uniform(min = params_list$density_lower,
                                                                      max = params_list$density_upper, dim = 1))
    parameters$capture_probability <- lapply(seq_len(replicates),
                                             function(i) greta::uniform(min = params_list$capture_lower,
                                                                        max = params_list$capture_upper,
                                                                        dim = c(classes, 1)))
    
  }
  
  integrated_process <- list(parameters = parameters,
                             mu_initial = mu_initial,
                             type = type,
                             structure = structure,
                             density_dependence = density_dependence,
                             classes = classes,
                             replicates = replicates,
                             replicate_id = replicate_id)
  
  integrated_process
  
}

#' @rdname integrated_process
#'
#' @export
#' 
#' @examples
#'
#' # check if an object is an integrated_process object
#'   
#' \dontrun{
#' is.integrated_process(model)
#' }

is.integrated_process <- function (model) {
  inherits(model, 'integrated_process')
}

#' @rdname integrated_process
#'
#' @export
#'
#' @examples
#' 
#' # Print information about an 'integrated_process' object
#'
#' \dontrun{
#' print(x)
#' }

print.integrated_process <- function (x, ...) {
  cat(paste0('This is an integrated_process object\n'))
}

#' @rdname integrated_process
#'
#' @export
#'
#' @examples
#' 
#' # Summarise an 'integrated_process' object
#'
#' \dontrun{
#' summary(x)
#' }

summary.integrated_process <- function (object, ...) {
  
  NULL
  
}


# internal function: create integrated_process object
as.integrated_process <- function (model) {
  as_class(model, name = 'integrated_process', type = 'list')
}

# internal function: create stage-structured matrix model with faster array setup
stage <- function (classes, params) {
  
  # survival prior
  if (is.matrix(params$surv_mat_param1)) {
    if (ncol(params$surv_mat_param1) != nrow(params$surv_mat_param1)) {
      stop('number of columns in matrix survival parameters must match number of rows',
           call. = FALSE)
    }
    if (ncol(params$surv_mat_param1) != classes) {
      stop('number of columns and rows in matrix survival parameters must match number of classes',
           call. = FALSE)
    }
    surv_mat_param1 <- params$surv_mat_param2
  } else {
    if (length(params$surv_mat_param1) != 1) {
      stop('matrix survival parameters must be a matrix or scalar value',
           call. = FALSE)
    }
    surv_mat_param1 <- matrix(params$surv_mat_param1, nrow = classes, ncol = classes)
  }
  if (is.matrix(params$surv_mat_param2)) {
    if (ncol(params$surv_mat_param2) != nrow(params$surv_mat_param2)) {
      stop('number of columns in matrix survival parameters must match number of rows',
           call. = FALSE)
    }
    if (ncol(params$surv_mat_param2) != classes) {
      stop('number of columns and rows in matrix survival parameters must match number of classes',
           call. = FALSE)
    }
    surv_mat_param2 <- params$surv_mat_param2
  } else {
    if (length(params$surv_mat_param1) != 1) {
      stop('matrix survival parameters must be a matrix or scalar value',
           call. = FALSE)
    }
    surv_mat_param2 <- matrix(1e6, nrow = classes, ncol = classes)
    diag(surv_mat_param2) <- params$surv_mat_param2
    for (j in seq_len(classes)) {
      if (j < classes) {
        surv_mat_param2[j + 1, j] <- params$surv_mat_param2
        if (j < (classes - 1)) {
          surv_mat_param2[j + 2, j] <- params$surv_mat_param2
        }
      } 
    }
  }
  
  # unpack survival parameters
  surv_param1 <- rep(1, classes)
  surv_param2 <- rep(1, classes)
  if (!is.null(params$surv_param1)) {
    if (length(params$surv_param1) == 1)
      params$surv_param1 <- rep(params$surv_param1, classes)
    surv_param1 <- params$surv_param1
  }
  if (!is.null(params$surv_param2)) {
    if (length(params$surv_param2) == 1)
      params$surv_param2 <- rep(params$surv_param2, classes)
    surv_param2 <- params$surv_param2
  }

  # standardise survival matrix
  survival_vec <- greta::beta(shape1 = surv_param1, shape2 = surv_param2)
  survival_tmp <- greta::beta(shape1 = surv_mat_param1, shape2 = surv_mat_param2)
  survival <- greta::sweep(survival_tmp, 2, colSums(survival_tmp), '/')
  
  # fecundity prior
  fec_min <- matrix(0.0, nrow = classes, ncol = classes)
  fec_min[1, classes] <- params$fec_lower
  fec_max <- matrix(0.0001, nrow = classes, ncol = classes)
  fec_max[1, classes] <- params$fec_upper
  fec_max <- greta::as_data(fec_max)
  fecundity <- fec_min + (fec_max - fec_min) * greta::uniform(min = 0, max = 1,
                                                              dim = c(classes, classes)) 
  
  # collate outputs
  params <- list(survival = survival,
                 fecundity = fecundity,
                 survival_vec = survival_vec)
  
  params
  
}

# internal function: create age-structured matrix model with faster array setup
age <- function (classes, params) {
  
  # survival prior
  surv_max <- matrix(0.0001, nrow = classes, ncol = classes)
  diag(surv_max) <- 1
  for (j in seq_len(classes)) {
    if (j < classes) {
      surv_max[j + 1, j] <- 1
    }
  }
  
  # standardise survival matrix
  surv_max <- greta::as_data(surv_max)
  survival_vec <- greta::uniform(min = 0, max = 1, dim = c(1, classes))
  survival <- surv_max * greta::uniform(min = 0, max = 1,
                                        dim = c(classes, classes))
  survival <- greta::sweep(survival, 2, colSums(survival), '/')
  survival <- greta::sweep(survival, 2, survival_vec, '*')
  
  # fecundity prior
  fec_min <- matrix(0.0, nrow = classes, ncol = classes)
  fec_min[1, classes] <- params$fec_lower
  fec_max <- matrix(0.0001, nrow = classes, ncol = classes)
  fec_max[1, classes] <- params$fec_upper
  fec_max <- greta::as_data(fec_max)
  fecundity <- fec_min + (fec_max - fec_min) * greta::uniform(min = 0, max = 1,
                                                              dim = c(classes, classes)) 
  
  # collate outputs
  params <- list(survival = survival,
                 fecundity = fecundity)
  
  params
  
}

# internal function: create unstructured matrix model with faster array setup
unstructured <- function (classes, params) {
  
  # survival prior
  surv_max <- matrix(0.0001, nrow = classes, ncol = classes)
  diag(surv_max) <- 1
  surv_max[lower.tri(surv_max)] <- 1
  
  
  # standardise survival matrix
  surv_max <- greta::as_data(surv_max)
  survival_vec <- greta::uniform(min = 0, max = 1, dim = c(1, classes))
  survival <- surv_max * greta::uniform(min = 0, max = 1,
                                        dim = c(classes, classes))
  survival <- greta::sweep(survival, 2, colSums(survival), '/')
  survival <- greta::sweep(survival, 2, survival_vec, '*')
  
  # fecundity prior
  fec_min <- matrix(0.0, nrow = classes, ncol = classes)
  fec_min[1, seq_len(classes)[-1]] <- rep(params$fec_lower, times = (classes - 1))
  fec_max <- matrix(0.0001, nrow = classes, ncol = classes)
  fec_max[1, seq_len(classes)[-1]] <- rep(params$fec_upper, times = (classes - 1))
  fec_max <- greta::as_data(fec_max)
  fecundity <- fec_min + (fec_max - fec_min) * greta::uniform(min = 0, max = 1,
                                                              dim = c(classes, classes)) 
  
  # collate outputs
  params <- list(survival = survival,
                 fecundity = fecundity)
  
  params
  
}

# internal function: create an IPM evaluated at `classes`
integrated_ipm <- function (classes, gp_tol) {
  
  # ipm transition kernel
  ipm_len <- greta::lognormal(mean = 0.0, sd = 1.0, dim = 1)
  ipm_sigma <- greta::lognormal(mean = 0.0, sd = 3.0, dim = 1)
  ipm_kern <- greta.gp::rbf(lengthscales = ipm_len, variance = ipm_sigma)
  ipm_mean <- greta::ilogit(greta.gp::gp(x = seq_len(classes),
                                         kernel = ipm_kern,
                                         tol = gp_tol))
  ipm_sigma_main <- greta::lognormal(mean = 0.0, sd = 3.0, dim = 1)
  
  # survival kernel
  surv_len <- greta::lognormal(mean = 0.0, sd = 1.0, dim = 1)
  surv_sigma <- greta::lognormal(mean = 0.0, sd = 3.0, dim = 1)
  surv_kern <- greta.gp::rbf(lengthscales = surv_len, variance = surv_sigma)
  ipm_surv <- greta::ilogit(greta.gp::gp(x = seq_len(classes),
                                         kernel = surv_kern,
                                         tol = gp_tol))
  
  # fecundity kernel
  fec_len <- greta::lognormal(mean = 0.0, sd = 1.0, dim = 1)
  fec_sigma <- greta::lognormal(mean = 0.0, sd = 3.0, dim = 1)
  fec_kern <- greta.gp::rbf(lengthscales = fec_len, variance = fec_sigma)
  ipm_fec <- exp(greta.gp::gp(x = seq_len(classes),
                              kernel = fec_kern,
                              tol = gp_tol))
  
  # convert flattened vector into matrix
  bins_mat <- matrix(rep(seq_len(classes), times = classes), ncol = classes)
  ipm_dist <- greta::sweep(greta::as_data(bins_mat), 2, ipm_mean, '-')
  ipm <- exp(-(ipm_dist * ipm_dist) / (2 * ipm_sigma_main))
  
  # standardise so that columns sum to one
  ipm <- greta::sweep(ipm, 2, greta::colSums(ipm), '/')
  
  params <- list(survival = ipm,
                 survival_vec = ipm_surv,
                 fecundity = ipm_fec)
  
  params
  
}
