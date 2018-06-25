#' An integrated process object
#'
#' @description A \code{integrated_process} object contains a component data object for integrated data analysis
#' 
#' @rdname integrated_process
#' 
#' @param type something
#' @param structure something
#' @param classes something
#' @param density_dependence something
#' @param replicates number of distinct populations
#' @param ... additional arguments
#'
#' @return An object of class \code{integrated_process}, which can be passed to
#'    \link[integrated]{define_integrated_data} and \link[integrated]{build_integrated_model}
#' 
#' @export
#' 
#' @import greta
#' @import greta.gp
#' 
#' @examples
#' 
#' library(integrated)
#' 
#' # prepare an example model
#' data <- define_integrated_process(type = 'MPM')

define_integrated_process <- function (type, classes,
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
                      density_lower = 0,
                      density_upper = 5)
  params_list[names(params)] <- params
  
  if (type == 'IPM') {
    
    parameters <- list()
    parameters$transitions <- build_integrated_ipm(classes = classes,
                                                   replicates = replicates,
                                                   gp_tol = 1e-5)
    
    mu_initial <- lapply(seq_len(replicates),
                         function(x) greta::lognormal(meanlog = 0.0,
                                                      sdlog = 2.0,
                                                      dim = classes))
    
    structure <- 'IPM'

  }
  
  if (type == 'MPM') {
    
    # create transition matrix    
    params <- get(structure)(classes = classes,
                             replicates = replicates,
                             params = params_list)
    parameters <- list(transitions = vector('list', length = replicates),
                       standard_deviations = params$standard_deviations)
    mu_initial <- lapply(seq_len(replicates),
                         function(x) greta::lognormal(meanlog = 0.0,
                                                      sdlog = 4.0,
                                                      dim = classes))
    parameters$survival <- lapply(seq_len(replicates),
                                  function(i) array(params$survival[, , i],
                                                    dim = c(classes, classes)))
    parameters$fecundity <- lapply(seq_len(replicates),
                                  function(i) array(params$fecundity[, , i],
                                                    dim = c(classes, classes)))
    parameters$density_parameter <- greta::uniform(min = params_list$density_lower,
                                                   max = params_list$density_upper, dim = 1)
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
stage <- function (classes, replicates, params) {
  
  # survival prior
  surv_max <- array(0.0001, dim = c(classes, classes, replicates))
  for (i in seq_len(replicates)) {
    
    diag(surv_max[, , i]) <- 1
    
    for (j in seq_len(classes)) {
      if (j < classes) {
        surv_max[j + 1, j, i] <- 1
        if (j < (classes - 1)) {
          surv_max[j + 2, j, i] <- 1
        }
      } 
    }
    
  } 
  surv_max <- greta::as_data(surv_max)
  survival <- surv_max * greta::uniform(min = 0, max = 1,
                                        dim = dim(surv_max))
  
  # fecundity prior
  fec_min <- array(0.0, dim = c(classes, classes, replicates))
  fec_min[1, classes, ] <- rep(params$fec_lower, times = replicates)
  fec_max <- array(0.0001, dim = c(classes, classes, replicates))
  fec_max[1, classes, ] <- rep(params$fec_upper, times = replicates)
  fec_max <- greta::as_data(fec_max)
  fecundity <- fec_min + (fec_max - fec_min) * greta::uniform(min = 0, max = 1,
                                                              dim = dim(fec_max)) 
  
  # collate outputs
  params <- list(survival = survival,
                 fecundity = fecundity)
  
  params
  
}

# internal function: create age-structured matrix model with faster array setup
age <- function (classes, replicates, params) {
  
  # survival prior
  surv_max <- array(0.0001, dim = c(classes, classes, replicates))
  for (i in seq_len(replicates)) {
    
    diag(surv_max[, , i]) <- 1
    
    for (j in seq_len(classes)) {
      if (j < classes) {
        surv_max[j + 1, j, i] <- 1
      }
    }

  } 
  surv_max <- greta::as_data(surv_max)
  survival <- surv_max * greta::uniform(min = 0, max = 1,
                                         dim = dim(surv_max))
  
  # fecundity prior
  fec_min <- array(0.0, dim = c(classes, classes, replicates))
  fec_min[1, classes, ] <- rep(params$fec_lower, times = replicates)
  fec_max <- array(0.0001, dim = c(classes, classes, replicates))
  fec_max[1, classes, ] <- rep(params$fec_upper, times = replicates)
  fec_max <- greta::as_data(fec_max)
  fecundity <- fec_min + (fec_max - fec_min) * greta::uniform(min = 0, max = 1,
                                                              dim = dim(fec_max)) 
  
  # collate outputs
  params <- list(survival = survival,
                 fecundity = fecundity)
  
  params
  
}

# internal function: create unstructured matrix model with faster array setup
unstructured <- function (classes, replicates, params) {
  
  # survival prior
  surv_max <- array(0.0001, dim = c(classes, classes, replicates))
  for (i in seq_len(replicates)) {
    diag(surv_max[, , i]) <- 1
    surv_max[, , i][lower.tri(surv_max[, , i])] <- 1
    
  } 
  surv_max <- greta::as_data(surv_max)
  survival <- surv_max * greta::uniform(min = 0, max = 1,
                                        dim = dim(surv_max))
  
  # fecundity prior
  fec_min <- array(0.0, dim = c(classes, classes, replicates))
  fec_min[1, seq_len(classes)[-1], ] <- rep(params$fec_lower, times = replicates)
  fec_max <- array(0.0001, dim = c(classes, classes, replicates))
  fec_max[1, seq_len(classes)[-1], ] <- rep(params$fec_upper, times = (replicates * (classes - 1)))
  fec_max <- greta::as_data(fec_max)
  fecundity <- fec_min + (fec_max - fec_min) * greta::uniform(min = 0, max = 1,
                                                              dim = dim(fec_max)) 
  
  # collate outputs
  params <- list(survival = survival,
                 fecundity = fecundity)
  
  params
  
}

# internal function: create an IPM evaluated at `classes`
build_integrated_ipm <- function (classes, replicates, gp_tol) {
  
  params <- vector('list', length = replicates)
  
  for (i in seq_len(replicates)) {
    
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
    
    # standardise to incorporate survival per column
    ipm <- greta::sweep(ipm, 2, ipm_surv, '*')
    
    # add fecundity elements to kernel
    ipm[1, ] <- ipm[1, ] + t(ipm_fec)
    
    params[[i]] <- ipm
    
  }
  
  params
   
}
