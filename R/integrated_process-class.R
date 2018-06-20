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
#' 
#' @examples
#' 
#' library(integrated)
#' 
#' # prepare an example model
#' data <- define_integrated_process(type = 'MPM')

define_integrated_process <- function (type, structure, classes, density_dependence, replicates = 1) {
  
  if (!(type %in% c('IPM', 'MPM'))) {
    stop('type must be one of IPM or MPM')
  }
  
  if (type == 'IPM') {
    stop('IPM models not currently implemented')
  }
  
  if (type == 'MPM') {
    
    # set density dependence prior
    density_parameter <- greta::uniform(min = 0.8, max = 1.2, dim = 1)

    # create transition matrix    
    params <- get(structure)(classes = classes, replicates = replicates)
    parameters <- list(transitions = vector('list', length = replicates),
                       standard_deviations = params$standard_deviations)
    # mu_initial <- vector('list', length = replicates)
    mu_initial <- lapply(seq_len(replicates), function(x) greta::lognormal(meanlog = 0.0, sdlog = 2.0, dim = classes))
    parameters$transitions <- lapply(seq_len(replicates), function(x) matrix(params$transitions[, , x],
                                                                             nrow = classes, ncol = classes))
    # for (i in seq_len(replicates)) {

      # convert paramters to transition matrices
      # parameters$transitions[[i]] <- array(do.call('c', lapply(params$transitions, function(x) x[i])),
      #                                      dim = c(classes, classes))

      # initial population abundances
      # mu_initial[[i]] <- greta::lognormal(meanlog = 0.0, sdlog = 2.0, dim = classes)
      
    # }
    
    
  }
  
  integrated_process <- list(parameters = parameters,
                             mu_initial = mu_initial,
                             type = type,
                             structure = structure,
                             classes = classes,
                             density_dependence = density_dependence,
                             density_parameter = density_parameter,
                             replicates = replicates)
  
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

# internal function: create stage-structured matrix model
stage <- function(classes, replicates) {
  
  demo_sd <- greta::lognormal(mean = rep(0.0, 3), sd = c(3.0, 3.0, 0.0001), dim = 3)
  
  array_mean <- array(data = -10, dim = c(classes, classes, replicates))
  array_mean_fec <- array(data = 0, dim = c(classes, classes, replicates))
  for (i in seq_len(classes)) {
    array_mean[i, i, ] <- rep(0, replicates)
    if (i < classes) {
      array_mean[(i + 1), i, ] <- rep(0, replicates)
      if (i < (classes - 1)) {
        array_mean[(i + 2), i, ] <- rep(0, replicates)
      }
    }
  }
  
  array_sd <- array(rep(demo_sd[1], (classes * classes * replicates)),
                    dim = dim(array_mean))
  array_sd_fec <- array(rep(c(rep(demo_sd[3], (classes - 1)), demo_sd[2],
                              rep(demo_sd[3], ((classes - 1) * classes))), replicates),
                        dim = dim(array_mean))
  
  # define parameters  
  surv_params <- greta::ilogit(greta::normal(mean = array_mean,
                                             sd = array_sd,
                                             dim = dim(array_mean)))
  fec_params <- greta::lognormal(meanlog = array_mean_fec,
                                 sdlog = array_sd_fec,
                                 dim = dim(array_mean_fec))
  
  # return outputs
  transitions <- surv_params + fec_params
  
  params <- list(transitions = transitions,
                 standard_deviations = demo_sd)
  
  params
  
}

# internal function: create age-structured matrix model
age <- function(classes, replicates) {
  
  # hyperpriors for sds
  demo_sd <- greta::lognormal(mean = 0.0, sd = 3.0, dim = (classes + 1))
  
  # fecundity and survival priors
  params <- vector("list", length = (classes * classes))
  for (i in 1:(classes - 1)) {
    params[[i]] <- greta::uniform(min = 0.0, max = 0.0001, dim = replicates)
  }
  params[[classes]] <- greta::lognormal(mean = 0.0, sd = demo_sd[1], dim = replicates)
  for (i in seq_len(classes)[-1]) {
    for (j in seq_len(classes)) {
      if (j == (i - 1)) {
        params[[((i - 1) * classes) + j]] <- greta::ilogit(greta::normal(mean = 0.0, sd = demo_sd[i],
                                                                         dim = replicates))
      } else {
        if ((i == classes) & (j == classes)) {
          params[[((i - 1) * classes) + j]] <- greta::ilogit(greta::normal(mean = 0.0, sd = demo_sd[(classes + 1)],
                                                                           dim = replicates))
        } else {
          params[[((i - 1) * classes) + j]] <- greta::uniform(min = 0.0, max = 0.0001, dim = replicates)
        }
      }
    }
  }
  
  # collate outputs
  out <- list(transitions = params, standard_deviations = demo_sd)
  
  # return outputs  
  out
  
}

# internal function: create unstructured matrix model
unstructured <- function(classes, replicates) {
  
  # hyperpriors for sds
  demo_sd <- greta::lognormal(mean = 0.0, sd = 3.0, dim = (classes + 1))
  
  # fecundity and survival priors
  params <- vector("list", length = (classes * classes))
  params[[1]] <- greta::ilogit(greta::normal(mean = 0.0, sd = demo_sd[1], dim = replicates))
  for (i in 2:nstage) {
    params[[i]] <- greta::lognormal(mean = 0.0, sd = demo_sd[2], dim = replicates)
  }
  for (i in (classes + 1):length(params)) {
    params[[i]] <- greta::ilogit(greta::normal(mean = 0.0,
                                               sd = demo_sd[(floor((i - 1) / classes) + 2)],
                                               dim = replicates))
  }
  
  # collate outputs
  out <- list(transitions = params,
              standard_deviations = demo_sd)
  
  # return outputs
  out
  
}
