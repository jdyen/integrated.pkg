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
    mu_initial <- lapply(seq_len(replicates), function(x) greta::lognormal(meanlog = 0.0, sdlog = 2.0, dim = classes))
    parameters$transitions <- lapply(seq_len(replicates), function(x) array(params$transitions[, , x],
                                                                             dim = c(classes, classes)))

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
  
  demo_sd <- greta::lognormal(mean = rep(0.0, 0.0, -50), sd = c(3.0, 3.0, 0.000001), dim = 3)
  
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
  
  demo_sd <- greta::lognormal(mean = rep(0.0, 0.0, -50), sd = c(3.0, 3.0, 0.000001), dim = 3)
  
  array_mean <- array(data = -10, dim = c(classes, classes, replicates))
  array_mean_fec <- array(data = 0, dim = c(classes, classes, replicates))
  array_mean[classes, classes, ] <- rep(0, replicates)
  for (i in seq_len(classes - 1)) {
    array_mean[(i + 1), i, ] <- rep(0, replicates)
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

# internal function: create unstructured matrix model
unstructured <- function(classes, replicates) {
  
  demo_sd <- greta::lognormal(mean = rep(0.0, 0.0, -50), sd = c(3.0, 3.0, 0.000001), dim = 3)
  
  array_mean <- array(data = 0, dim = c(classes, classes, replicates))
  array_mean_fec <- array(data = 0, dim = c(classes, classes, replicates))
  array_sd <- array(rep(demo_sd[1], (classes * classes * replicates)),
                    dim = dim(array_mean))
  array_sd_fec <- array(rep(c(rep(demo_sd[2], classes),
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
