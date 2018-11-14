#' @name integrated_process
#' @title create integrated process objects
#'
#' @description An \code{integrated_process} object contains the underlying
#'   process model for an integrated population analysis
#' 
#' @param classes something
#' @param density something
#' @param replicates number of distinct populations
#' @param params named list of parameters (see details for information on setting prior distributions)
#' @param ... additional arguments
#' @param x an \code{integrated_process} object
#' @param object an \code{integrated_process} object
#'
#' @details something. Prior distributions can be specified as single-dimensional
#'   greta distribution, e.g., \code{normal(0, 1)}. Link functions and transformations
#'   can be specified directly in-line, e.g., \code{ilogit(normal(0, 1))} specifies
#'   normal priors with a mean of zero and a standard deviation of one, transformed
#'   with an inverse-logit link.
#'
#' @return An object of class \code{integrated_process}, which can be used to create
#'    \link[integrated]{integrated_data} and \link[integrated]{integrated_model} objects
#' 
#' @import greta
#' 
#' @examples
#' \dontrun {
#' 
#' library(integrated)
#' 
#' # a really basic age-structured model with five age classes
#' process <- leslie(5, density = "none")
#' 
#' # setting custom priors
#' process <- leslie(5, density = "none",
#'                   params = list(survival = iprobit(normal(0, 1)),
#'                                 fecundity = exp(normal(0, 1))))
#' }


#' @export
#' @rdname integrated_process
#' 
leslie <- function(classes, density = "none", params = list()) {
  
  # set type
  type <- "leslie"
  
  # initialise model parameters
  param_list <- list(density = uniform(0, 1),
                     reproductive = max(classes),
                     survival = beta(1, 1),
                     fecundity = normal(0, 10, truncation = c(0, Inf)),
                     initials = normal(0, 10, truncation = c(0, Inf)),
                     random = normal(0, 10, truncation = c(0, Inf)))
  param_list[names(params)] <- params

  # do the parameters have reasonable bounds?
  survival_bounds <- extract_bounds(param_list$survival)
  fecundity_bounds <- extract_bounds(param_list$survival)
  
  # warn if not
  if (survival_bounds[1] < 0 | survival_bounds[2] > 1)
    warning("the prior for survival has bounds outside of [0, 1]; is this reasonable?", call. = FALSE)
  if (fecundity_bounds[1] < 0)
    warning("the prior for fecundity has a lower bound less than 0; is this reasonable?", call. = FALSE)
  
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
                     reproductive = max(classes),
                     survival = beta(1, 1),
                     growth = beta(1, 1),
                     fecundity = normal(0, 10, truncation = c(0, Inf)),
                     initials = normal(0, 10, truncation = c(0, Inf)),
                     random = normal(0, 10, truncation = c(0, Inf)))
  param_list[names(params)] <- params
  
  # do the parameters have reasonable bounds?
  survival_bounds <- extract_bounds(param_list$survival)
  growth_bounds <- extract_bounds(param_list$growth)
  fecundity_bounds <- extract_bounds(param_list$survival)
  
  # warn if not
  if (survival_bounds[1] < 0 | survival_bounds[2] > 1)
    warning("the prior for survival has bounds outside of [0, 1]; is this reasonable?", call. = FALSE)
  if (growth_bounds[1] < 0 | growth_bounds[2] > 1)
    warning("the prior for growth has bounds outside of [0, 1]; is this reasonable?", call. = FALSE)
  if (fecundity_bounds[1] < 0)
    warning("the prior for fecundity has a lower bound less than 0; is this reasonable?", call. = FALSE)
  
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
                     reproductive = max(classes),
                     survival = beta(1, 1),
                     growth = beta(1, 1),
                     fecundity = normal(0, 10, truncation = c(0, Inf)),
                     initials = normal(0, 10, truncation = c(0, Inf)),
                     random = normal(0, 10, truncation = c(0, Inf)))
  param_list[names(params)] <- params
  
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
is.integrated_process <- function(x) {
  inherits(model, "integrated_process")
}

#' @export
#' @rdname integrated_process
#' 
print.integrated_process <- function(x, ...) {
  cat(paste0("This is an integrated_process object\n"))
}

#' @export
#' @rdname integrated_process
#' 
summary.integrated_process <- function(object, ...) {
  
  NULL
  
}

# internal function: create integrated_process object
as.integrated_process <- function (model) {
  as_class(model, name = "integrated_process", type = "list")
}
