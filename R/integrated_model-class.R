#' An integrated model object
#'
#' @description A \code{integrated_model} object contains a compiled greta model for integrated data analysis
#' 
#' @rdname integrated_model
#' 
#' @param ... \code{integrated_data} objects created with \link[integrated]{build_integrated_data}
#'
#' @return An object of class \code{integrated_model}, which can be passed to
#'    \link[greta]{mcmc} and has associated `print`, `plot`, and `summary` methods
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
#' model <- build_integrated_model()
#'                         
#' \dontrun{                 
#' # summarise fitted model
#' model
#' summary(model)
#' plot(model)
#' }

build_integrated_model <- function (...) {
  
  component_data <- list(...)
  
  greta_model <- NULL
  
  out <- list(greta_model = greta_model,
              component_data = component_data)
  
  as.integrated_model(out)

}

#' @rdname integrated_model
#'
#' @export
#' 
#' @examples
#'
#' # check if an object is an integrated_model object
#'   
#' \dontrun{
#' is.integrated_model(model)
#' }

is.integrated_model <- function (model) {
  inherits(model, 'integrated_model')
}

#' @rdname integrated_model
#'
#' @export
#'
#' @examples
#' 
#' # Print information about an 'integrated_model' object
#'
#' \dontrun{
#' print(x)
#' }

print.integrated_model <- function (x, ...) {
  cat(paste0('This is an integrated_model object\n'))
}

#' @rdname integrated_model
#'
#' @export
#'
#' @examples
#' 
#' # Plot an 'integrated_model' object
#'
#' \dontrun{
#' plot(x)
#' }

plot.integrated_model <- function (x, ...) {

  plot(x$greta_model, ...)

}

#' @rdname integrated_model
#'
#' @export
#'
#' @examples
#' 
#' # Summarise an 'integrated_model' object
#'
#' \dontrun{
#' summary(x)
#' }

summary.integrated_model <- function (object, ...) {
  
  NULL
  
}


# internal function: create integrated_model object
as.integrated_model <- function (model) {
  as_class(model, name = 'integrated_model', type = 'list')
}
