#' An integrated model object
#'
#' @description A \code{integrated_model} object contains a compiled greta model for integrated data analysis
#' 
#' @rdname integrated_model
#' 
#' @param ... \code{integrated_data} objects created with \link[integrated]{define_integrated_data}
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

build_integrated_model <- function (integrated_process, ...) {
  
  data_modules <- list(...)
  
  for (i in seq_along(data_modules)) {
    
    data_tmp <- data_modules[[i]]
    
    if (!(data_tmp$process_link %in% c('abundance', 'growth'))) {
      stop('only abundance and growth modules are currently implemented')
    }
    
    if (data_tmp$process_link == 'abundance') {
      data <- do.call('c', data_tmp$data)
      greta::distribution(data) <- greta::poisson(data_tmp$data_module)
    }
    
    if (data_tmp$process_link == 'growth') {
      
      for (i in seq_along(data_tmp$data_module)) {
        for (j in seq_along(data_tmp$data_module[[i]])) {
          greta::distribution(data_tmp$data_module[[i]][[j]]) <-
            greta::multinomial(size = sum(data_tmp$data_module[[i]][[j]]),
                               prob = t(integrated_process$parameters$transitions[[i]][j, ]),
                               dim = 1)
        }
      }
      
    }
    
  } 
  
  do.call('c', integrated_process$parameters$transitions)
  
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
