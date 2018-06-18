#' An integrated data object
#'
#' @description A \code{integrated_model} call adds a component data object to process model
#'   created with \link[integrated]{define_integrated_process}
#' 
#' @rdname integrated_data
#' 
#' @param data something
#' @param process_model
#' @param process_link
#' @param observation_model something else
#' @param ... additional arguments
#'
#' @return An object of class \code{integrated_data}, which contains information on the data module and
#'   can be passed to a \link[integrated]{integrated_model}
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
#' data <- define_integrated_data()
#'                         
#' \dontrun{                 
#' # summarise data module
#' model
#' summary(model)
#' plot(model)
#' }



## COULD BE "BUILD/DEFINE integrated_data" if we can return greta_arrays
add_integrated_data <- function (data,
                                 process_model,
                                 process_link,
                                 observation_model = 'naive') {
  
  if (!(process_link %in% c())) {
    stop('process_link must be a known process module')
  } 
  
  if (process_link == 'growth') {
    add_growth_module(data = data,
                      process_model = process_model,
                      observation_model = observation_model)
  } 
  
  if (process_link == 'abundance') {
    add_abundance_module(data = data,
                         process_model = process_model,
                         observation_model = observation_model)
  } 
  
  if (process_link == 'mark_recapture') {
    add_mark_recapture_module(data = data,
                              process_model = process_model,
                              observation_model = observation_model)
  } 
  
  if (process_link == 'size_abundance') {
    add_size_abundance_module(data = data,
                              process_model = process_model,
                              observation_model = observation_model)
  }  
  
  if (process_link == 'biomass') {
    add_biomass_module(data = data,
                       process_model = process_model,
                       observation_model = observation_model)
  }  
  
  if (process_link == 'community') {
    add_community_module(data = data,
                         process_model = process_model,
                         observation_model = observation_model)
  }  
  
  data_module <- list(data = data,
                      process_model = process_model,
                      process_link = process_link,
                      observation_model = observation_model)
  
  as.integrated_data(data_module)
  
}

#' @rdname integrated_data
#'
#' @export
#' 
#' @examples
#'
#' # check if an object is an integrated_data object
#'   
#' \dontrun{
#' is.integrated_data(model)
#' }

is.integrated_data <- function (model) {
  inherits(model, 'integrated_data')
}

#' @rdname integrated_data
#'
#' @export
#'
#' @examples
#' 
#' # Print information about an 'integrated_data' object
#'
#' \dontrun{
#' print(x)
#' }

print.integrated_data <- function (x, ...) {
  cat(paste0('This is an integrated_data object\n'))
}

#' @rdname integrated_data
#'
#' @export
#'
#' @examples
#' 
#' # Plot an 'integrated_data' object
#'
#' \dontrun{
#' plot(x)
#' }

plot.integrated_data <- function (x, ...) {

  plot(x$greta_model, ...)

}

#' @rdname integrated_data
#'
#' @export
#'
#' @examples
#' 
#' # Summarise an 'integrated_data' object
#'
#' \dontrun{
#' summary(x)
#' }

summary.integrated_data <- function (object, ...) {
  
  NULL
  
}


# internal function: create integrated_data object
as.integrated_data <- function (model) {
  as_class(model, name = 'integrated_data', type = 'list')
}

# internal function: build growth data module
build_growth_module <- function (data, process_model, observation_model) {
  
  for (i in seq_len(process_model$replicates)) {
    for (j in seq_len(ncol(data))) {
      data_tmp <- matrix(data[, k], ncol = process_model$classes)
      distribution(data_tmp) <- multinomial(size = sum(data_tmp),
                                            prob = process_model$parameters$transitions[[i]][, j],
                                            dim = 1)
    }
  }
  
}

# internal function: build abundance data module
build_abundance_module <- function (data, process_model, observation_model) {
  
  # create output lists
  mu_iterated <- vector("list", length = process_model$replicates)
  
  for (i in seq_len(process_model$replicates)) {
    mu_iterated[[i]] <- iterate_state(t(process_model$transitions[[i]]),
                                      process_model$mu_initial[[i]],
                                      process_model$density_parameter,
                                      seq_len(ncol(data[[i]])),
                                      dens_form = process_model$density_dependence)
  }
  
  mu_flattened <- do.call("c", mu_iterated)
  distribution(data) <- poisson(mu_flattened)
  
}

# internal function: build mark-recapture data module
build_mark_recapture_module <- function (data, process_model, observation_model) {
  
  NULL
  
}

# internal function: build size-abundance data module
build_size_abundance_module <- function (data, process_model, observation_model) {
  
  data_module <- NULL
  
  data_module
  
}

# internal function: build biomass data module
build_biomass_module <- function (data, process_model, observation_model) {
  
  data_module <- NULL
  
  data_module
  
}

# internal function: build community data module
build_community_module <- function (data, process_model, observation_model) {
  
  data_module <- NULL
  
  data_module
  
}
