#' An integrated data object
#'
#' @description A \code{integrated_model} call adds a component data object to process model
#'   created with \link[integrated]{define_integrated_process}
#' 
#' @rdname integrated_data
#' 
#' @param data something
#' @param process_model something
#' @param process_link something
#' @param observation_model something else
#' @param ... additional arguments
#'
#' @return An object of class \code{integrated_data}, which contains information on the data module and
#'   can be passed to a \link[integrated]{integrated_model}
#' 
#' @export
#' 
#' @import greta
#' @import gretaDynamics
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
define_integrated_data <- function (data,
                                 process_model,
                                 process_link,
                                 observation_model = 'naive') {
  
  if (!(process_link %in% c('growth', 'abundance',
                            'mark_recapture', 'size_abundance',
                            'biomass', 'community'))) {
    stop('process_link must be a known process module')
  } 
  
  if (process_link == 'growth') {
    data_module <- define_growth_module(data = data,
                      process_model = process_model,
                      observation_model = observation_model)
  } 
  
  if (process_link == 'abundance') {
    data_module <- define_abundance_module(data = data,
                         process_model = process_model,
                         observation_model = observation_model)
  } 
  
  if (process_link == 'mark_recapture') {
    data_module <- define_mark_recapture_module(data = data,
                              process_model = process_model,
                              observation_model = observation_model)
  } 
  
  if (process_link == 'size_abundance') {
    data_module <- define_size_abundance_module(data = data,
                              process_model = process_model,
                              observation_model = observation_model)
  }  
  
  if (process_link == 'biomass') {
    data_module <- define_biomass_module(data = data,
                       process_model = process_model,
                       observation_model = observation_model)
  }  
  
  if (process_link == 'community') {
    data_module <- define_community_module(data = data,
                         process_model = process_model,
                         observation_model = observation_model)
  }  
  
  data_module <- list(data_module = data_module,
                      data = data,
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
define_growth_module <- function (data, process_model, observation_model) {
  
  size_data <- vector('list', length = replicates)
  for (i in seq_len(process_model$replicates)) {
    size_data_tmp <- vector('list', length = ncol(data))
    for (j in seq_len(ncol(data))) {
      size_data_tmp[[j]] <- matrix(data[, j], ncol = process_model$classes)
    }
    size_data[[i]] <- do.call('c', size_data_tmp)
  }
  
  size_data

}

# internal function: build abundance data module
define_abundance_module <- function (data, process_model, observation_model) {
  
  # create output lists
  mu_iterated <- vector("list", length = process_model$replicates)
  
  for (i in seq_len(process_model$replicates)) {
    mu_iterated[[i]] <- iterate_state(t(process_model$parameters$transitions[[i]]),
                                      process_model$mu_initial[[i]],
                                      process_model$density_parameter,
                                      seq_len(ncol(data[[i]])),
                                      dens_form = process_model$density_dependence)
  }
  
  mu_flattened <- do.call('c', mu_iterated)
  
  mu_flattened
  
}

# internal function: build mark-recapture data module
define_mark_recapture_module <- function (data, process_model, observation_model) {
  
  NULL
  
}

# internal function: build size-abundance data module
define_size_abundance_module <- function (data, process_model, observation_model) {
  
  data_module <- NULL
  
  data_module
  
}

# internal function: build biomass data module
define_biomass_module <- function (data, process_model, observation_model) {
  
  data_module <- NULL
  
  data_module
  
}

# internal function: build community data module
define_community_module <- function (data, process_model, observation_model) {
  
  data_module <- NULL
  
  data_module
  
}
