#' An integrated model object
#'
#' @description An \code{integrated_model} object returns a \link[greta]{greta_array} that can be
#'  passed to \link[greta]{model}
#' 
#' @rdname integrated_model
#' 
#' @param ... \code{integrated_data} objects created with \link[integrated]{integrated_data}
#'
#' @return An object of class \code{greta_array}, which can be passed to
#'    \link[greta]{model} and has associated `print`, `plot`, and `summary` methods
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
#' model <- integrated_model()
#'                         
#' \dontrun{                 
#' # summarise fitted model
#' model
#' summary(model)
#' plot(model)
#' }

integrated_model <- function (integrated_process, ...) {
  
  data_modules <- list(...)
  
  # initialise mu values
  mu_param <- NULL
  
  for (i in seq_along(data_modules)) {
    
    data_tmp <- data_modules[[i]]
    
    if (!(data_tmp$process_link %in% c('stage_abundance', 'individual_growth', 'stage_recapture'))) {
      stop('only abundance, growth, and mark_recapture modules are currently implemented',
           call. = FALSE)
    }
    
    if (data_tmp$process_link == 'stage_abundance') {
      
      # it's easy if data and process contain the same number of stages
      if (all(integrated_process$classes == sapply(data_tmp$data, nrow))) {
        
        # flatten the data set into a vector
        data <- do.call('c', data_tmp$data)
        
        # poisson likelihood for observed abundances      
        greta::distribution(data) <- greta::poisson(data_tmp$data_module)
        
        # store mu values to include in params vector
        if (is.null(mu_param)) {
          mu_param <- data_tmp$data_module
        } else {
          mu_param <- c(mu_param, data_tmp$data_module)
        }
        
      } else {
        
        # we need to be more careful because data are binned more
        #  coarsely than the integrated process
        index <- NULL
        for (i in seq_along(data_tmp$data)) {
          
          # we want to keep the existing greta_arrays for each element
          #   where the process and daata have the same number of stages
          index_max <- ifelse(i == 1, 0, max(index))
          if (nrow(data_tmp$data[[i]]) == integrated_process$classes) {
            index <- c(index, seq_len(length(data_tmp$data[[i]])))
          } else {
            index <- c(index,
                       floor(seq(index_max + 1,
                                 index_max + nrow(data_tmp$data[[i]]) * ncol(data_tmp$data[[i]]),
                                 length = integrated_process$classes * ncol(data_tmp$data[[i]]))))
          }
          
        }
        
        # flatten the data into a vector
        data <- do.call('c', data_tmp$data)
        
        # collapse the process abundances into fewer classes
        collapsed_data_module <- tapply(data_tmp$data_module, index, 'sum')
        
        # poisson likelihood for observed abundances      
        greta::distribution(data) <- greta::poisson(collapsed_data_module)
        
        # store mu values to include in params vector
        if (is.null(mu_param)) {
          mu_param <- collapsed_data_module
        } else {
          mu_param <- c(mu_param, collapsed_data_module)
        }
        
      }
      
    }
    
    if (data_tmp$process_link == 'individual_growth') {
      
      # if there is more than one observed data set
      for (i in seq_along(data_tmp$data_module)) {
        
        # if there are multiple data elements and only one process matrix
        if (integrated_process$replicates == 1) {
          
          greta::distribution(data_tmp$data_module[[i]]) <-
            greta::multinomial(size = rowSums(data_tmp$data_module[[i]]),
                               prob = t(integrated_process$parameters$survival[[1]]))
          
        } else {
          
          # otherwise there must be one process matrix for each data element
          greta::distribution(data_tmp$data_module[[i]]) <-
            greta::multinomial(size = rowSums(data_tmp$data_module[[i]]),
                               prob = t(integrated_process$parameters$survival[[integrated_process$replicate_id[i]]]))
          
        }
        
      }
      
    }
    
    if (data_tmp$process_link == 'stage_recapture') {
      
      for (i in seq_along(data_tmp$data_module$count)) {
        
        # use separate process models if they exist
        if (integrated_process$replicates > 1) {
          
          greta::distribution(data_tmp$data_module$count[[i]]) <-
            greta::multinomial(size = rowSums(data_tmp$data_module$count[[i]]),
                               prob = t(integrated_process$parameters$survival[[integrated_process$replicate_id[i]]]))
          greta::distribution(data_tmp$data_module$count2[[i]]) <-
            greta::binomial(size = data_tmp$data_module$total[[i]],
                            prob = integrated_process$parameters$survival_vec[[integrated_process$replicate_id[i]]])

        } else {  
          
          greta::distribution(data_tmp$data_module$count[[i]]) <-
            greta::multinomial(size = rowSums(data_tmp$data_module$count[[i]]),
                               prob = t(integrated_process$parameters$survival[[1]]))
          greta::distribution(data_tmp$data_module$count2[[i]]) <-
            greta::binomial(size = data_tmp$data_module$total[[i]],
                            prob = integrated_process$parameters$survival_vec[[1]])

        }
        
      }
      
    }
    
  } 
  
  c(do.call('c', integrated_process$parameters$survival),
    do.call('c', integrated_process$parameters$fecundity),
    do.call('c', integrated_process$parameters$survival_vec),
    do.call('c', integrated_process$parameters$capture_probability),
    do.call('c', integrated_process$mu_initial),
    do.call('c', integrated_process$parameters$density_parameter),
    mu_param)
  
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
