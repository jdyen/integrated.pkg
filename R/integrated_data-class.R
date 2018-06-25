#' An integrated data object
#'
#' @description A \code{integrated_model} call adds a component data object to process model
#'   created with \link[integrated]{define_integrated_process}
#' 
#' @rdname integrated_data
#' 
#' @param data something
#' @param integrated_process something
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
#' @import greta.dynamics
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

define_integrated_data <- function (data,
                                    integrated_process, 
                                    process_link,
                                    observation_model = 'naive') {
  
  if (!(process_link %in% c('growth', 'abundance',
                            'mark_recapture', 'size_abundance',
                            'biomass', 'community'))) {
    stop('process_link must be a known process module',
         call. = FALSE)
  }  
  
  if (process_link == 'growth') {
    
    # convert matrix or data.frame data to a list
    if (is.matrix(data) | is.data.frame(data)) {
      
      # check if data are formatted correctly
      if (ncol(data) != nrow(data)) {
        data <- make_growth_data_matrix(data = data,
                                        classes = integrated_process$classes)
      }
      
      data <- list(data)
      
    }
    
    # data should be a list or matrix
    if (!is.list(data)) {
      stop('growth data must be a matrix, data.frame, or list of matrices or data.frames',
           call. = FALSE)
    }
    
    # error if the number of data classes doesn't match the process
    if (!all(integrated_process$classes == sapply(data, nrow))) {
      stop('number of classes in growth data must match number of classes in integrated_process')
    }
    
    if (integrated_process$replicates > 1) {      
      if (length(data) > 1) {
        if (length(data) != length(integrated_process$replicate_id)) {
          stop(paste0('growth data have ', length(data), ' elements but ',
                      ' integrated process has ', integrated_process$replicates,
                      ' replicates'),
               call. = FALSE)
        }
      } else {
        stop(paste0('one growth data set should be supplied for each of the ',
                    integrated_process$replicates, ' process replicates'),
             call. = FALSE)
      }
    }

    # create data module from growth data matrix
    data_module <- define_growth_module(data = data,
                                        integrated_process = integrated_process, 
                                        observation_model = observation_model)
  }   
  
  if (process_link == 'abundance') {
    
    # check data format
    if (is.matrix(data) | is.data.frame(data)) {
      data <- list(data)
    }
    
    # error if not a list or matrix
    if (!is.list(data)) {
      stop('abundance data must be a list with classes in rows and samples in columns',
           call. = FALSE) 
    }
    
    # if there is more than one replicate
    if (integrated_process$replicates > 1) {
      # check that there is one data element for each replicate
      if (length(data) != length(integrated_process$replicate_id)) {
        stop(paste0('abundance data have ', length(data), ' elements but ',
                    ' integrated_process contains ', integrated_process$replicates,
                    ' replicates'),
             call. = FALSE)
      }
    }
    
    # this won't work if there are more classes in data than in the process model
    if (max(sapply(data, nrow)) > integrated_process$classes) {
      stop(paste0('there are up to ', max(sapply(data, nrow)), ' classes in the data set ',
                  'but only ', integrated_process$classes, ' classes in integrated_process'),
           call. = FALSE)
    }
    
    # create data module from list data
    data_module <- define_abundance_module(data = data,
                                           integrated_process = integrated_process, 
                                           observation_model = observation_model)
  }  
  
  if (process_link == 'mark_recapture') {
    data_module <- define_mark_recapture_module(data = data,
                                                integrated_process = integrated_process, 
                                                observation_model = observation_model)
  }  
  
  if (process_link == 'size_abundance') {
    data_module <- define_size_abundance_module(data = data,
                                                integrated_process = integrated_process, 
                                                observation_model = observation_model)
  }   
  
  if (process_link == 'biomass') {
    data_module <- define_biomass_module(data = data,
                                         integrated_process = integrated_process, 
                                         observation_model = observation_model)
  }   
  
  if (process_link == 'community') {
    data_module <- define_community_module(data = data,
                                           integrated_process = integrated_process, 
                                           observation_model = observation_model)
  }   
  
  data_module <- list(data_module = data_module,
                      data = data,
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
define_growth_module <- function (data, integrated_process, observation_model) {
  
  size_data <- vector('list', length = length(data))
  for (i in seq_along(data)) {
    size_data[[i]] <- lapply(seq_len(ncol(data[[i]])),
                             function(index) greta::as_data(matrix(data[[i]][, index],
                                                                   ncol = ncol(data[[i]]))))
  } 
  
  size_data
  
} 

# internal function: build abundance data module
define_abundance_module <- function (data, integrated_process, observation_model) {
  
  # create output lists
  mu_iterated <- vector("list", length = length(data))
  
  # use separate process models if they exist
  if (integrated_process$replicates > 1) {
    
    for (i in seq_along(integrated_process$replicate_id)) {
      mu_iterated[[i]] <- iterate_state((sweep(integrated_process$parameters$survival[[integrated_process$replicate_id[i]]],
                                               2, integrated_process$parameters$survival_vec[[integrated_process$replicate_id[i]]],
                                               '*') +
                                           integrated_process$parameters$fecundity[[integrated_process$replicate_id[i]]]),
                                        integrated_process$mu_initial[[integrated_process$replicate_id[i]]],
                                        integrated_process$parameters$density_parameter,
                                        seq_len(ncol(data[[i]])),
                                        integrated_process$density_dependence)
      
    } 
  } else {
    
    # fit all elements of data to the same process model
    for (i in seq_len(length(data))) {
      mu_iterated[[i]] <- iterate_state((sweep(integrated_process$parameters$survival[[1]],
                                               2, integrated_process$parameters$survival_vec[[1]],
                                               '*') +
                                           integrated_process$parameters$fecundity[[1]]),
                                        integrated_process$mu_initial[[1]],
                                        integrated_process$parameters$density_parameter,
                                        seq_len(ncol(data[[i]])),
                                        integrated_process$density_dependence)
    
    } 
  }
  
  mu_flattened <- do.call('c', mu_iterated)
  
  mu_flattened
  
} 

# internal function: build mark-recapture data module
define_mark_recapture_module <- function (data, integrated_process, observation_model) {
  
  NULL
  
}

# internal function: build size-abundance data module
define_size_abundance_module <- function (data, integrated_process, observation_model) {
  
  data_module <- NULL
  
  data_module
  
}

# internal function: build biomass data module
define_biomass_module <- function (data, integrated_process, observation_model) {
  
  data_module <- NULL
  
  data_module
  
}

# internal function: build community data module
define_community_module <- function (data, integrated_process, observation_model) {
  
  data_module <- NULL
  
  data_module
  
}

# internal function: convert growth data in long format to matrix used by multinomial
make_growth_data_matrix <- function(data, classes) {
  
  nbreaks <- classes + 1
  
  # calculate size_now and size_next
  size_now <- NULL
  size_next <- NULL
  id <- NULL
  sizemax <- max(tapply(data$growth, data$id, sum, na.rm = TRUE),
                 na.rm = TRUE)
  for (i in seq_along(unique(data$id))) {
    data_sub <- data[data$id == unique(data$id)[i], ]
    size_tmp <- cumsum(data_sub$growth[order(data_sub$year)])
    for (j in seq_len((length(size_tmp) - 1))) {
      size_now <- c(size_now, size_tmp[j])
      size_next <- c(size_next, size_tmp[j + 1])
      id <- c(id, unique(data_sub$id)[i])
    }
  }
  data_clean <- data.frame(id = id,
                           size_now = (size_now / sizemax),
                           size_next = (size_next / sizemax)) 
  
  # calculate breaks
  break_set <- c(0, quantile(data_clean$size_now,
                             p = seq(0.1, 0.9, length = (nbreaks - 2))), 1)
  
  label_set <- seq_len(length(break_set) - 1)
  data_clean$bin_now <- as.numeric(cut(data_clean$size_now,
                                       breaks = break_set,
                                       labels = label_set))
  data_clean$bin_next <- as.numeric(cut(data_clean$size_next,
                                        breaks = break_set,
                                        labels = label_set)) 
  
  # calculate transition probabilities  
  growth_matrix <- matrix(0, nrow = (nbreaks - 1), ncol = (nbreaks - 1))
  for (i in seq_len(nrow(data_clean))) {
    xind <- data_clean$bin_next[i]
    yind <- data_clean$bin_now[i]
    growth_matrix[xind, yind] <- growth_matrix[xind, yind] + 1
  }
  
  growth_matrix
  
}

# internal function: calculate structured capture history from long-format data
calculate_capture_history <- function(data, classes) {
  
  nbreaks <- classes + 1
  
  # remove NAs in fish IDs
  if (any(is.na(data$id))) {
    data <- data[!is.na(data$id), ]
  }
  
  # work out how many times each tagged fish was captured
  recaptures <- tapply(rep(1, nrow(data)), data$id, sum)
  
  # filter to fish that were recaptured at least once
  recaptures <- recaptures[recaptures > 1]
  
  # pull out the sample dates (years at this stage; go to season/month perhaps)
  times <- sort(unique(data$time))
  
  # prepare an output matrix with one row for each fish captured more than once
  capture_history <- matrix(NA, nrow = length(recaptures), ncol = length(times))
  size_history <- matrix(NA, nrow = length(recaptures), ncol = length(times))
  
  # add sample dates and fishids to output matrix
  colnames(capture_history) <- colnames(size_history) <- times
  rownames(capture_history) <- rownames(size_history) <- names(recaptures)
  
  # for each fish, work out which years it was caught
  for (i in seq_along(recaptures)) {
    
    # subset data to a single fish
    data_sub <- data[data$id == names(recaptures)[i], ]
    
    # calculate capture history
    capture_history[i, ] <- ifelse(times %in% data_sub$time, 1, 0)
    
    # calculate size at each recapture
    size_tmp <- tapply(data_sub$size, data_sub$time, mean)
    size_history[i, times %in% names(time_tmp)] <- size_tmp 
    
  }
  
  # calculate breaks
  break_set <- c(0, quantile(data$size,
                             p = seq(0.1, 0.9, length = (nbreaks - 2))), 1)
  
  
}

