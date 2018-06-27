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
#' @param settings a named list of settings passed to data formatting functions (see details)
#' @param ... additional arguments
#'
#' @return An object of class \code{integrated_data}, which contains information on the data module and
#'   can be passed to a \link[integrated]{integrated_model}
#' 
#' @details Do something. The settings list can be used to specify how the data are binned, either with
#'   specific breaks for binning or with the number of breaks to use. If these are not provided, the 
#'   functions use the \code{classes} argument in the \code{integrated_process} to determine the number
#'   of bins (\code{nbreaks = classes + 1}).
#' 
#' @export
#' 
#' @import greta
#' @import tensorflow
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
                                    observation_model = 'naive',
                                    settings = list()) {
  
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
                                        classes = integrated_process$classes,
                                        settings = settings)
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
      
      # want the data to be a matrix with a size column or a matrix with classes in rows and samples in columns
      if (!('size' %in% colnames(data))) {
        
        data <- list(data)
        
      } else {
        
        if (!all(c('size', 'time', 'site') %in% colnames(pop_data))) {
          stop('data with a size column must also have site and time columns',
               call. = FALSE)
        }
        
        # need to collapse data into appopriate size classes
        data <- make_pop_data_matrix(data = data,
                                     classes = integrated_process$classes,
                                     settings = settings)
        
      }
    }
    
    # error if not a list or matrix
    if (!is.list(data)) {
      stop('abundance data must contain a size column or be a list with classes in rows and samples in columns',
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
    
    # turn data into a list and check that each element has the correct columns
    if (is.matrix(data) | is.data.frame(data)) {
      
      data <- list(data)
      
    } else {
      
      if (!is.list(data)) {
        stop('mark_recapture data must be a matrix, data.frame, or list of matrices or data.frames',
             call. = FALSE) 
      }
      
    }

    # check data format
    for (i in seq_along(data)) {
      if (!('size' %in% colnames(data[[i]]))) {
        stop('mark_recapture models require size measurements at each recapture',
             call. = FALSE)
      }
      
      if (!all(c('size', 'id', 'time') %in% colnames(data[[i]]))) { 
        stop('mark_recapture data should be in long format with size, id, and time columns',
             call. = FALSE)
      }
    }
    
    # create capture histories (binary and structured)
    data <- calculate_capture_history(data = data,
                                      classes = integrated_process$classes,
                                      settings = settings)
    
    # create data module
    data_module <- define_mark_recapture_module(data = data,
                                                integrated_process = integrated_process, 
                                                observation_model = observation_model)
    
  }  
  
  if (process_link == 'population_abundance') {
    data_module <- define_population_abundance_module(data = data,
                                                      integrated_process = integrated_process, 
                                                      observation_model = observation_model)
  }    
  
  if (process_link == 'population_biomass') {
    data_module <- define_population_biomass_module(data = data,
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
                      observation_model = observation_model,
                      settings = settings)
  
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
  
  history <- vector('list', length = length(data))
  unique_history <- vector('list', length = length(data))
  count <- vector('list', length = length(data))
  
  for (i in seq_along(data)) {
    
    # reduce capture histories to a list without pre/post capture information
    history[[i]] <- vector('list', length = nrow(data[[i]]$structured))
    ntime <- ncol(data[[i]]$structured)
    for (j in seq_along(history[[i]])) {
      data_tmp <- data[[i]]$structured[j, which.max(data[[i]]$binary[j, ]):(ntime - which.max(rev(data[[i]]$binary[j, ])) + 1)]
      history[[i]][[j]] <- data_tmp
    }
    
    # calculate unique CMR histories and counts of each
    unique_history_vec <- unique(history[[i]])
    unique_history[[i]] <- vector('list', length = length(unique_history_vec))
    count[[i]] <- rep(NA, length = length(unique_history_vec))
    for (j in seq_along(unique_history_vec)) {
      count[j] <- sum(sapply(history[[i]], function(x) ifelse(length(x) == length(unique_history_vec[[j]]),
                                                              all(x == unique_history_vec[[j]]),
                                                              FALSE))) 
      mat_tmp <- c(t(matrix(0, nrow = length(unique_history_vec[[j]]), ncol = integrated_process$classes)))
      mat_tmp[seq(1, length(unique_history_vec[[j]]) * integrated_process$classes,
                  by = integrated_process$classes)[seq_len(length(unique_history_vec[[j]]))] +
                ifelse(unique_history_vec[[j]] == 0, 1, unique_history_vec[[j]]) - 1] <- unique_history_vec[[j]]
      unique_history[[i]][[j]] <- matrix(ifelse(mat_tmp > 0, 1, 0), ncol = integrated_process$classes,
                                         byrow = TRUE)
    }   
    count[[i]] <- matrix(count[[i]], nrow = 1)
    
  } 
  
  # initialise main outputs
  probs <- vector('list', length = length(data))
  
  # use separate process models if they exist
  if (integrated_process$replicates > 1) {
    
    capture_probability <- vector('list', length = length(data))
    for (i in seq_along(integrated_process$replicate_id)) {
      
      # define observation matrix
      capture_probability[[i]] <- greta::uniform(min = 0.5, max = 0.9,
                                                 dim = c(1, integrated_process$classes))
      probs[[i]] <- calculate_history_probability(history = unique_history[[i]],
                                                  capture_probability = capture_probability[[i]],
                                                  parameters = greta::sweep(integrated_process$parameters$survival[[integrated_process$replicate_id[i]]],
                                                                            2, integrated_process$parameters$survival_vec[[integrated_process$replicate_id[i]]], '*'))
      
    }
    
  } else {
    
    # just a single set of parameters
    capture_probability <- list(greta::uniform(min = 0.5, max = 0.9,
                                               dim = c(integrated_process$classes, 1)))
    
    # fit all elements of data to the same process model
    for (i in seq_along(data)) {
      
      probs[[i]] <- calculate_history_probability(history = unique_history[[i]],
                                                  capture_probability = capture_probability[[1]],
                                                  parameters = greta::sweep(integrated_process$parameters$survival[[1]],
                                                                            2, integrated_process$parameters$survival_vec[[1]], '*'))
      
    }
    
  }
  
  # collate outputs
  cmr_module <- list(history = history,
                     unique = unique_history,
                     capture_probability = capture_probability,
                     count = count,
                     probs = probs)

  cmr_module
  
}

# internal function: build population abundance data module
define_population_abundance_module <- function (data, integrated_process, observation_model) {
  
  data_module <- NULL
  
  data_module
  
}

# internal function: build population biomass data module
define_population_biomass_module <- function (data, integrated_process, observation_model) {
  
  data_module <- NULL
  
  data_module
  
}

# internal function: build community data module
define_community_module <- function (data, integrated_process, observation_model) {
  
  data_module <- NULL
  
  data_module
  
}

# internal function: convert growth data in long format to matrix used by multinomial
make_growth_data_matrix <- function(data, classes, settings) {
  
  # unpack settings
  matrix_set <- list(nbreaks = NULL,
                     breaks = NULL)
  matrix_set[names(settings)] <- settings
  
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
  if (is.null(matrix_set$breaks)) {
    if (is.null(matrix_set$nbreaks)) {
      matrix_set$nbreaks <- classes + 1
    }
    break_set <- c(0, quantile(data_clean$size_now,
                               p = seq(0.1, 0.9, length = (matrix_set$nbreaks - 2))), 1)
  } else {
    break_set <- matrix_set$breaks
    if (is.null(matrix_set$nbreaks)) {
      matrix_set$nbreaks <- length(break_set)
    }
    if (matrix_set$nbreaks != length(break_set)) {
      stop(paste0(length(break_set), ' breaks have been specified but nbreaks = ',
                  matrix_set$nbreaks),
           call. = FALSE)
    }
  }
  
  label_set <- seq_len(length(break_set) - 1)
  data_clean$bin_now <- as.numeric(cut(data_clean$size_now,
                                       breaks = break_set,
                                       labels = label_set))
  data_clean$bin_next <- as.numeric(cut(data_clean$size_next,
                                        breaks = break_set,
                                        labels = label_set)) 
  
  # calculate transition probabilities  
  growth_matrix <- matrix(0, nrow = (matrix_set$nbreaks - 1), ncol = (matrix_set$nbreaks - 1))
  for (i in seq_len(nrow(data_clean))) {
    xind <- data_clean$bin_next[i]
    yind <- data_clean$bin_now[i]
    growth_matrix[xind, yind] <- growth_matrix[xind, yind] + 1
  }
  
  growth_matrix
  
}

# internal function: convert population data in long format to a list of structured matrices
make_pop_data_matrix <- function(data, classes, settings) {
  
  # unpack settings
  matrix_set <- list(nbreaks = classes + 1,
                     breaks = NULL)
  matrix_set[names(settings)] <- settings
  
  # calculate breaks
  if (is.null(matrix_set$breaks)) {
    if (is.null(matrix_set$nbreaks)) {
      matrix_set$nbreaks <- classes + 1
    }
    break_set <- c(0, quantile(data$size,
                               p = seq(0.1, 0.9, length = (matrix_set$nbreaks - 2))), 1)
  } else {
    break_set <- matrix_set$breaks
    if (is.null(matrix_set$nbreaks)) {
      matrix_set$nbreaks <- length(break_set)
    }
    if (matrix_set$nbreaks != length(break_set)) {
      stop(paste0(length(break_set), ' breaks have been specified but nbreaks = ',
                  matrix_set$nbreaks),
           call. = FALSE)
    }
  }  
  
  # loop through sites and calculate hist in each time and site
  popdata <- vector('list', length = length(unique(data$site)))
  for (i in seq_along(unique(data$site))) {
    data_sub <- data[(data$site == unique(data$site)[i]), ]
    hist_tmp <- tapply(data$size, data$time, function(x) hist(x, breaks = break_set, plot = FALSE)$counts)
    popdata[[i]] <- matrix(unlist(hist_tmp), nrow = (length(break_set) - 1))
  }
  
  names(popdata) <- unique(data$site)
  
  popdata
  
}

# internal function: calculate structured capture history from long-format data
calculate_capture_history <- function(data, classes, settings) {
  
  # unpack settings
  matrix_set <- list(nbreaks = NULL,
                     breaks = NULL)
  matrix_set[names(settings)] <- settings
  
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
    size_history[i, times %in% names(size_tmp)] <- size_tmp 
    
  }
  
  # calculate breaks
  if (is.null(matrix_set$breaks)) {
    if (is.null(matrix_set$nbreaks)) {
      matrix_set$nbreaks <- classes + 1
    }
    break_set <- c(0, quantile(data$size,
                               p = seq(0.1, 0.9, length = (matrix_set$nbreaks - 2))), 1)
  } else {
    break_set <- matrix_set$breaks
    if (is.null(matrix_set$nbreaks)) {
      matrix_set$nbreaks <- length(break_set)
    }
    if (matrix_set$nbreaks != length(break_set)) {
      stop(paste0(length(break_set), ' breaks have been specified but nbreaks = ',
                  matrix_set$nbreaks),
           call. = FALSE)
    }
  }
  
  # put the size history into categories based on breaks
  size_history_binned <- matrix(cut(size_history, breaks = break_set, labels = FALSE),
                                ncol = ncol(capture_history))
  size_history_binned <- ifelse(is.na(size_history_binned), 0, size_history_binned)
  
  # collate outputs
  capture_history <- list(binary = capture_history,
                          structured = size_history_binned)
  
  # return outputs
  capture_history
  
}

# internal function: create greta_array containing probabilities of CMR histories
calculate_history_probability <- function(history, capture_probability, parameters) {

  history_mat <- greta::as_data(do.call('rbind', history))

  observed <- apply(history_mat, 1, function(x) any(x != 0))
  
  state_vector <- parameters %*% t(history_mat)
  
  state_vector[, observed] <- greta::sweep(state_vector[, observed],
                                           1, capture_probability, '*')
  state_vector[, !observed] <- greta::sweep(state_vector[, !observed],
                                            1, (1 - capture_probability), '*')
  
  id_vec <- rep(seq_along(history),
                times = sapply(history, function(x) nrow(x) * ncol(x)))
  probs <- greta::tapply(c(state_vector), id_vec, 'prod')
  
}
