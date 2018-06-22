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
    stop('process_link must be a known process module')
  }  
  
  if (process_link == 'growth') {
    
    # check if data are a list or just a single matrix
    if (is.matrix(data) | is.data.frame(data)) {
      
      # check if data are formatted correctly
      if (ncol(data) != nrow(data)) {
        data <- make_growth_data_matrix(data = data,
                                        classes = integrated_process$classes)
      }
      
      data <- list(data)
      
    } else {

      if (is.list(data)) {
      if ((length(data) > 1) & (integrated_process$replicates > 1)) {      
        if (length(data) != integrated_process$replicates) {
          stop('growth data with more than one element should have one matrix for each replicate')
        }
      }
      } else {
        stop('growth data must be a matrix, data.frame, or list of matrices or data.frames')
      }
      
    }
    
    # create data module from growth data matrix
    data_module <- define_growth_module(data = data,
                                        integrated_process = integrated_process, 
                                        observation_model = observation_model)
  }   
  
  if (process_link == 'abundance') {
    
    # check data format
    if (!is.list(data)) {
      if (is.matrix(data) | is.data.frame(data)) {
        data <- list(data)
      } else {
        stop('abundance data must be a list with classes in rows and samples in columns')
      }
    }
    
    # if there is more than one replicate
    if (integrated_process$replicates > 1) {
      # check that there is one data element for each replicate
      if (length(data) != integrated_process$replicates) {
        stop('abundance data should be a list with one entry for each replicate in integrated_process')
      }
    }
    
    # this won't work if there are more classes in data than in the process model
    if (max(sapply(data, nrow)) > integrated_process$classes) {
      stop(paste0('there are ', max(sapply(data, nrow)), ' classes in the data set ',
                  'but only ', integrated_process$classes, ' classes in integrated_process'))
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
  mu_iterated <- vector("list", length = integrated_process$replicates)
  
  # use separate process models if they exist
  if (integrated_process$replicates > 1) {
    
    for (i in seq_len(integrated_process$replicates)) {
      mu_iterated[[i]] <- iterate_state((integrated_process$parameters$survival[[i]] +
                                           integrated_process$parameters$fecundity[[i]]),
                                        integrated_process$mu_initial[[i]],
                                        seq_len(ncol(data[[i]])))
      
    } 
  } else {
    
    # fit all elements of data to the same process model
    for (i in seq_len(length(data))) {
      mu_iterated[[i]] <- iterate_state((integrated_process$parameters$survival[[1]] +
                                           integrated_process$parameters$fecundity[[1]]),
                                        integrated_process$mu_initial[[1]],
                                        seq_len(ncol(data[[i]])))
    
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
      id <- c(id, unique(data_sub$id)[1])
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

