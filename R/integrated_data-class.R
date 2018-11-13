#' @name integrated_data
#' @title create integrated data objects
#' 
#' @description \code{integrated_data} defines a data object with appropriate likelihood based on
#'  a process model defined with \link[integrated]{integrated_process}
#' 
#' @param data something
#' @param process something
#' @param bias something else
#' @param settings a named list of settings passed to data formatting functions (see details)
#' @param ... additional arguments
#' @param x
#' @param object
#'
#' @details Do something. The settings list can be used to specify how the data are binned, either with
#'   specific breaks for binning or with the number of breaks to use. If these are not provided, the 
#'   functions use the \code{classes} argument in the \code{integrated_process} to determine the number
#'   of bins (\code{nbreaks = classes + 1}).
#' 
#' @return An object of class \code{integrated_data}, which contains information on the data module and
#'   can be passed to \link[integrated]{integrated_model}
#' 
#' @import greta
#' 
#' @examples
#' \dontrun{
#' 
#' library(integrated)
#' 
#' # prepare an example model
#' data <- add_data()
#'                         
#' # summarise data module
#' model
#' summary(model)
#' plot(model)
#' }

#' @export
#' @rdname integrated_data
#' 
abundance <- function(data,
                      process,
                      bias,
                      settings = list()) {
  
  # prepare abund matrix
  
  # match dims and classes to work out if it's pop or class

  # use process to determine whether it's stage or age (warn)
  
  data_module <- abundance_module(data = data,
                                  process = process, 
                                  bias = bias)
  
  # return outputs
  as.integrated_data(data_module)
  
}

#' @export
#' @rdname integrated_data
#' 
cmr <- function(data,
                process,
                bias,
                settings = list()) {
  
  # prepare cmr histories
  
  # return outputs
  as.integrated_data(data_module)
  
} 


#' @export
#' @rdname integrated_data
#' 
growth <- function(data,
                   process,
                   bias,
                   settings = list()) {
  
  # prepare cmr histories
  
  # return outputs
  as.integrated_data(data_module)
  
}   

#' @export
#' @rdname integrated_data
#' 
community <- function(data,
                      process,
                      bias,
                      settings = list()) {
  
  # prepare cmr histories
  
  # return outputs
  as.integrated_data(data_module)
  
}  

function () {

  if (!(process_link %in% c("individual_growth",
                            "age_abundance",
                            "stage_abundance",
                            "age_recapture",
                            "stage_recapture",
                            "population_abundance",
                            "population_biomass",
                            "community"))) {
    stop("process_link must be a known process module",
         call. = FALSE)
  }  
  
  if (process_link == "individual_growth") {
    
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
      stop("growth data must be a matrix, data.frame, or list of matrices or data.frames",
           call. = FALSE)
    }
    
    # error if the number of data classes doesn't match the process
    if (!all(integrated_process$classes == sapply(data, nrow))) {
      stop("number of classes in growth data must match number of classes in integrated_process")
    }
    
    if (integrated_process$replicates > 1) {      
      if (length(data) > 1) {
        if (length(data) != length(integrated_process$replicate_id)) {
          stop(paste0("growth data have ", length(data), " elements but ",
                      " integrated process has ", integrated_process$replicates,
                      " replicates"),
               call. = FALSE)
        }
      } else {
        stop(paste0("one growth data set should be supplied for each of the ",
                    integrated_process$replicates, " process replicates"),
             call. = FALSE)
      }
    }
    
    # create data module from growth data matrix
    data_module <- define_individual_growth_module(data = data,
                                                   integrated_process = integrated_process, 
                                                   observation_model = observation_model)
    
  }    
  
  if (process_link == "age_abundance") {
    
    # create data module from age abundance data matrix
    data_module <- define_age_abundance_module(data = data,
                                               integrated_process = integrated_process, 
                                               observation_model = observation_model)
    
  }
  
  if (process_link == "stage_abundance") {
    
    # check data format
    if (is.matrix(data) | is.data.frame(data)) {
      
      # want the data to be a matrix with a size column or a matrix with classes in rows and samples in columns
      if (!("size" %in% colnames(data))) {
        
        data <- list(data)
        
      } else {
        
        if (!all(c("size", "time", "site") %in% colnames(pop_data))) {
          stop("data with a size column must also have site and time columns",
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
      stop("abundance data must contain a size column or be a list with classes in rows and samples in columns",
           call. = FALSE) 
    }
    
    # if there is more than one replicate
    if (integrated_process$replicates > 1) {
      # check that there is one data element for each replicate
      if (length(data) != length(integrated_process$replicate_id)) {
        stop(paste0("abundance data have ", length(data), " elements but ",
                    " integrated_process contains ", integrated_process$replicates,
                    " replicates"),
             call. = FALSE)
      }
    }
    
    # this won't work if there are more classes in data than in the process model
    if (max(sapply(data, nrow)) > integrated_process$classes) {
      stop(paste0("there are up to ", max(sapply(data, nrow)), " classes in the data set ",
                  "but only ", integrated_process$classes, " classes in integrated_process"),
           call. = FALSE)
    }
    
    # create data module from list data
    data_module <- define_stage_abundance_module(data = data,
                                                 integrated_process = integrated_process, 
                                                 observation_model = observation_model)
  }  
  
  if (process_link == "age_recapture") {
    
    data_module <- define_age_recapture_module(data = data,
                                               integrated_process = integrated_process,
                                               observation_model = observation_model)
    
  }
  
  if (process_link == "stage_recapture") {
    
    # turn data into a list and check that each element has the correct columns
    if (is.matrix(data) | is.data.frame(data)) {
      
      data <- list(data)
      
    } else {
      
      if (!is.list(data)) {
        stop("stage_recapture data must be a matrix, data.frame, or list of matrices or data.frames",
             call. = FALSE) 
      }
      
    }
    
    # check data format
    for (i in seq_along(data)) {
      if (!("size" %in% colnames(data[[i]]))) {
        stop("stage_recapture models require size measurements at each recapture",
             call. = FALSE)
      }
      
      if (!all(c("size", "id", "time") %in% colnames(data[[i]]))) { 
        stop("stage_recapture data should be in long format with size, id, and time columns",
             call. = FALSE)
      }
    }
    
    # create capture histories (binary and structured)
    data <- calculate_capture_history(data = data,
                                      classes = integrated_process$classes,
                                      settings = settings)
    
    # create data module
    data_module <- define_stage_recapture_module(data = data,
                                                 integrated_process = integrated_process, 
                                                 observation_model = observation_model)
    
  }    
  
  if (process_link == "population_abundance") {
    
    data_module <- define_population_abundance_module(data = data,
                                                      integrated_process = integrated_process, 
                                                      observation_model = observation_model)
    
  }     
  
  if (process_link == "population_biomass") {
    
    data_module <- define_population_biomass_module(data = data,
                                                    integrated_process = integrated_process, 
                                                    observation_model = observation_model)
    
  }     
  
  if (process_link == "community") {
    
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

#' @export
#' @rdname integrated_data
#' 
is.integrated_data <- function(x) {
  inherits(model, integrated_data)
}

#' @export
#' @rdname integrated_data
#' 
print.integrated_data <- function(x, ...) {
  cat(paste0(This is an integrated_data object\n))
}

#' @export
#' @rdname integrated_data
#' 
plot.integrated_data <- function(x, ...) {
  
  plot(x$greta_model, ...)
  
}

#' @export
#' @rdname integrated_data
#' 
summary.integrated_data <- function(object, ...) {
  
  NULL
  
}


# internal function: create integrated_data object
as.integrated_data <- function(model) {
  as_class(model, name = integrated_data, type = list)
}

# internal function: build growth data module
define_individual_growth_module <- function (data, integrated_process, observation_model) {
  
  size_data <- vector(list, length = length(data))
  for (i in seq_along(data)) {
    size_data[[i]] <- sapply(seq_len(ncol(data[[i]])),
                             function(index) greta::as_data(matrix(data[[i]][, index],
                                                                   ncol = ncol(data[[i]]))))
  } 
  
  size_data
  
} 

# internal function: build age abundance data module
define_age_abundance_module <- function (data, integrated_process, observation_model) {
  
  data_module <- NULL
  
  data_module
  
}

# internal function: build stage abundance data module
define_stage_abundance_module <- function (data, integrated_process, observation_model) {
  
  # create output lists
  mu_iterated <- vector("list", length = length(data))
  
  # use separate process models if they exist
  if (integrated_process$replicates > 1) {
    
    for (i in seq_along(integrated_process$replicate_id)) {
      mu_iterated[[i]] <- iterate_state(t(greta::sweep(integrated_process$parameters$survival[[integrated_process$replicate_id[i]]],
                                                       2, integrated_process$parameters$survival_vec[[integrated_process$replicate_id[i]]],
                                                       "*") +
                                            integrated_process$parameters$fecundity[[integrated_process$replicate_id[i]]]),
                                        integrated_process$mu_initial[[integrated_process$replicate_id[i]]],
                                        integrated_process$parameters$density_parameter[[integrated_process$replicate_id[i]]],
                                        seq_len(ncol(data[[i]])),
                                        integrated_process$density_dependence)
      
    } 
  } else {  
    
    # fit all elements of data to the same process model
    for (i in seq_len(length(data))) {
      mu_iterated[[i]] <- iterate_state(t(greta::sweep(integrated_process$parameters$survival[[1]],
                                                       2, integrated_process$parameters$survival_vec[[1]],
                                                       "*") +  
                                            integrated_process$parameters$fecundity[[1]]),
                                        integrated_process$mu_initial[[1]],
                                        integrated_process$parameters$density_parameter[[1]],
                                        seq_len(ncol(data[[i]])),
                                        integrated_process$density_dependence)
    
    } 
  }
  
  mu_flattened <- do.call("c", mu_iterated)
  
  mu_flattened
  
} 

# internal function: build age mark-recapture data module
define_age_recapture_module <- function (data, integrated_process, observation_model) {
  
  data_module <- NULL
  
  data_module
  
}

# internal function: build stage mark-recapture data module
define_stage_recapture_module <- function (data, integrated_process, observation_model) {
  
  history <- vector("list", length = length(data))
  unique_history <- vector("list", length = length(data))
  count <- vector("list", length = length(data))
  count2 <- vector("list", length = length(data))
  total <- vector("list", length = length(data))

  for (i in seq_along(data)) {
    
    # reduce capture histories to a list without pre/post capture information
    history[[i]] <- vector("list", length = nrow(data[[i]]$structured))
    ntime <- ncol(data[[i]]$structured)
    for (j in seq_along(history[[i]])) {
      data_tmp <- data[[i]]$structured[j, which.max(data[[i]]$binary[j, ]):(ntime - which.max(rev(data[[i]]$binary[j, ])) + 1)]
      history[[i]][[j]] <- data_tmp
    }
    
    # calculate unique CMR histories and counts of each
    unique_history_vec <- unique(history[[i]])
    unique_history[[i]] <- vector("list", length = length(unique_history_vec))
    count[[i]] <- matrix(0, nrow = integrated_process$classes, ncol = integrated_process$classes)
    for (j in seq_along(unique_history_vec)) {
      mat_tmp <- c(t(matrix(0, nrow = length(unique_history_vec[[j]]), ncol = integrated_process$classes)))
      mat_tmp[seq(1, length(unique_history_vec[[j]]) * integrated_process$classes,
                  by = integrated_process$classes)[seq_len(length(unique_history_vec[[j]]))] +
                ifelse(unique_history_vec[[j]] == 0, 1, unique_history_vec[[j]]) - 1] <- unique_history_vec[[j]]
      unique_history[[i]][[j]] <- matrix(ifelse(mat_tmp > 0, 1, 0), ncol = integrated_process$classes,
                                         byrow = TRUE)
      
      for (k in seq_len(nrow(unique_history[[i]][[j]]))[-1]) {
        ind1 <- which(unique_history[[i]][[j]][(k - 1), ] != 0)
        ind2 <- which(unique_history[[i]][[j]][k, ] != 0)
        if (length(ind1) & length(ind2)) {
          count[[i]][ind1, ind2] <- count[[i]][ind1, ind2] + 1
        }
        if (length(ind1) & !(length(ind2))) {
          if (ind1 < (integrated_process$classes - 1)) {
            ind1_set <- ind1:(ind1 + 2)
          } else {
            if (ind1 < integrated_process$classes) {
              ind1_set <- ind1:(ind1 + 1)
            } else {
              ind1_set <- ind1
            }
          }
          count[[i]][ind1, ind1_set] <- count[[i]][ind1, ind1_set] + 1
        }
        if (!(length(ind1)) & length(ind2)) {
          if (ind2 < (integrated_process$classes - 1)) {
            ind2_set <- (ind2 - 2):ind2
          } else {
            if (ind2 < integrated_process$classes) {
              ind2_set <- (ind2 - 1):ind2
            } else {
              ind2_set <- ind2
            }
          }
          for (kk in seq_along(ind2_set)) {
            count[[i]][ind2_set[kk], ind2] <- count[[i]][ind2_set[kk], ind2] + 1
          }
        }
      }
      
    }

    sizes_lived <- apply(data[[i]]$structured, 1, count_stages_lived, classes = integrated_process$classes)
    total[[i]] <- apply(sizes_lived, 1, sum)
    sizes_survived <- apply(data[[i]]$structured, 1, count_stages_survived, classes = integrated_process$classes)
    count2[[i]] <- apply(sizes_survived, 1, sum)
    total[[i]] <- ifelse(count2[[i]] == 0, total[[i]] + 1, total[[i]])
    count2[[i]] <- ifelse(count2[[i]] == 0, count2[[i]] + 1, count2[[i]])
      
  } 
  
  # collate outputs
  cmr_module <- list(history = unique_history,
                     count = count,
                     count2 = count2,
                     total = total)

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
