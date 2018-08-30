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
  
  #### ADD history of never recaptured (1 if captured at a later point)
  ## apply(structured, 1, function(x) max(which(x != 0)))
  ## THIS gives index, which can be used to create new mat
  ##   and also to identify size at last capture
  ### ALL OF THIS assumes we focus on apparent survival
  ### ideally would have a latent state for "alive but not recaptured"

  ## DO THIS FIRST: others can be post-processed
  #### ADD individuals sampled at one time point only
  
  # unpack settings
  matrix_set <- list(nbreaks = NULL,
                     breaks = NULL)
  matrix_set[names(settings)] <- settings
  
  capture_history <- vector('list', length = length(data))
  for (i in seq_along(data)) {
    
    # remove NAs in fish IDs
    if (any(is.na(data[[i]]$id))) {
      data[[i]] <- data[[i]][!is.na(data[[i]]$id), ]
    }
    
    # work out how many times each tagged fish was captured
    recaptures <- tapply(rep(1, nrow(data[[i]])), data[[i]]$id, sum)
    
    # filter to fish that were recaptured at least once
    recaptures <- recaptures[recaptures > 1]
    
    # pull out the sample dates (years at this stage; go to season/month perhaps)
    times <- sort(unique(data[[i]]$time))
    
    # prepare an output matrix with one row for each fish captured more than once
    capture_history_tmp <- matrix(NA, nrow = length(recaptures), ncol = length(times))
    size_history <- matrix(NA, nrow = length(recaptures), ncol = length(times))
    
    # add sample dates and fishids to output matrix
    colnames(capture_history_tmp) <- colnames(size_history) <- times
    rownames(capture_history_tmp) <- rownames(size_history) <- names(recaptures)
    
    # for each fish, work out which years it was caught
    for (j in seq_along(recaptures)) {
      
      # subset data to a single fish
      data_sub <- data[[i]][data[[i]]$id == names(recaptures)[j], ]
      
      # calculate capture history
      capture_history_tmp[j, ] <- ifelse(times %in% data_sub$time, 1, 0)
      
      # calculate size at each recapture
      size_tmp <- tapply(data_sub$size, data_sub$time, mean)
      size_history[j, times %in% names(size_tmp)] <- size_tmp 
      
    }
    
    # calculate breaks
    if (is.null(matrix_set$breaks)) {
      if (is.null(matrix_set$nbreaks)) {
        matrix_set$nbreaks <- classes + 1
      }
      break_set <- c(0, quantile(data[[i]]$size,
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
                                  ncol = ncol(capture_history_tmp))
    size_history_binned <- ifelse(is.na(size_history_binned), 0, size_history_binned)
    
    # collate outputs
    capture_history[[i]] <- list(binary = capture_history_tmp,
                                 structured = size_history_binned)
    
  }
  
  # return outputs
  capture_history
  
}
