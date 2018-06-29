#' @import tensorflow

# define op function to setup tensorflow functions
op <- greta::.internals$nodes$constructors$op

# internal function: create greta_array containing probabilities of CMR histories
calculate_history_probability <- function(history, capture_probability, parameters) {
  
  # index used to ID recaptured and not-recaptured individuals
  nrows <- sapply(history, nrow)
  nstage <- ncol(parameters)
  
  # collapse the list of matrices into one big matrix  
  history_mat <- do.call('rbind', history)
  
  # separate the final observation of each individual from earlier observations
  end_index <- cumsum(nrows)
  start <- history_mat[-end_index, ]
  end <- history_mat[end_index, ]
  
  # create unique index for each stage and each individual
  stage_id <- rep(seq_len(nstage), times = nrow(start)) +
    rep(seq(from = 0, to = (nstage * sum(nrows > 1) - 1), by = nstage),
        times = (nstage * (nrows[nrows > 1] - 1)))
  stage_id <- as.integer(stage_id) - 1L
  
  # for each individual, work out when it was observed and unobserved
  obs_tmp <- apply(start, 1, function(x) any(x != 0))
  observed <- which(obs_tmp) - 1L
  unobserved <- which(!obs_tmp) - 1L
  
  # create indices to sort and reorder intermediate and final outputs
  sort <- match(seq_along(obs_tmp), c(which(obs_tmp), which(!obs_tmp))) - 1L
  final <- match(seq_along(nrows), c(which(nrows > 1), which(nrows == 1))) - 1L
  
  # identify individuals observed once and those observed multiple times
  single <- which(nrows == 1) - 1L
  multi <- which(nrows > 1) - 1L
  
  dimfun <- function (x) c(length(history), 1)
  
  greta:::op('calculate_history_probability',
             capture_probability, parameters,
             operation_args = list(start = start,
                                   end = end,
                                   nstage = nstage,
                                   sort = sort,
                                   stage_id = stage_id,
                                   observed = observed,
                                   unobserved = unobserved,
                                   single = single,
                                   multi = multi,
                                   final = final),
             tf_operation = 'integrated:::tf_calculate_history_probability',
             dimfun = dimfun)
  
}

# internal: tf function to calculate cmr history probabilities
tf_calculate_history_probability <- function(capture_probability,
                                             parameters,
                                             start, end,
                                             nstage, sort, stage_id,
                                             observed, unobserved,
                                             single, multi, final) {
  
  # calculate probabilities of all possible states at times t>1 given state at time 1
  state <- tf$matmul(parameters, tf$transpose(tf$constant(start, dtype = tf$float32)))
  
  # multiply transition probabilities by probability of being captured in a given state
  state_observed <- do.call('*', list(tf$transpose(tf$gather(tf$transpose(state), observed)),
                                      capture_probability)) 
  state_unobserved <- do.call('*', list(tf$transpose(tf$gather(tf$transpose(state), unobserved)),
                                        1.0 - capture_probability)) 
  states <- tf$concat(list(state_observed, state_unobserved), axis = 1L) 
  
  # reorder vector so that observed and unobserved elements are in the correct
  #   places rather than just concatenated
  states <- tf$transpose(tf$gather(tf$transpose(states), sort))
  
  # calculate products of each stage within each individual's capture history
  stage_prob <- tf$unsorted_segment_prod(
    tf$reshape(tf$transpose(states), shape = c(length(states), 1L)),
    stage_id, length(unique(stage_id))) 
  stage_prob <- tf$transpose(tf$reshape(stage_prob, shape = c(nstage, length(multi))))
  
  # separate vectors and matrices for individuals that were and were not recaptured
  single_prob <- do.call('*',
                         list(tf$transpose(tf$gather(tf$constant(end, dtype = tf$float32), single)),
                              capture_probability)) 
  multi_prob <- stage_prob * tf$gather(tf$constant(end, dtype = tf$float32), multi)
  
  # calculate probs of recaptures/detection for all individuals (out of order)  
  probs_tmp <- tf$concat(list(tf$reduce_sum(tf$transpose(multi_prob), axis = 0L),
                              tf$reduce_sum(single_prob, axis = 0L)),
                         axis = 0L)
  
  # reorder so that the single observations get reinserted at correct point
  probs <- tf$gather(probs_tmp, final)
  
  probs
  
}
