# internal helper functions for integrated R package

# count stages that an individual successfully transitions out of
count_stages_survived <- function(x, classes) {
  out <- rep(0, classes)
  xpos <- x[x > 0]
  xsurv <- xpos[xpos < max(xpos)]
  out[seq_len(classes) %in% xsurv] <- 1
  tmp_var <- runif(classes)
  out[tmp_var < seq(0.3, 1, length = classes)] <-
    ifelse(seq_len(classes) %in% xpos, 1, 0)[tmp_var < seq(0.3, 1, length = classes)]
  out
}

# count stages that an individual is observed in
count_stages_lived <- function(x, classes) {
  out <- rep(0, classes)
  x <- x[x > 0]
  out[seq_len(classes) %in% x] <- 1
  out
}

# calculate stage at last observation
calculate_final_stage <- function(x) {
  x[max(which(x != 0))]
}

# set an object class
as_class <- function (object, name, type = c("function", "list")) {
  
  type <- match.arg(type)
  stopifnot(inherits(object, type))
  class(object) <- c(name, class(object))
  
  object
  
}

# test if values are equal to zero with a non-zero tolerance
almost_equal <- function(x, y, tolerance = 1e-10) {
  
  diff <- abs(x - y)
  mag <- pmax(abs(x), abs(y))
  
  ifelse(mag > tolerance, (diff / mag) <= tolerance, diff <= tolerance)
  
}

# test if multiple values are all equal
all_equal <- function (..., tolerance = 1e-10) {
  
  test <- list(...)
  
  if (length(test) > 1) {
    out <- rep(NA, length(test))
    for (i in seq_along(test)) {
      out[i] <- almost_equal(test[[i]], test[[1]], tolerance = tolerance)
    }
  } else {
    out <- TRUE
  }
  
  all(out)
  
}
