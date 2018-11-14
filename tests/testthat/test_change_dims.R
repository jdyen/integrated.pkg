context("change_dims function")

dist_list <- greta:::distribution_classes_module
dist_names <- names(dist_list)

discrete <- c("bernoulli_distribution",
              "binomial_distribution",
              "negative_binomial_distribution",
              "poisson_distribution",
              "multinomial_distribution",
              "categorical_distribution",
              "dirichlet_multinomial_distribution")

too_hard <- "wishart_distribution"

# test operations as well
valid_ops <- c("ilogit", "iprobit", "icloglog", "icauchit",
               "log1pe", "imultilogit",
               "add", "subtract", "multiply", "divide", "power",
               "log", "exp", "log1p", "expm1",
               "cos", "sin", "tan", "acos", "asin", "atan")


test_that("errors on discrete distributions", {
  
}) 

test_that("errors on wishart distribution", {
}) 

test_that("changing dimensions of probability distributions works", {
  
  
  for (i in seq_along(dist_list)) {
    
    if (!(dist_names[i] %in% c(not_allowed, too_hard))) {
      
      dist_tmp <- strsplit(dist_names[i], "_")[[1]]
      dist_tmp <- paste(dist_tmp[seq_len(length(dist_tmp) - 1)], collapse = "_")
      
      arg_list <- names(formals(dist_list[[i]]$public_methods$initialize))
      
      nargs <- sum(!(arg_list %in% c("dim", "truncation", "n_realisations", "dimension")))
      
      param_list <- replicate(nargs, 0, simplify = FALSE)
      
      names(param_list) <- arg_list[seq_len(nargs)]
      
      if (dist_names[i] == "uniform_distribution")
        param_list$max <- 1
      
      new_dim <- c(3, 3)
      if (dist_names[i] == "lkj_correlation_distribution") {
        param_list$eta <- 1
        new_dim <- 3
      }
      if (dist_names[i] == "multivariate_normal_distribution") {
        param_list$mean <- matrix(rep(0, 5), nrow = 1)
        param_list$Sigma <- diag(5)
        new_dim <- 3
      }
      if (dist_names[i] == "dirichlet_distribution") {
        param_list$alpha <- matrix(rep(1, 5), nrow = 1)
        new_dim <- 3
      }
      
      init <- do.call(dist_tmp, param_list) 
      
      expanded <- change_dims(init, new_dim = new_dim)
      
    }
    
  }
  
  
}) 

test_that("initial nodes removed correctly", {
  
}) 

test_that("changing dimensions of operation nodes works", {
  
  for (i in seq_along(valid_ops)) {
    
    op <- valid_ops[i]
    if (op %in% c("add", "subtract", "multiply", "divide", "power")) {
      op <- switch(op,
                   "add" = "+",
                   "subtract" = "-",
                   "multiply" = "*",
                   "divide" = "/",
                   "power" = "^")
      dist_tmp <- list(normal(0, 1), normal(0, 1))
    } else {
      dist_tmp <- list(normal(0, 1))
    }
    
    init <- do.call(op, dist_tmp)
    expanded <- change_dims(init, new_dim = new_dim)
    
  }
  
}) 

test_that("operation nodes reapplied correctly", {
  
}) 
