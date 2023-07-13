#######################################################################
############_______________________________________________############
############                                               ############
############         Polya Tree Process Procedures         ############
############_______________________________________________############
############                                               ############
#######################################################################

#Load dependencies
library(parallel)
library(dplyr)
library(Matrix)

#Identify the possible quantile partitions for theta's on a given layer 
##theta_i => [limits[1,i], limits[2,i]) U [limits[2,i], limits[3,i])
theta.part <- function(layer){
  limits <- matrix(NA, nrow = 2^(layer-1), ncol = 3)
  for(i in 1:(2^(layer-1))){
    limits[i,] <- c((2*i - 2)/2^(layer),
                    (2*i - 1)/2^(layer),
                    (2*i)/2^(layer))
  }
  return(limits)
}

#Specification of the process' architecture
PT.prior <- function(qdist, layers, vartype = 'continuous', v = 1){
  #Wrong vartype
  if(!(vartype %in% c('continuous', 'discrete', 'singular'))){
    stop("Please specify an adequate vartype (continuous, discrete or singular).")
  }
  
  #Layer of a given theta
  layer <- numeric()
  for(i in 1:layers){
    layer <- c(layer, rep(i, 2^(i-1)))
  }
  
  #Quantiles of the intervals
  quant.part <- qdist(do.call(rbind, sapply(1:layers,theta.part)))
  
  #Parameter values
  if(vartype == 'continuous') alpha <- v*layer^2
  else if(vartype == 'discrete') alpha <- v*1/2^layer
  else alpha <- rep(v*1, length(layer))
  
  #Partially specified Polya tree
  PT.model <- data.frame(layer, quant.part, alpha, alpha)
  names(PT.model) <- c("layer", "i1", "i2", "i3", "alpha0", "alpha1")
  
  return(PT.model)
}

#PT posterior
PT.posterior <- function(Y, #Data
                         qdist = NULL, #Quantile function
                         layers = ceiling(log(length(Y), base = 2)), #Number of layers
                         vartype = 'continuous', #Variable type
                         v = 1, #proximity to the centering distribution
                         PT.model = NULL, #Previously specified Polya tree
                         cores = 1 #Number of cores for parallel processing
){
  if(is.null(qdist) & is.null(PT.model)){
    stop("Please specify either the quantile function or the PT model")
  }
  
  #Generate the Polya tree if not previously specified
  if(is.null(PT.model)) PT.model <- PT.prior(qdist, layers, vartype, v)
  
  last_layer <- max(PT.model$layer)
  PT_ind <-which(PT.model$layer == last_layer) 
  all_ints <- c(as.numeric(t(as.matrix(PT.model[PT_ind, 2:3]))), as.numeric(tail(PT.model,1)[4]))
  
  mat.up <- matrix(0, nrow = 2^(last_layer - 1), ncol = 2)
  
  if(cores == 1) where_y <- sapply(Y, function(x) findInterval(x, all_ints) - 1)
  else{
    cl <- makeCluster(cores)
    where_y <- parSapply(cl, Y, function(x, all_ints) findInterval(x, all_ints) - 1, all_ints = all_ints)
    stopCluster(cl)
  }
  row_ind <- as.integer(where_y/2) + 1
  col_ref <- where_y%%2 + 1
  ind_n <- data.frame(row_ind, col_ref) %>%
    group_by(row_ind, col_ref) %>%
    summarise(num = n(), .groups = 'drop')
  
  for(i in 1:dim(ind_n)[1]){
    my_ind <- ind_n[i,]
    mat.up[my_ind$row_ind, my_ind$col_ref] <- my_ind$num
  }
  
  PT.model[PT_ind, c(5,6)] <- PT.model[PT_ind, c(5,6)] + mat.up
  
  if(last_layer > 1){
    for(layer in (last_layer - 1):1){
      mat.up <- matrix(rowSums(mat.up), ncol = 2, byrow = T)
      PT_ind <- which(PT.model$layer == layer) 
      PT.model[PT_ind, c(5,6)] <- PT.model[PT_ind, c(5,6)] + mat.up
    }
  }
  return(PT.model)
}

#Sampling distributions with unconditional log-prob
rPT <- function(n, PT.model, seed = NULL){
  par_num <- dim(PT.model)[1]
  layers <- unique(PT.model$layer)
  last_layer <- max(layers)
  sam_mat <- array(NA, c(par_num, 2, n))
  set.seed(42)
  for(k in 1:n){
    theta <- rbeta(par_num, PT.model$alpha0, PT.model$alpha1)
    log_prob <- log_theta <- log(as.vector(rbind(theta, 1-theta)))
    for(layer in 1:(last_layer-1)){
      ind <- (sum(2^(1:layer)) + 1):(2*par_num)
      theta_rep <- log_theta[1:(2*par_num - sum(2^(last_layer + 1 - 1:layer)))]
      log_prob[ind] <- log_prob[ind] + rep(theta_rep, each = 2^layer)
    }
    sam_mat[,,k] <- matrix(log_prob, ncol = 2, byrow = T)
  }
  return(sam_mat)
}

#Sampling observations
rPTdist <- function(n, PT.model, dist.samples, pdist, qdist, seed = NULL){
  last_layer <- max(PT.model$layer)
  ind <- which(PT.model$layer == last_layer)
  sam_obs <- matrix(NA, nrow = n, ncol = dim(dist.samples)[3])
  for(j in 1:dim(dist.samples)[3]){
    probs <- exp(as.vector(t(dist.samples[ind,,j])))
    set.seed(seed)
    sam_ind <- sample(0:(length(probs)-1), n, replace = T, prob = probs)
    round_sam <- as.integer(sam_ind/2)
    row_ind <- ind[1] + round_sam
    col_ref <- sam_ind%%2 ==0
    sam_int <- matrix(NA, nrow = length(row_ind), ncol = 2)
    sam_int[which(col_ref),] <- as.matrix(PT.model[row_ind[col_ref], 2:3])
    sam_int[which(!col_ref),] <- as.matrix(PT.model[row_ind[!col_ref], 3:4])
    sam_obs[,j] <- qdist(runif(n, pdist(sam_int[,1]), pdist(sam_int[,2])))
  }
  return(sam_obs)
}

#Density function
dPTdist <- function(X, PT.model, dist.sample, ddist, pdist, log = F){
  last_layer <- max(PT.model$layer)
  ind <- which(PT.model$layer == last_layer)
  all_ints <- c(as.numeric(t(as.matrix(PT.model[ind, 2:3]))), as.numeric(tail(PT.model,1)[4]))
  den_x <- numeric()
  len_x <- numeric()
  for(i in 1:length(X)){
    x <- X[i]
    where_x <- findInterval(x, all_ints) - 1
    row_ind <- ind[1] + as.integer(where_x/2)
    col_ref <- where_x%%2 + 1
    reg_x <- as.numeric(PT.model[row_ind, 1:2 + col_ref])
    den_x[i] <- dist.sample[row_ind, col_ref] +
      base::log(ddist(x)) - base::log(pdist(reg_x[2]) - pdist(reg_x[1]))
    len_x[i] <- length(den_x)
  }
  if(!log) den_x <- exp(den_x)
  return(den_x)
}

#Distribution function
pPTdist <- function(X, PT.model, dist.sample, pdist, log = F){
  last_layer <- max(PT.model$layer)
  ind <- which(PT.model$layer == last_layer)
  all_ints <- c(as.numeric(t(as.matrix(PT.model[ind, 2:3]))), as.numeric(tail(PT.model,1)[4]))
  prob_x <- numeric()
  for(i in 1:length(X)){
    x <- X[i]
    where_x <- findInterval(x, all_ints) - 1
    row_ind <- ind[1] + as.integer(where_x/2)
    col_ref <- where_x%%2 + 1
    reg_x <- as.numeric(PT.model[row_ind, 1:2 + col_ref])
    prob_x[i] <- (row_ind > ind[1])*sum(exp(dist.sample[ind[1]:(row_ind-1), ])) +
      (col_ref == 2)*exp(dist.sample[row_ind, 1]) + 
      exp(dist.sample[row_ind, col_ref])*(pdist(x) - pdist(reg_x[1]))/(pdist(reg_x[2]) - pdist(reg_x[1]))
  }
  if(log) prob_x <- log(prob_x)
  return(prob_x)
}

#Quantile function
qPTdist <- function(p, PT.model, dist.sample, pdist, qdist){
  last_layer <- max(PT.model$layer)
  ind <- which(PT.model$layer == last_layer)
  cum_probs <- cumsum(exp(t(dist.sample[ind,])))
  q_p <- numeric()
  for(i in 1:length(p)){
    prob <- p[i]
    ind_less <- tail(which(cum_probs < prob), 1)
    row_ind <- min(sum(c(ind[1], as.integer((ind_less + 1)/2))),
                   dim(PT.model)[1])
    col_ref <- sum(c(ind_less, 0))%%2 + 1
    reg_p <- as.numeric(PT.model[row_ind, 1:2 + col_ref])
    true_p <- sum(c((prob - cum_probs[ind_less + (col_ref == 2)])*(pdist(reg_p[2]) - pdist(reg_p[1])), pdist(reg_p[1])))
    q_p[i] <- qdist(true_p)
  }
  return(q_p)
}