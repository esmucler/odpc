is.odpc <- function(object, ...) {
  # This function checks whether an object is of class odpc
  # First check if object is a list.
  # Second check if the list has the correct entries
  # Third check if the entries have the correct classes
  # Fourth check if the entries have the correct dims
  if (any(!inherits(object, "list"), !inherits(object, "odpc"))) {
    return(FALSE)
  } else if (any(is.null(object$f), is.null(object$B), is.null(object$alpha),
                 is.null(object$a), is.null(object$mse),
                 is.null(object$k1), is.null(object$k2), is.null(object$call),
                 is.null(object$conv))) {
    return(FALSE)
  } else if (any(!inherits(object$mse, "numeric"), !inherits(object$a, "numeric"),
                 !inherits(object$B, "matrix"), !inherits(object$call, "call"), !inherits(object$conv, "logical"),
                 all(!inherits(object$f,"numeric"), !inherits(object$f, "ts"), !inherits(object$f, "xts")),
                 all(!inherits(object$k1, "numeric"), !inherits(object$k1, "integer")),
                 all(!inherits(object$k2, "numeric"), !inherits(object$k2, "integer")),
                 !inherits(object$alpha, "numeric"))) {
    return(FALSE)
  } else if (any(length(object$a) != (object$k1 + 1) * dim(object$B)[2],
                 dim(object$B)[1] != object$k2 + 1, length(object$alpha) != dim(object$B)[2])
             ) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

construct.odpc <- function(out, data) {
  #This function constructs an object of class odpc.
  #INPUT
  # out: the output of odpc_priv
  # data: the data matrix passed to odpc_priv
  #OUTPUT
  # An object of class odpc, that is, a list with entries:
  # a: vector to construct the principal component
  # alpha: vector of intercepts corresponding to the principal component
  # B: matrix of loadings corresponding to the principal component
  # k: number of lags used
  # f: principal component
  # mse:  mean (in T and m) squared error 
  # conv: logical. Did the iterations converge?
  # expart: proportion of the variance explained
  # call: the matched call
  # conv: logical. Did the iterations converge?
  
  out$f <- as.vector(out$f)
  out$B <- as.matrix(out$B)
  colnames(out$B) <- colnames(data)
  out$a <- as.vector(out$a)
  out$alpha <- as.vector(out$alpha)
  out$res <- NULL
  out$criter <- NULL
  class(out) <- append(class(out), "odpc")
  return(out)
}

fitted.odpc <- function(object, ...) {
  # Returns the fitted values of a odpc object
  matF <- getMatrixFitted(object$f, object$k2, object$k2)
  fitted <- matF %*% rbind(as.vector(object$alpha), as.matrix(object$B))
  colnames(fitted) <- colnames(object$B)
  return(fitted)
}

plot.odpc <- function(x, which = "Component", which_load = 0, ...) {
  #Plots a odpc object
  #INPUT
  # x: An object of class odpc, the result of odpc or one of the entries 
  # of the result of auto.odpc
  # which: String. Indicates what to plot, either 'Component' or 'Loadings'
  # Default is 'Component'
  # which_load: Lag number indicating which loadings should be plotted. 
  # Only used if which = 'Loadings'. Default is 0.
  if (!is.odpc(x)) {
    stop("x should be of class odpc")
  }
  if (!which %in% c("Component", "Loadings")) {
    stop("which should be either Component or Loadings ")
  }
  if (!inherits(which_load, "numeric")) {
    stop("which_load should be numeric")
  } else if (any(!(which_load == floor(which_load)), which_load < 0, which_load > nrow(x$B) - 1)) {
    stop("which_load should be a non-negative integer, at most equal to the number of lags")
  }
  if (which == "Component"){
    plot(x$f, type = "l", main = "Principal Component", ...) 
  } else if (which == "Loadings"){
    plot(x$B[which_load + 1, ], type = "l", main = c(paste(which_load, "lag loadings")), ...)
  }
}


print.odpc <- function(x, ...) {
  #Print method for the odpc class
  if (!is.odpc(x)) {
    stop("x should be of class odpc")
  }
  y <- list(x)
  class(y) <- c('list', 'odpcs')
  print(y)
}

is.odpcs <- function(object, ...) {
  #This function checks whether an object is of class odpcs,
  #that is, if each of its entries is a list of class odpc
  if (any(!inherits(object, "odpcs"), !inherits(object, "list"))) {
    return(FALSE)
  } else {
    return(all(sapply(object, is.odpc)))
  }
}

construct.odpcs <- function(out, data, fn_call) {
  #This function constructs an object of class odpcs that is, a list of length equal to 
  #the number of computed components. The i-th entry of this list is an object of class odpc.
  #INPUT
  # out: the output of auto.odpc
  # data: the data matrix passed to auto.odpc
  # fn_call: the original call to auto.odpc
  #OUTPUT
  # An object of class odpcs, that is, a list where each entry is an object of class odpc.
  
  out <- lapply(out, function(x, fn_call){ x$call <- fn_call; return(x)}, fn_call)
  out <- lapply(out, construct.odpc, data)
  class(out) <- append(class(out), "odpcs")
  return(out)
}

components <- function(object, ...){
  # Generic function for getting components out of an object
  UseMethod("components", object)
}


components.odpcs <- function(object, which_comp = 1) {
  # This function extracts the desired components from an odpcs object
  #INPUT
  # object: the output of odpc
  # which_comp: Integer vector. Indicates which components to get
  #OUTPUT
  # A list with the desired components as entries
  
  # if (!is.odpcs(object)) {
  #   stop("object should be of class odpcs")
  # }
  if (all(!inherits(which_comp, "numeric"), !inherits(which_comp, "integer"))) {
    stop("which_comp should be numeric")
  } else if (any(!(which_comp == floor(which_comp)), which_comp <= 0, which_comp > length(object))) {
    stop("The entries of which_comp should be positive integers, at most equal to the number of components")
  }
  object <- object[which_comp]
  comps <- lapply(object, function(object){ object$f })
  names(comps) <- paste("Component number", which_comp)
  return(comps)
}

fitted.odpcs <- function(object, num_comp = 1, ...) {
  # Returns the fitted values of a odpcs object using components 1,...,num_comp
  if (all(!inherits(num_comp, "numeric"), !inherits(num_comp, "integer"))) {
    stop("num_comp should be numeric")
  } else if (any(!(num_comp == floor(num_comp)), num_comp <= 0, num_comp > length(object))) {
    stop("num_comp should be a positive integer, at most equal to the number of components")
  }
  object <- object[1:num_comp]
  Nmin <- length(object[[num_comp]]$f) - object[[num_comp]]$k2
  fits <- 0
  aux <- 0
  for (i in 1:length(object)){
    aux <- fitted(object[[i]])
    N <- nrow(aux)
    fits <- fits + aux[(1 + N - Nmin):N,]
  }
  colnames(fits) <- colnames(object[[1]]$B)
  return(fits)
}


forecast.odpcs <- function(object, h, Z = NULL, add_residuals = FALSE, ...){
  # This function computes h-steps ahead forecast from a odpcs object
  # INPUT
  # object: an object of class odpcs
  # h: steps ahead to forecast
  # Z: original data.
  # add_residuals: logical? should the forecasts of the residuals be added?
  # ...: further arguments to be passed to auto.arima
  ncomp <- length(object)
  comps <- components.odpcs(object, 1:ncomp)
  fore <- 0
  
  fores_comps <- lapply(comps, function(x, h, ...) { auto <- auto.arima(x, ...)
                                                     return(forecast(auto, h)$mean)
                                                     }, h, ...)
  for (i in 1:ncomp){
    comps[[i]] <- c(comps[[i]], fores_comps[[i]])
    matF <- getMatrixFore(comps[[i]], object[[i]]$k2, h)
    fore <- fore + matF %*% rbind(as.vector(object[[i]]$alpha), as.matrix(object[[i]]$B))
  }
    
  
  if (add_residuals){
    k1s <- sapply(object, function(x) { return(x$k1) })
    k2s <- sapply(object, function(x) { return(x$k2) })
    k_tot <- sum(k1s + k2s)
    residuals <- Z[(k_tot + 1):nrow(Z), ] - fitted(object, num_comp = length(object))
    fore_res <- apply(residuals, 2, function(x, h, ...){ 
      auto <- auto.arima(x, ...)
      return(forecast(auto, h)$mean)
    }, h, ...)
    fore <- fore + fore_res
  }
  return(fore)
}

print.odpcs <- function(x, ...) {
  #   Print method for the odpcs class
  if (!is.odpcs(x)) {
    stop("x should be of class odpcs")
  }
  k1s <- sapply(x, function(x){ round(x$k1, 3) })
  k2s <- sapply(x, function(x){ round(x$k2, 3) })
  mses <- sapply(x, function(x){ round(x$mse, 3) })
  mat <- cbind(k1s, k2s, mses)
  nums <- paste(1:length(x))
  nums <- sapply(nums, function(x){ paste("Component", x)})
  colnames(mat) <- c("k1", "k2", "MSE")
  tab <- data.frame(mat, row.names = nums)
  print(tab)
}

new_window_object <- function(object_new, object_base){
  output <- list()
  num_k <- length(object_new)
  for (k in 1:num_k){
    output[[k]] <- merge_comps(object_new[[k]], object_base)
  }
  return(output)
}

merge_comps <- function(object_new, object_base){
  window_size <- length(object_new)
  output <- lapply(1:window_size, function(w) { c(object_base[[w]], object_new[[w]])  })
  return(output)
}

build_data_field <- function(object=NULL, Z=NULL, window_size=NULL, h=NULL){
  # Build a list of data-sets from the original data set (Z) or the residuals
  # in object
  if(is.null(object)){
    N <- nrow(Z)
    data_field <- lapply(1:window_size, function(w) { Z[1:(N - h - w + 1),] } )
  } else {
    num_comp <- length(object[[1]]) #we get the residuals from the last component
    data_field <- lapply(object, function(fit) { fit[[num_comp]]$res })
  }
  return(data_field)
}

get_best_fit <- function(object, Z, h, window_size, ...){
  # Find the k (indexing object) that gives the smallest h-forecast
  # mse evaluated using a rolling window of size window_size
  mses <- get_ave_mse(object, Z, h, window_size)
  ind_opt <- which.min(mses)
  opt_comp <- object[[ind_opt]] # this has window_size entries
  return(list(opt_comp=opt_comp, opt_mse = min(mses), opt_ind = ind_opt))
}
 
get_ave_mse <- function(object, Z, h, window_size, ...){
  # Get the mean mse of h-forecast mses from object: a list of lists, first level
  # is the candidate ks, second level the computation over the rolling window
  # of size window_size, each entry of the secnd level is a list (a la of odpcs)
  m <- ncol(Z)
  num_ks <- length(object)
  mses <- matrix(NA, window_size, num_ks) #Columns of the matrix correpond to different ks
  mses <- sapply(1:num_ks, function(ind){ get_vector_mses(object[[ind]], Z, h, window_size) })
  mses <- mses/m
  mses <- apply(mses, 2, mean)
  return(mses)
}

get_vector_mses <- function(object, Z, h, window_size){
  # Get a vector of h-forecast mses from object: a list one entry per computation over the rolling window
  # of size window_size.
  N <- nrow(Z)
  mses <- sapply(1:window_size, function(ind){ get_fore_mse(object[[ind]], Z[N - ind + 1,], h) })
  return(mses)
}

get_fore_mse <- function(comp, test, h, ...){
  fore <- forecast.odpcs(comp, h=h)
  fore <- fore[h,]
  mse <- sum((test - fore)^2)
  return(mse)
}

convert_rename <- function(fits){
  # The output of grid_odpc is complicated. It is a list, with one entry per element in k_list.
  # In the next level we get a matrix with one column, each row being a list corresponding to
  # a window. In there, there's another matrix with one column, each entry being the elements of a 
  # componente. This function takes this thing and outputs a list of lists. First level for elements
  # in k_list, next windows, in here are the components
  marrano <- lapply(fits, convert_rename_window)
  return(marrano)
}

convert_rename_window <- function(fits_by_w){
  fits_by_w <- lapply(fits_by_w[,1], convert_rename_comp)
  return(fits_by_w)
}

convert_rename_comp <- function(comp, wrap=TRUE){
  new_comp <- list(alpha=as.numeric(comp[,1][[1]]), B = as.matrix(comp[,1][[2]]), k2 = as.numeric(comp[,1][[3]]),
                   mse = as.numeric(comp[,1][[4]]), f = as.numeric(comp[,1][[5]]), res = as.matrix(comp[,1][[6]]),
                   k1=as.numeric(comp[,1][[7]]), criter=as.numeric(comp[,1][[8]]), conv = as.logical(comp[,1][[9]]),
                   a=as.numeric(comp[,1][[10]]))
  if(wrap){
    new_comp <- list(new_comp)
  }
  return(new_comp)
}

get_best_fit_crit <- function(object, Z, num_comp){
  crits <- sapply(object, function(comp){ get_crit(comp=comp, Z=Z, num_comp)})
  ind_opt <- which.min(crits)
  opt_comp <- object[[ind_opt]]
  opt_crit <- min(crits)
  return(list(opt_comp=opt_comp, opt_crit = opt_crit, opt_ind = ind_opt))
}

get_crit <- function(comp, Z, num_comp){
  m <- ncol(Z)
  T_c <- nrow(comp[[1]]$res)
  mse <- comp[[1]]$mse
  k <- comp[[1]]$k1
  crit <- T_c * log(mse) +  num_comp * (k + 1) * log(min(T_c, m)) 
}