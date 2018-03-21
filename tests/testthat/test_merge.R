library(odpc)
context("Test loop for building optimal odpcs over rolling window")

set.seed(1234)
N <- 50 #length of series
m <- 50 #number of series
f <- rnorm(N + 1)
x_small <- matrix(0, N, m)
u <- matrix(rnorm(N * m), N, m)
for (i in 1:m) {
  x_small[, i] <- 10 * sin(2 * pi * (i/m)) * f[1:N] + 10 * cos(2 * pi * (i/m)) * f[2:(N + 1)] + u[, i]
}

h <- 2
window_size <- 4
k_list <- c(1, 2, 3)
tol <- 1e-2
niter_max <- 20
ncores <- 1
method <- 2

data_field_orig <- build_data_field(Z=x_small, window_size = window_size, h = h)

fits <- grid_odpc(data_field = data_field_orig, k_list=k_list, window_size=window_size, tol=tol,
                  niter_max=niter_max, method=method)

best_fit <- get_best_fit(fits, Z=x_small, h=h, window_size = window_size)
opt_comp <- best_fit$opt_comp
new_best_mse <- best_fit$opt_mse
new_opt_k <- k_list[best_fit$opt_ind]
num_comp <- 0
ks <- c(new_opt_k)
num_comp <- length(ks)

while (num_comp <= 2){
  # data field is now formed by the residuals from the fit using the current optimal componentss
  data_field <- build_data_field(opt_comp) 
  
  # compute another component using the previous fitted ones
  fits <- grid_odpc(data_field = data_field, k_list=k_list, window_size=window_size, tol=tol,
                    niter_max=niter_max, method=method)
  # append to current components the new fitted ones
  extended_fits <- new_window_object(fits, opt_comp)
  
  best_fit <- get_best_fit(extended_fits, Z=x_small, h=h, window_size = window_size)
  
  opt_comp <- best_fit$opt_comp
  old_best_mse <- new_best_mse
  new_best_mse <- best_fit$opt_mse
  new_opt_k <- k_list[best_fit$opt_ind]
  ks <- c(ks, new_opt_k)
  num_comp <- length(ks)
}

for (t in 1:window_size){
    comp <- opt_comp[[t]] 
    comp <- odpc:::construct.odpcs(comp, data=x_small, fn_call=match.call(mean))
    N_t <- nrow(data_field_orig[[t]])
    min_recon <- 2 * sum(ks)+1
    dim(fitted(comp, num_comp = 3))
    mse_3comp <- mean((data_field_orig[[t]][min_recon:N_t,] - fitted(comp, num_comp = 3))**2)
    test_that(paste0('Correct mse for w=', t), {
    expect_lte(abs(1- mse_3comp/comp[[3]]$mse), 2e-2)
    })
}


