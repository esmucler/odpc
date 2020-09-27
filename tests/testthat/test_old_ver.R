library(odpc)
context("Compare odpc output with stored values")

T <- 30 #length of series
m <- 10 #number of series
set.seed(1234)
f <- rnorm(T + 1)
x <- matrix(0, T, m)
u <- matrix(rnorm(T * m), T, m)
for (i in 1:m) {
  x[, i] <- 10 * sin(2 * pi * (i/m)) * f[1:T] + 10 * cos(2 * pi * (i/m)) * f[2:(T + 1)] + u[, i]
}
# from odpc 1.0.0
old_mse <- 0.775981
old_B <- rbind(c(0.6047292, 0.2579100, -0.2383648, -0.6071289, -0.72405177, -0.5911932, -0.2396958,  0.2089657,  0.5788413,  0.7264755127),
               c(0.4129345, 0.7099759,  0.7018187,  0.4096304,  0.04051445, -0.4245995, -0.7128675, -0.7038352, -0.4241004, -0.0000175697))
old_a <- c(0.281353097, -0.168680377, -0.365183874, -0.294046511, -0.262466187, -0.244247851,  0.038762303, -0.387626573,  0.416784899,  0.155071336,
          -0.090215542, -0.035780053,  0.131756690, -0.185041501, -0.062290283,  0.005003871, -0.227599030,  0.006187320,  0.268459624, -0.075306478)

fit <- odpc(x, ks = c(1))

test_that('Equality with stored a', {
  expect_equal(fit[[1]]$a, old_a, tolerance=1e-2)
})
test_that('Equality with stored B', {
  expect_equal(fit[[1]]$B, old_B, tolerance=1e-2)
})
test_that('Equality with stored mse', {
  expect_equal(fit[[1]]$mse, old_mse, tolerance=1e-2)
})

fit_grad <- odpc(x, ks = c(1), method='gradient')
old_grad_mse <- 0.8019487

test_that('Equality with stored mse for method gradient', {
  expect_equal(fit_grad[[1]]$mse, old_grad_mse, tolerance=1e-2)
})

N <- 50
m <- 30
set.seed(1234)
Z <- matrix(0, N, m)
u <- matrix(rnorm(N * m), N, m)
b1 <- rnorm(m) * 0.2
a1 <- rep(0, m)
a1[1] <- 1/sqrt(2)
a1[3] <- 1/sqrt(2)

for (t in 2:N){
  Z[t, ] <- b1 * (sum(a1 * Z[t - 1, ]) ) + u[t, ]
}

old_mse_sparse_var <- 0.8271775
old_lambda_sparse_var <- 0.2514091
sparse_fit <- crit.sparse_odpc(Z=Z, k_list=1, ncores=1)
 
test_that('Equality with stored mse for VAR model for single component sparse odpc', {
  expect_equal(sparse_fit[[1]]$mse, old_mse_sparse_var, tolerance=1e-2)
})

test_that('Equality with stored lambda for VAR model for single component sparse odpc', {
  expect_equal(sparse_fit[[1]]$lambda, old_lambda_sparse_var, tolerance=1e-4)
})


old_mse_sparse_dfm <- 0.7941179
old_lambda_sparse_dfm <- 1.375486
sparse_fit <- crit.sparse_odpc(Z=x, k_list=1, ncores=1)

test_that('Equality with stored mse for DFM model for single component sparse odpc', {
  expect_equal(sparse_fit[[1]]$mse, old_mse_sparse_dfm, tolerance=1e-2)
})

test_that('Equality with stored lambda for DFM model for single component sparse odpc', {
  expect_equal(sparse_fit[[1]]$lambda, old_lambda_sparse_dfm, tolerance=1e-4)
})

T <- 30 #length of series
m <- 10 #number of series
set.seed(1234)
f <- rnorm(T + 1)
f2 <- rnorm(T + 1)
x <- matrix(0, T, m)
u <- matrix(rnorm(T * m), T, m)
for (i in 1:m) {
  x[, i] <- 10 * sin(2 * pi * (i/m)) * f[1:T] + 10 * cos(2 * pi * (i/m)) * f[2:(T + 1)] + 
    (i/m) * f2[1:T] + f2[2:(T + 1)]+ u[, i]
}

old_mse_1_sparse_dfm_two_fact <- 2.67426
old_lambda_1_sparse_dfm_two_fact <- 0.1478385
old_mse_2_sparse_dfm_two_fact <-0.6866781
old_lambda_2_sparse_dfm_two_fact <- 0.08996326
sparse_fit <- crit.sparse_odpc(Z=x, k_list=1, ncores=1, max_num_comp=2)

test_that('Equality with stored mse for two factor DFM model for first component sparse odpc', {
  expect_equal(sparse_fit[[1]]$mse, old_mse_1_sparse_dfm_two_fact, tolerance=1e-2)
})

test_that('Equality with stored lambda for two factor DFM model for first component sparse odpc', {
  expect_equal(sparse_fit[[1]]$lambda, old_lambda_1_sparse_dfm_two_fact, tolerance=1e-2)
})

test_that('Equality with stored mse for two factor DFM model for second component sparse odpc', {
  expect_equal(sparse_fit[[2]]$mse, old_mse_2_sparse_dfm_two_fact, tolerance=1e-2)
})

test_that('Equality with stored lambda for two factor DFM model for second component sparse odpc', {
  expect_equal(sparse_fit[[2]]$lambda, old_lambda_2_sparse_dfm_two_fact, tolerance=1e-3)
})