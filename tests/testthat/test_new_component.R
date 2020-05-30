library(odpc)
context("Test that function to get component from new data and fitted a works")

set.seed(1234)
N <- 200 #length of series
m <- 10 #number of series
set.seed(1234)
f <- rnorm(N + 1)
x_small <- matrix(0, N, m)
u <- matrix(rnorm(N * m), N, m)
for (i in 1:m) {
  x_small[, i] <- 10 * sin(2 * pi * (i/m)) * f[1:N] + 10 * cos(2 * pi * (i/m)) * f[2:(N + 1)] + u[, i]
}

odpc_fit_1 <- odpc(x_small, k=1)
new_f_1 <- get_new_comp(odpc_fit_1[[1]]$a, rolled_data=x_small, k_max=1)

odpc_fit_3 <- odpc(x_small, k=3)
new_f_3 <- get_new_comp(odpc_fit_3[[1]]$a, rolled_data=x_small, k_max=3)

test_that(paste0("New component matches old one for k=1"), {
  expect_equal(odpc_fit_1[[1]]$f, new_f_1, tolerance=1e-3)
})

test_that(paste0("New component matches old one for k=3"), {
  expect_equal(odpc_fit_3[[1]]$f, new_f_3, tolerance=1e-3)
})
