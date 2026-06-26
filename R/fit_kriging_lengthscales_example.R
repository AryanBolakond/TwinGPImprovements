library(rkriging)

source("R/glgp.R")

set.seed(123)

dim_ <- 2
n_train <- 100
n_test <- 10

X_train <- matrix(runif(n_train * dim_), ncol = dim_)
y_train <- sin(2.0 * pi * X_train[, 1]) + X_train[, 2]^2
X_test <- matrix(runif(n_test * dim_), ncol = dim_)
y_test <- sin(2.0 * pi * X_test[, 1]) + X_test[, 2]^2

xy <- cbind(X_train, y_train)
x_test <- X_test

gIndices <- seq_len(12L)
predIndices <- 13L:20L
lNum <- 6L
theta <- 0.3
fixed_alpha <- 2.0

fit <- rkriging::Fit.Kriging(
  X = X_train[gIndices, , drop = FALSE],
  y = y_train[gIndices],
  interpolation = TRUE,
  fit = TRUE,
  model = "OK",
  kernel.parameters = list(
    type = "Gaussian",
    lengthscale = rep(0.25, dim_),
    lengthscale.lower.bound = rep(0.05, dim_),
    lengthscale.upper.bound = rep(2.0, dim_)
  ),
  nlopt.parameters = list(
    algorithm = "NLOPT_LN_SBPLX",
    maxeval = 25L
  )
)

rkriging_params <- rkriging::Get.Kriging.Parameters(fit)
rkriging_lengthscales <- rkriging_params$lengthscale

# rkriging Gaussian: exp(-sum(h^2 / (2 * ell^2))).
# GLGP global:       exp(-sum(a * abs(h)^alpha)).
a <- 1.0 / (2.0 * rkriging_lengthscales^fixed_alpha)
global_params <- list(
  a = a,
  alpha = fixed_alpha,
  nugget = 1e-7
)

result <- glgp(
  xy = xy,
  x_test = x_test,
  gIndices = gIndices,
  theta = theta,
  predIndices = predIndices,
  lNum = lNum,
  nugget = FALSE,
  # global_params = global_params  # Uncomment this line to use parameters from rkriging
)

stopifnot(
  identical(names(result), c("mu", "sigma", "global_mu", "global_sigma")),
  length(result$mu) == n_test,
  length(result$global_mu) == n_test,
  all(is.finite(unlist(result)))
)

print(result)
sqrt(mean((y_test - result$global_mu)^2))
sqrt(mean((y_test - result$mu)^2))
