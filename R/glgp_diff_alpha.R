eig_sym_asc <- function(A) {
  e <- eigen(A, symmetric = TRUE)
  ord <- order(e$values)
  list(values = e$values[ord], vectors = e$vectors[, ord, drop = FALSE])
}

stabilized_inverse <- function(values, vectors, n) {
  kappa <- max(values) / min(values)
  delta <- max(0.0, min(values) * (kappa - exp(20.0)) / (exp(20.0) - 1.0))
  adj_values <- values + rep(delta, n)
  vectors %*% diag(1.0 / adj_values, nrow = n, ncol = n) %*% t(vectors)
}

stabilized_logdet <- function(values, n, delta_vec) {
  sum(log(values + delta_vec))
}

wg_kernel <- function(dist, theta, l_) {
  r <- sqrt(dist) / theta
  pmax(0.0, 1.0 - r)^(l_ + 1.0) * ((l_ + 1.0) * r + 1.0)
}

create_gp <- function(xy, x_test, gIndices, theta, predIndices, lNum, leaf_size, nugget = FALSE) {
  xy <- as.matrix(xy)
  x_test <- as.matrix(x_test)
  gIndices <- as.integer(gIndices)
  predIndices <- as.integer(predIndices)
  
  dim_ <- ncol(xy) - 1L
  gNum_ <- length(gIndices)
  l_ <- floor(dim_ / 2.0) + 2.0
  
  gParams_ <- rep(0.0, 2L * dim_ + 1L)
  yg_ <- xy[gIndices, dim_ + 1L]
  
  searchable_idx <- setdiff(seq_len(nrow(xy)), gIndices)
  searchable_coords <- xy[searchable_idx, seq_len(dim_), drop = FALSE]
  
  env <- new.env(parent = emptyenv())
  env$xy <- xy
  env$x_test <- x_test
  env$dim_ <- dim_
  env$gIndices <- gIndices
  env$gNum_ <- gNum_
  env$lNum_ <- as.integer(lNum)
  env$predIndices <- predIndices
  env$theta_ <- theta
  env$leaf_size_ <- leaf_size
  env$nugget_ <- nugget
  env$l_ <- l_
  env$gParams_ <- gParams_
  env$lam_ <- NA_real_
  env$nug_ <- NA_real_
  env$nug_local_ <- NA_real_
  env$yg_ <- yg_
  env$Rg_ <- matrix(0.0, gNum_, gNum_)
  env$Rl_ <- matrix(0.0, gNum_, gNum_)
  env$Ainv_ <- matrix(0.0, gNum_, gNum_)
  env$oneVecG_ <- rep(1.0, gNum_)
  env$oneVecL_ <- rep(1.0, lNum)
  env$oneVecGL_ <- rep(1.0, gNum_ + lNum)
  env$searchable_idx <- searchable_idx
  env$searchable_coords <- searchable_coords
  
  env
}

find_neighbors <- function(gp, query_row, k) {
  if (length(gp$searchable_idx) == 0L) {
    stop("No searchable points remain after removing global indices.")
  }
  query_mat <- matrix(query_row[seq_len(gp$dim_)], nrow = 1L)
  k_eff <- min(k, nrow(gp$searchable_coords))
  nn <- FNN::get.knnx(gp$searchable_coords, query_mat, k = k_eff)
  gp$searchable_idx[nn$nn.index[1L, ]]
}

find_RgRl <- function(gp) {
  dim_ <- gp$dim_
  gNum_ <- gp$gNum_
  gIndices <- gp$gIndices
  gParams_ <- gp$gParams_
  theta_ <- gp$theta_
  l_ <- gp$l_
  
  Rg <- matrix(0.0, gNum_, gNum_)
  Rl <- matrix(0.0, gNum_, gNum_)
  
  for (i in seq_len(gNum_)) {
    for (j in i:gNum_) {
      if (i == j) {
        Rg[i, j] <- 1.0
        Rl[i, j] <- 1.0
      } else {
        row_i <- gp$xy[gIndices[i], seq_len(dim_)]
        row_j <- gp$xy[gIndices[j], seq_len(dim_)]
        value <- 0.0
        dist <- 0.0
        for (d in seq_len(dim_)) {
          temp <- abs(row_i[d] - row_j[d])
          value <- value - gParams_[d] * temp^gParams_[dim_ + 1L]
          dist <- dist + temp^2
        }
        w <- wg_kernel(dist, theta_, l_)
        Rl[i, j] <- w
        Rl[j, i] <- w
        Rg[i, j] <- exp(value)
        Rg[j, i] <- Rg[i, j]
      }
    }
  }
  
  gp$Rg_ <- Rg
  gp$Rl_ <- Rl
  invisible(gp)
}

find_Ainv <- function(gp, lam, nugget) {
  gNum_ <- gp$gNum_
  A <- (1.0 - lam) * gp$Rg_ + lam * gp$Rl_
  diag(A) <- diag(A) + nugget
  
  eig <- eig_sym_asc(A)
  gp$Ainv_ <- stabilized_inverse(eig$values, eig$vectors, gNum_)
  invisible(gp)
}

get_nllg <- function(gp, gParams) {
  dim_ <- gp$dim_
  gNum_ <- gp$gNum_
  gIndices <- gp$gIndices
  yg_ <- gp$yg_
  oneVecG_ <- gp$oneVecG_
  
  nugget <- exp(gParams[2L * dim_ + 1L])
  Rg <- matrix(0.0, gNum_, gNum_)
  
  for (i in seq_len(gNum_)) {
    for (j in i:gNum_) {
      if (i == j) {
        Rg[i, j] <- 1.0 + nugget
      } else {
        row_i <- gp$xy[gIndices[i], seq_len(dim_)]
        row_j <- gp$xy[gIndices[j], seq_len(dim_)]
        value <- 0.0
        for (d in seq_len(dim_)) {
          temp <- abs(row_i[d] - row_j[d])
          value <- value - gParams[d] * temp^gParams[dim_ + 1L]
        }
        Rg[i, j] <- exp(value)
        Rg[j, i] <- Rg[i, j]
      }
    }
  }
  
  eig <- eig_sym_asc(Rg)
  kappa <- max(eig$values) / min(eig$values)
  delta <- max(0.0, min(eig$values) * (kappa - exp(20.0)) / (exp(20.0) - 1.0))
  delta_vec <- rep(delta, gNum_)
  
  RgInv <- eig$vectors %*% diag(1.0 / (eig$values + delta_vec), gNum_, gNum_) %*% t(eig$vectors)
  RgLogDet <- stabilized_logdet(eig$values, gNum_, delta_vec)
  
  mu <- (matrix(colSums(RgInv), nrow = 1L) %*% yg_) / sum(RgInv)
  yg_mu <- yg_ - oneVecG_ * mu[1, 1]
  tau2 <- (1.0 / gNum_) * t(yg_mu) %*% RgInv %*% yg_mu
  
  gNum_ * log(tau2[1, 1]) + RgLogDet
}

predict_point <- function(gp, ind, lam, nugget, test = FALSE, return_sigma = FALSE) {
  dim_ <- gp$dim_
  gNum_ <- gp$gNum_
  lNum_ <- gp$lNum_
  gIndices <- gp$gIndices
  gParams_ <- gp$gParams_
  theta_ <- gp$theta_
  l_ <- gp$l_
  yg_ <- gp$yg_
  Ainv_ <- gp$Ainv_
  oneVecG_ <- gp$oneVecG_
  oneVecL_ <- gp$oneVecL_
  oneVecGL_ <- gp$oneVecGL_
  
  nn <- if (test) lNum_ else lNum_ + 1L
  query_row <- if (test) {
    gp$x_test[ind, ]
  } else {
    gp$xy[ind, ]
  }
  
  index <- find_neighbors(gp, query_row, nn)
  local_start <- nn - lNum_ + 1L
  local_indices <- index[local_start:nn]
  
  yl <- gp$xy[local_indices, dim_ + 1L]
  y <- c(yg_, yl)
  
  D <- matrix(0.0, lNum_, lNum_)
  for (u in seq_len(lNum_)) {
    for (v in u:lNum_) {
      if (u == v) {
        D[u, v] <- 1.0 + nugget
      } else {
        row_u <- gp$xy[local_indices[u], seq_len(dim_)]
        row_v <- gp$xy[local_indices[v], seq_len(dim_)]
        value <- 0.0
        dist <- 0.0
        for (d in seq_len(dim_)) {
          temp <- abs(row_u[d] - row_v[d])
          value <- value - gParams_[d] * temp^gParams_[dim_ + 1L]
          dist <- dist + temp^2
        }
        w <- wg_kernel(dist, theta_, l_)
        D[u, v] <- (1.0 - lam) * exp(value) + lam * w
        D[v, u] <- D[u, v]
      }
    }
  }
  
  B <- matrix(0.0, gNum_, lNum_)
  for (u in seq_len(gNum_)) {
    for (v in seq_len(lNum_)) {
      row_u <- gp$xy[gIndices[u], seq_len(dim_)]
      row_v <- gp$xy[local_indices[v], seq_len(dim_)]
      value <- 0.0
      dist <- 0.0
      for (d in seq_len(dim_)) {
        temp <- abs(row_u[d] - row_v[d])
        value <- value - gParams_[d] * temp^gParams_[dim_ + 1L]
        dist <- dist + temp^2
      }
      w <- wg_kernel(dist, theta_, l_)
      B[u, v] <- (1.0 - lam) * exp(value) + lam * w
    }
  }
  
  CAinv <- t(B) %*% Ainv_
  S <- D - CAinv %*% B
  
  eigS <- eig_sym_asc(S)
  kappa <- max(eigS$values) / min(eigS$values)
  delta <- max(0.0, min(eigS$values) * (kappa - exp(20.0)) / (exp(20.0) - 1.0))
  delta_vec <- oneVecL_ * delta
  Sinv <- eigS$vectors %*% diag(1.0 / (eigS$values + delta_vec), lNum_, lNum_) %*% t(eigS$vectors)
  
  CAinvTSinv <- t(CAinv) %*% Sinv
  Rinv <- matrix(0.0, gNum_ + lNum_, gNum_ + lNum_)
  Rinv[seq_len(gNum_), seq_len(gNum_)] <- Ainv_ + CAinvTSinv %*% CAinv
  Rinv[seq_len(gNum_), gNum_ + seq_len(lNum_)] <- -CAinvTSinv
  Rinv[gNum_ + seq_len(lNum_), seq_len(gNum_)] <- -t(CAinvTSinv)
  Rinv[gNum_ + seq_len(lNum_), gNum_ + seq_len(lNum_)] <- Sinv
  
  row <- query_row[seq_len(dim_)]
  rVec <- numeric(gNum_ + lNum_)
  for (u in seq_len(gNum_ + lNum_)) {
    row_u <- if (u <= gNum_) {
      gp$xy[gIndices[u], seq_len(dim_)]
    } else {
      gp$xy[local_indices[u - gNum_], seq_len(dim_)]
    }
    value <- 0.0
    dist <- 0.0
    for (d in seq_len(dim_)) {
      temp <- abs(row_u[d] - row[d])
      value <- value - gParams_[d] * temp^gParams_[dim_ + 1L]
      dist <- dist + temp^2
    }
    w <- wg_kernel(dist, theta_, l_)
    rVec[u] <- (1.0 - lam) * exp(value) + lam * w
  }
  
  mu <- (matrix(colSums(Rinv), nrow = 1L) %*% y) / sum(Rinv)
  y_mu <- y - oneVecGL_ * mu[1, 1]
  prediction <- mu + t(rVec) %*% Rinv %*% y_mu
  pred <- prediction[1, 1]
  
  if (!return_sigma) {
    return(pred)
  }
  
  tau2 <- (1.0 / (gNum_ + lNum_)) * t(y_mu) %*% Rinv %*% y_mu
  quad <- (t(rVec) %*% Rinv %*% rVec)[1, 1]
  sigma <- sqrt(tau2[1, 1] * max(1e-7, 1.0 + nugget - quad))
  list(mu = pred, sigma = sigma)
}

predict_point_global <- function(gp, ind, nugget, test = FALSE) {
  predict_point(gp, ind, lam = 0, nugget = nugget, test = test, return_sigma = FALSE)
}

get_mse <- function(gp, lam, nugget) {
  dim_   <- gp$dim_
  nugget <- (1.0 - lam) * gp$gParams_[2L * dim_ + 1L] + lam * nugget
  find_Ainv(gp, lam, nugget)
  
  mse <- 0.0
  for (i in seq_along(gp$predIndices)) {
    ind  <- gp$predIndices[i]
    pred <- predict_point(gp, ind, lam, nugget, test = FALSE, return_sigma = FALSE)
    mse  <- mse + (gp$xy[ind, dim_ + 1L] - pred)^2
  }
  mse
}

#' Recycle a scalar or per-dimension bound vector to length \code{dim_}.
normalize_hyperparam_bound <- function(bound, dim_, name) {
  bound <- as.numeric(bound)
  if (length(bound) == 1L) {
    return(rep(bound, dim_))
  }
  if (length(bound) != dim_) {
    stop(sprintf(
      "%s must have length 1 or %d (got %d).",
      name, dim_, length(bound)
    ), call. = FALSE)
  }
  bound
}

#' Convert physical lengthscale bounds to GLGP global-kernel coefficient bounds.
#'
#' GLGP uses \code{exp(-sum(a * abs(h)^alpha))}; with fixed \code{alpha}, this
#' matches a Gaussian kernel when \code{a = 1 / (2 * lengthscale^alpha)}.
a_bounds_from_lengthscales <- function(lengthscale_lower, lengthscale_upper, alpha) {
  lengthscale_lower <- as.numeric(lengthscale_lower)
  lengthscale_upper <- as.numeric(lengthscale_upper)
  alpha <- as.numeric(alpha)
  if (any(lengthscale_lower <= 0) || any(lengthscale_upper <= 0)) {
    stop("lengthscale bounds must be positive.", call. = FALSE)
  }
  if (any(lengthscale_lower > lengthscale_upper)) {
    stop("lengthscale_lower cannot exceed lengthscale_upper.", call. = FALSE)
  }
  list(
    theta_lower = 1.0 / (2.0 * lengthscale_upper^alpha),
    theta_upper = 1.0 / (2.0 * lengthscale_lower^alpha)
  )
}

set_gParams_from_list <- function(gp, global_params) {
  dim_ <- gp$dim_
  required <- c("a", "alpha", "nugget")
  missing <- setdiff(required, names(global_params))
  if (length(missing) > 0L) {
    stop(
      "global_params must be a list with components: ",
      paste(required, collapse = ", "),
      ". Missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  
  a <- as.numeric(global_params$a)
  if (length(a) == 1L) {
    a <- rep(a, dim_)
  }
  if (length(a) != dim_) {
    stop(sprintf("global_params$a must have length 1 or %d.", dim_), call. = FALSE)
  }
  
  gParams <- numeric(dim_ + 2L)
  gParams[seq_len(dim_)] <- a
  gParams[dim_ + 1L] <- as.numeric(global_params$alpha)
  gParams[dim_ + 2L] <- as.numeric(global_params$nugget)
  gp$gParams_ <- gParams
  find_RgRl(gp)
  invisible(gp)
}

#' Estimate global-kernel hyperparameters by maximum likelihood.
#'
#' @param gp GP environment from \code{create_gp()}.
#' @param theta_lower Lower bound(s) on global length-scale coefficients
#'   \code{gParams[1:d]} (scalar or length-\code{d} vector).
#' @param theta_upper Upper bound(s) on global length-scale coefficients
#'   (scalar or length-\code{d} vector).
#' @param alpha_lower Lower bound on the power-exponential exponent \code{alpha}.
#' @param alpha_upper Upper bound on the power-exponential exponent \code{alpha}.
estimate_gParams <- function(
    gp,
    theta_lower = 1e-7,
    theta_upper = 1000,
    alpha_lower = 1.0,
    alpha_upper = 2.0
) {
  dim_ <- gp$dim_
  nugget_ <- gp$nugget_
  
  lb_theta <- normalize_hyperparam_bound(theta_lower, dim_, "theta_lower")
  ub_theta <- normalize_hyperparam_bound(theta_upper, dim_, "theta_upper")
  if (any(lb_theta >= ub_theta)) {
    stop("Each theta_lower value must be strictly less than theta_upper.", call. = FALSE)
  }
  if (alpha_lower >= alpha_upper) {
    stop("alpha_lower must be strictly less than alpha_upper.", call. = FALSE)
  }
  
  lb <- c(lb_theta, alpha_lower, log(1e-7))
  ub <- c(ub_theta, alpha_upper, 0.0)
  
  opt_dim <- if (nugget_) 2L * dim_ + 1L else 2L * dim_
  nugget_init <- if (nugget_) log(1e-3) else log(1e-7)
  
  rho <- 1.0 / sqrt(rep(dim_, dim_))
  if (alpha_lower == alpha_upper) {
    alpha_grid <- alpha_lower
  } else {
    alpha_grid <- exp(seq(log(alpha_upper), log(alpha_lower), length.out = 11L))
  }
  
  nllg_values <- vapply(alpha_grid, function(a) {
    gParams <- numeric(2L * dim_ + 1L)
    gParams[seq_len(dim_)] <- pmin(pmax(a * rho, lb_theta), ub_theta)
    gParams[dim_ + seq_len(dim_)] <- pmin(pmax(1.95, alpha_lower), alpha_upper)
    gParams[2L * dim_ + 1L] <- nugget_init
    get_nllg(gp, gParams)
  }, numeric(1))
  
  min_index <- which.min(nllg_values)
  num_threads <- max(1L, parallel::detectCores(logical = TRUE))
  num_opt <- min(11L, max(3L, num_threads))
  factor <- exp(seq(-0.5, 0.5, length.out = num_opt))
  max_eval <- as.integer(min(500.0, 100.0 * log(1.0 + dim_)))
  
  opt_results <- vector("list", num_opt)
  for (i in seq_len(num_opt)) {
    gParams <- numeric(2L * dim_ + 1L)
    gParams[seq_len(dim_)] <- pmin(
      pmax(alpha_grid[min_index] * rho * factor[i], lb_theta),
      ub_theta
    )
    gParams[dim_ + seq_len(dim_)] <- pmin(pmax(1.95, alpha_lower), alpha_upper)
    gParams[2L * dim_ + 1L] <- nugget_init
    
    res <- nloptr::nloptr(
      x0 = gParams[seq_len(opt_dim)],
      eval_f = function(x) {
        params <- gParams
        params[seq_len(opt_dim)] <- x
        get_nllg(gp, params)
      },
      lb = lb[seq_len(opt_dim)],
      ub = ub[seq_len(opt_dim)],
      opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = max_eval)
    )
    
    gParams[seq_len(opt_dim)] <- res$solution
    opt_results[[i]] <- list(nllg = res$objective, gParams = gParams)
  }
  
  best <- opt_results[[which.min(vapply(opt_results, `[[`, numeric(1), "nllg"))]]
  gParams <- best$gParams
  gParams[2L * dim_ + 1L] <- exp(gParams[2L * dim_ + 1L])
  gp$gParams_ <- gParams
  find_RgRl(gp)
  invisible(gp)
}

estimate_sParams <- function(gp) {
  dim_ <- gp$dim_
  nugget_ <- gp$nugget_
  
  lb <- c(log(1e-7), log(1e-7))
  ub <- c(log(0.999), 0.0)
  opt_dim <- if (nugget_) 2L else 1L
  max_eval <- 20L
  
  sParams <- c(log(1e-1), if (nugget_) log(1e-3) else log(1e-7))
  
  res <- nloptr::nloptr(
    x0 = sParams[seq_len(opt_dim)],
    eval_f = function(x) {
      lam <- exp(x[1L])
      nug <- if (length(x) >= 2L) exp(x[2L]) else exp(sParams[2L])
      get_mse(gp, lam, nug)
    },
    lb = lb[seq_len(opt_dim)],
    ub = ub[seq_len(opt_dim)],
    opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = max_eval)
  )
  
  sParams[seq_len(opt_dim)] <- res$solution
  gp$lam_ <- exp(sParams[1L])
  gp$nug_local_ <- exp(sParams[2L])
  gp$nug_ <- (1.0 - gp$lam_) * gp$gParams_[2L * dim_ + 1L] + gp$lam_ * gp$nug_local_
  invisible(gp)
}

gp_predict <- function(gp) {
  find_Ainv(gp, gp$lam_, gp$nug_)
  test_num <- nrow(gp$x_test)
  
  predictions <- numeric(test_num)
  sigmas <- numeric(test_num)
  
  for (i in seq_len(test_num)) {
    out <- predict_point(gp, i, gp$lam_, gp$nug_, test = TRUE, return_sigma = TRUE)
    predictions[i] <- out$mu
    sigmas[i] <- out$sigma
  }
  
  list(mu = predictions, sigma = sigmas)
}

#' Global-only predictions: \code{lam = 0}, global-kernel nugget.
gp_predict_global <- function(gp) {
  dim_ <- gp$dim_
  nugget <- gp$gParams_[2L * dim_ + 1L]
  find_Ainv(gp, lam = 0, nugget = nugget)
  test_num <- nrow(gp$x_test)
  
  predictions <- numeric(test_num)
  sigmas <- numeric(test_num)
  
  for (i in seq_len(test_num)) {
    out <- predict_point(gp, i, lam = 0, nugget = nugget, test = TRUE, return_sigma = TRUE)
    predictions[i] <- out$mu
    sigmas[i] <- out$sigma
  }
  
  list(mu = predictions, sigma = sigmas)
}

#' Local-only predictions: \code{lam = 1}, local-kernel nugget.
gp_predict_local <- function(gp) {
  nugget <- if (!is.null(gp$nug_local_) && is.finite(gp$nug_local_)) {
    gp$nug_local_
  } else {
    gp$nug_
  }
  find_Ainv(gp, lam = 1, nugget = nugget)
  test_num <- nrow(gp$x_test)
  
  predictions <- numeric(test_num)
  sigmas <- numeric(test_num)
  
  for (i in seq_len(test_num)) {
    out <- predict_point(gp, i, lam = 1, nugget = nugget, test = TRUE, return_sigma = TRUE)
    predictions[i] <- out$mu
    sigmas[i] <- out$sigma
  }
  
  list(mu = predictions, sigma = sigmas)
}

#' RMSE for hybrid, global-only, and/or local-only predictions.
#'
#' @param y Numeric vector of observed responses (same order as predictions).
#' @param result List returned by \code{glgp()} (must contain \code{mu}; optionally
#'   \code{global_mu} and \code{local_mu}).
#' @param y_mean, y_sd Optional mean / SD used to back-transform scaled predictions
#'   (\code{pred * y_sd + y_mean}). Defaults leave predictions on the same scale as
#'   \code{y}.
#' @return Named list with \code{rmse} (hybrid). When available, also
#'   \code{rmse_global} and \code{rmse_local}.
glgp_rmse <- function(y, result, y_mean = 0, y_sd = 1) {
  y <- as.numeric(y)
  backtransform <- function(mu) as.numeric(mu) * y_sd + y_mean
  out <- list(rmse = sqrt(mean((y - backtransform(result$mu))^2)))
  if (!is.null(result$global_mu)) {
    out$rmse_global <- sqrt(mean((y - backtransform(result$global_mu))^2))
  }
  if (!is.null(result$local_mu)) {
    out$rmse_local <- sqrt(mean((y - backtransform(result$local_mu))^2))
  }
  out
}

#' @param xy Training matrix (n x (d+1)), last column is response.
#' @param x_test Test inputs (m x (d+1)); response column ignored if present.
#' @param gIndices Integer indices of global design points (1-based).
#' @param theta Bandwidth parameter for Wendland kernel.
#' @param predIndices Indices used for MSE-based lambda/nugget tuning.
#' @param lNum Number of local neighbors.
#' @param nugget If TRUE, optimize nugget jointly in hyperparameter search.
#' @param leaf_size KD-tree leaf size (retained for API compatibility; FNN is used in R).
#' @param theta_lower Lower bound(s) on global length-scale coefficients passed to
#'   \code{estimate_gParams()} (scalar or length-\code{d} vector).
#' @param theta_upper Upper bound(s) on global length-scale coefficients (scalar or
#'   length-\code{d} vector).
#' @param alpha_lower Lower bound on the global power-exponential exponent.
#' @param alpha_upper Upper bound on the global power-exponential exponent.
#' @param global_params Optional list with components \code{a}, \code{alpha}, and
#'   \code{nugget}. When provided, global hyperparameters are fixed and
#'   \code{estimate_gParams()} is skipped.
#' @param predict_global If \code{TRUE}, also compute global-only (\code{lam = 0})
#'   and local-only (\code{lam = 1}) predictions so RMSEs can be separated into
#'   global vs local parts. Adds two extra prediction passes over test points.
#' @return List with \code{mu} and \code{sigma} (hybrid GLGP). When
#'   \code{predict_global = TRUE}, also includes \code{global_mu},
#'   \code{global_sigma}, \code{local_mu}, and \code{local_sigma}. Use
#'   \code{glgp_rmse(y_test, result)} to obtain hybrid / global / local RMSE.
glgp <- function(
    xy,
    x_test,
    gIndices,
    theta,
    predIndices,
    lNum,
    nugget = FALSE,
    leaf_size = 10L,
    theta_lower = 1e-7,
    theta_upper = 1000,
    alpha_lower = 1.0,
    alpha_upper = 2.0,
    global_params = NULL,
    predict_global = FALSE
) {
  gp <- create_gp(xy, x_test, gIndices, theta, predIndices, lNum, leaf_size, nugget)
  
  t0 <- proc.time()
  if (is.null(global_params)) {
    estimate_gParams(
      gp,
      theta_lower = theta_lower,
      theta_upper = theta_upper,
      alpha_lower = alpha_lower,
      alpha_upper = alpha_upper
    )
  } else {
    set_gParams_from_list(gp, global_params)
  }
  t_gParams <- (proc.time() - t0)["elapsed"]
  cat(sprintf("estimate_gParams : %.3f sec\n", t_gParams))
  
  t0 <- proc.time()
  estimate_sParams(gp)
  t_sParams <- (proc.time() - t0)["elapsed"]
  cat(sprintf("estimate_sParams : %.3f sec\n", t_sParams))
  
  t0 <- proc.time()
  result <- gp_predict(gp)
  if (predict_global) {
    global_result <- gp_predict_global(gp)
    local_result <- gp_predict_local(gp)
    result <- c(
      result,
      list(
        global_mu = global_result$mu,
        global_sigma = global_result$sigma,
        local_mu = local_result$mu,
        local_sigma = local_result$sigma
      )
    )
  }
  t_predict <- (proc.time() - t0)["elapsed"]
  cat(sprintf("gp_predict       : %.3f sec\n", t_predict))
  
  cat(sprintf("total            : %.3f sec\n", t_gParams + t_sParams + t_predict))
  
  result
}

glgp_cpp <- glgp
