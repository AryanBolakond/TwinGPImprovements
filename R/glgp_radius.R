# =============================================================================
# GLGP with radius-based local neighborhoods
# =============================================================================
# Sources glgp.R and overrides neighbor selection so that lNum is the
# *minimum* number of local neighbors within a data-driven radius, rather
# than a fixed neighbor count.
# =============================================================================

.glgp_radius_dir <- local({
  ofile <- tryCatch(
    normalizePath(sys.frames()[[1]]$ofile, winslash = "/"),
    error = function(e) NA_character_
  )
  if (!is.na(ofile)) dirname(ofile) else getwd()
})

source(file.path(.glgp_radius_dir, "glgp.R"), local = FALSE)

create_gp <- function(xy, x_test, gIndices, theta, predIndices, lNum, leaf_size, nugget = FALSE) {
  xy <- as.matrix(xy)
  x_test <- as.matrix(x_test)
  gIndices <- as.integer(gIndices)
  predIndices <- as.integer(predIndices)

  dim_ <- ncol(xy) - 1L
  gNum_ <- length(gIndices)
  l_ <- floor(dim_ / 2.0) + 2.0

  gParams_ <- rep(0.0, dim_ + 2L)
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
  env$yg_ <- yg_
  env$Rg_ <- matrix(0.0, gNum_, gNum_)
  env$Rl_ <- matrix(0.0, gNum_, gNum_)
  env$Ainv_ <- matrix(0.0, gNum_, gNum_)
  env$oneVecG_ <- rep(1.0, gNum_)
  env$local_radius_ <- NA_real_
  env$k_upper_ <- NA_integer_
  env$searchable_idx <- searchable_idx
  env$searchable_coords <- searchable_coords

  env
}

get_knn_search <- function(gp, query_row, k) {
  if (length(gp$searchable_idx) == 0L) {
    stop("No searchable points remain after removing global indices.", call. = FALSE)
  }
  query_mat <- matrix(query_row[seq_len(gp$dim_)], nrow = 1L)
  k_eff <- min(as.integer(k), nrow(gp$searchable_coords))
  nn <- FNN::get.knnx(gp$searchable_coords, query_mat, k = k_eff)
  list(
    indices = gp$searchable_idx[nn$nn.index[1L, ]],
    distances = nn$nn.dist[1L, ]
  )
}

find_neighbors <- function(gp, query_row, k) {
  get_knn_search(gp, query_row, k)$indices
}

#' Minimum squared-distance radius so every prediction location has enough neighbors.
#'
#' \code{lNum_} is treated as the minimum number of local neighbors. The returned
#' radius is the largest k-th neighbor distance (squared Euclidean) over all
#' training and test query locations, using \code{lNum_ + 1} for training points
#' that appear in the searchable set (leave-one-out).
compute_local_radius <- function(gp) {
  lMin <- gp$lNum_
  k_upper <- nrow(gp$searchable_coords)
  if (k_upper < lMin) {
    stop(sprintf(
      "Need at least %d searchable local neighbors, but only %d are available.",
      lMin, k_upper
    ), call. = FALSE)
  }

  max_radius <- 0.0

  k_min_for <- function(ind, test) {
    if (test) {
      lMin
    } else if (ind %in% gp$searchable_idx) {
      min(lMin + 1L, k_upper)
    } else {
      lMin
    }
  }

  radius_for_query <- function(query_row, k_min) {
    get_knn_search(gp, query_row, k_min)$distances[k_min]
  }

  for (i in seq_len(nrow(gp$x_test))) {
    max_radius <- max(max_radius, radius_for_query(gp$x_test[i, ], k_min_for(NA, TRUE)))
  }

  for (i in seq_len(nrow(gp$xy))) {
    max_radius <- max(max_radius, radius_for_query(gp$xy[i, ], k_min_for(i, FALSE)))
  }

  gp$local_radius_ <- max_radius
  gp$k_upper_ <- as.integer(k_upper)
  invisible(gp)
}

ensure_local_radius <- function(gp) {
  if (is.na(gp$local_radius_) || is.na(gp$k_upper_)) {
    compute_local_radius(gp)
  }
  invisible(gp)
}

select_local_neighbors <- function(gp, query_row, test, query_ind) {
  lMin <- gp$lNum_
  knn <- get_knn_search(gp, query_row, gp$k_upper_)

  in_radius <- knn$distances <= gp$local_radius_
  local_indices <- knn$indices[in_radius]

  if (!test) {
    local_indices <- setdiff(local_indices, query_ind)
  }

  if (length(local_indices) < lMin) {
    k_need <- if (!test && query_ind %in% gp$searchable_idx) lMin + 1L else lMin
    knn_fb <- get_knn_search(gp, query_row, k_need)
    local_indices <- knn_fb$indices
    if (!test && query_ind %in% gp$searchable_idx) {
      local_indices <- setdiff(local_indices, query_ind)
    }
    if (length(local_indices) < lMin) {
      stop(sprintf(
        "Could not find %d local neighbors for query index %s.",
        lMin, query_ind
      ), call. = FALSE)
    }
  }

  local_indices
}

predict_point <- function(gp, ind, lam, nugget, test = FALSE, return_sigma = FALSE) {
  ensure_local_radius(gp)
  dim_ <- gp$dim_
  gNum_ <- gp$gNum_
  gIndices <- gp$gIndices
  gParams_ <- gp$gParams_
  theta_ <- gp$theta_
  l_ <- gp$l_
  yg_ <- gp$yg_
  Ainv_ <- gp$Ainv_
  oneVecG_ <- gp$oneVecG_

  query_row <- if (test) {
    gp$x_test[ind, ]
  } else {
    gp$xy[ind, ]
  }

  query_ind <- if (test) NA_integer_ else ind
  local_indices <- select_local_neighbors(gp, query_row, test, query_ind)
  lNum_ <- length(local_indices)
  oneVecL_ <- rep(1.0, lNum_)
  oneVecGL_ <- rep(1.0, gNum_ + lNum_)

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

#' @param lNum Minimum number of local neighbors; a data-driven radius is chosen
#'   so every prediction location has at least this many neighbors within range.
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
  compute_local_radius(gp)

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
    result <- c(
      result,
      list(
        global_mu = global_result$mu,
        global_sigma = global_result$sigma
      )
    )
  }
  t_predict <- (proc.time() - t0)["elapsed"]
  cat(sprintf("gp_predict       : %.3f sec\n", t_predict))

  cat(sprintf("total            : %.3f sec\n", t_gParams + t_sParams + t_predict))

  result
}

glgp_cpp <- glgp
