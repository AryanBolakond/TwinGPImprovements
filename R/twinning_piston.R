# =============================================================================
# Twinning Start Strategy Comparison for GLGP — Piston Function
# =============================================================================
# Compares three ways of running twinning to select GLOBAL points for GLGP:
#
#   (A) Default      — twin() with a random starting point (u1 = NULL)
#   (B) Extreme Point — twin() seeded from the point farthest from the mean,
#                       passed as u1
#   (C) Maximin Multiplet — multiplet() splits data into k groups, then one
#                           representative per group is chosen as the maximin
#                           point within that group (maximizes minimum distance
#                           between any two chosen global points)
#
# All three produce the same number of global points (n_global).
# Each set of gIndices is passed into the full GLGP pipeline and compared
# on RMSE, NLPD, and runtime.
# =============================================================================

library(FNN)
library(nloptr)
library(twinning)

# Source GLGP — adjust path if needed
source("C:/Users/Aryan/Downloads/RoshanJosephResearch/glgp.R")

# Piston Simulation function implementation and bounds
# (credit to sfu.ca/~ssurjano/piston.html)
piston <- function(xx) {
  M  <- xx[1]; S  <- xx[2]; V0 <- xx[3]; k  <- xx[4]
  P0 <- xx[5]; Ta <- xx[6]; T0 <- xx[7]
  
  A      <- P0*S + 19.62*M - k*V0/S
  V      <- (S/(2*k)) * (sqrt(A^2 + 4*k*(P0*V0/T0)*Ta) - A)
  fact2  <- k + (S^2)*(P0*V0/T0)*(Ta/V^2)
  
  2 * pi * sqrt(M / fact2)
}

piston_bounds <- list(
  M  = c(30,    60),
  S  = c(0.005, 0.020),
  V0 = c(0.002, 0.010),
  k  = c(1000,  5000),
  P0 = c(90000, 110000),
  Ta = c(290,   296),
  T0 = c(340,   360)
)

unit_to_physical <- function(u) {
  lo <- sapply(piston_bounds, `[`, 1)
  hi <- sapply(piston_bounds, `[`, 2)
  lo + u * (hi - lo)
}
piston_unit <- function(u) piston(unit_to_physical(u))


# Simulation settings
set.seed(42)

n_train  <- 300
n_test   <- 150
dim_     <- 7
n_global <- 75    # target global points; r = n_train / n_global for twin()
lNum     <- 25
theta    <- 0.3
nugget   <- TRUE
r        <- max(2L, round(n_train / n_global))   # reduction ratio for twin()


# Generate data
X_train <- matrix(runif(n_train * dim_), ncol = dim_)
y_train <- apply(X_train, 1, piston_unit)

X_test  <- matrix(runif(n_test  * dim_), ncol = dim_)
y_test  <- apply(X_test,  1, piston_unit)

y_mean     <- mean(y_train)
y_sd       <- sd(y_train)
y_train_sc <- (y_train - y_mean) / y_sd

xy     <- cbind(X_train, y_train_sc)
x_test <- X_test


# METHOD A: DEFAULT — random starting point (u1 = NULL)
t0 <- proc.time()
gIndices_A <- twin(X_train, r = r, u1 = NULL, format_data = TRUE)
t_twin_A   <- (proc.time() - t0)["elapsed"]
cat(sprintf("[A] Default (random)   : %d global pts, twinning took %.3fs\n",
            length(gIndices_A), t_twin_A))

# METHOD B: EXTREME POINT — seed twin() from the point farthest from the mean
col_means  <- colMeans(X_train)
dists_mean <- rowSums(sweep(X_train, 2, col_means)^2)   # sq. dist from mean
u1_extreme <- which.max(dists_mean)
cat(sprintf("[B] Extreme point index: %d  (sq-dist from mean = %.4f)\n",
            u1_extreme, dists_mean[u1_extreme]))

t0 <- proc.time()
gIndices_B <- twin(X_train, r = r, u1 = u1_extreme, format_data = TRUE)
t_twin_B   <- (proc.time() - t0)["elapsed"]
cat(sprintf("[B] Extreme point      : %d global pts, twinning took %.3fs\n",
            length(gIndices_B), t_twin_B))

# -----------------------------------------------------------------------------
# METHOD C: MAXIMIN MULTIPLET
#
# We use k = n_global so we get exactly one representative per group.
# strategy = 1 works for any k; strategy = 2 is better but needs k = 2^m.
# -----------------------------------------------------------------------------
k_mult <- n_global
# Use strategy 2 if k is a power of 2, otherwise fall back to strategy 1
is_pow2   <- function(x) x >= 1 && bitwAnd(x, x - 1L) == 0L
strategy  <- if (is_pow2(k_mult)) 2L else 1L
cat(sprintf("[C] multiplet strategy : %d  (k=%d, power-of-2: %s)\n",
            strategy, k_mult, is_pow2(k_mult)))

t0 <- proc.time()
group_labels <- multiplet(X_train, k = k_mult, strategy = strategy,
                          format_data = TRUE)
t_mult <- (proc.time() - t0)["elapsed"]

# For each group, find the point that maximises min-distance to all points
# outside the group
X_scaled <- scale(X_train)    # use the same standardised space multiplet() uses

gIndices_C <- integer(k_mult)
for (g in seq_len(k_mult)) {
  in_g    <- which(group_labels == g)
  out_g   <- which(group_labels != g)
  
  if (length(in_g) == 1L) {
    gIndices_C[g] <- in_g
    next
  }
  # For each candidate in group g, compute min distance to any outside point
  D        <- as.matrix(FNN::get.knnx(X_scaled[out_g, , drop = FALSE],
                                      X_scaled[in_g,  , drop = FALSE],
                                      k = 1L)$nn.dist)
  min_dists <- D[, 1L]
  gIndices_C[g] <- in_g[which.max(min_dists)]
}
gIndices_C <- unique(gIndices_C)
t_twin_C   <- (proc.time() - t0)["elapsed"]
cat(sprintf("[C] Maximin multiplet  : %d global pts, total took %.3fs\n",
            length(gIndices_C), t_twin_C))

# Helper: build predIndices that never overlap with gIndices
make_predIndices <- function(gIndices, n_pred = 30L) {
  pool <- setdiff(seq_len(min(120L, n_train)), gIndices)
  pool[seq_len(min(n_pred, length(pool)))]
}

# Fit GLGP and return metrics + runtime
run_glgp <- function(label, gIndices, t_selection) {
  predIdx <- make_predIndices(gIndices)
  
  cat(sprintf("\n--- Fitting GLGP [%s] (%d global pts) ---\n",
              label, length(gIndices)))
  
  t0     <- proc.time()
  result <- glgp(
    xy          = xy,
    x_test      = x_test,
    gIndices    = gIndices,
    theta       = theta,
    predIndices = predIdx,
    lNum        = lNum,
    nugget      = nugget,
    leaf_size   = 20L
  )
  t_fit <- (proc.time() - t0)["elapsed"]
  
  mu_pred    <- result$mu    * y_sd + y_mean
  sigma_pred <- result$sigma * y_sd
  
  errors    <- y_test - mu_pred
  rmse      <- sqrt(mean(errors^2))
  mae       <- mean(abs(errors))
  r2        <- 1 - sum(errors^2) / sum((y_test - mean(y_test))^2)
  nlpd      <- mean(0.5*(errors/sigma_pred)^2 + log(sigma_pred) + 0.5*log(2*pi))
  lower95   <- mu_pred - 1.96 * sigma_pred
  upper95   <- mu_pred + 1.96 * sigma_pred
  cov95     <- mean(y_test >= lower95 & y_test <= upper95)
  
  list(
    label      = label,
    n_global   = length(gIndices),
    t_select   = t_selection,
    t_fit      = t_fit,
    t_total    = t_selection + t_fit,
    rmse       = rmse,
    mae        = mae,
    r2         = r2,
    nlpd       = nlpd,
    cov95      = cov95,
    mu_pred    = mu_pred,
    sigma_pred = sigma_pred
  )
}


# Run all three methods
res_A <- run_glgp("Default (random)",    gIndices_A, t_twin_A)
res_B <- run_glgp("Extreme point",       gIndices_B, t_twin_B)
res_C <- run_glgp("Maximin multiplet",   gIndices_C, t_twin_C)


# Summary table
results <- list(res_A, res_B, res_C)

summary_df <- data.frame(
  Method       = sapply(results, `[[`, "label"),
  N_Global     = sapply(results, `[[`, "n_global"),
  RMSE         = round(sapply(results, `[[`, "rmse"),  4),
  MAE          = round(sapply(results, `[[`, "mae"),   4),
  R2           = round(sapply(results, `[[`, "r2"),    4),
  NLPD         = round(sapply(results, `[[`, "nlpd"),  4),
  Coverage_95  = round(sapply(results, `[[`, "cov95") * 100, 1),
  t_total_s    = round(sapply(results, `[[`, "t_total"),  2)
)

cat("\n================= TWINNING STRATEGY COMPARISON =================\n")
print(summary_df, row.names = FALSE)
cat("=================================================================\n\n")

# Relative to Method A (default) as baseline
cat("Relative RMSE vs. Default baseline:\n")
cat(sprintf("  Extreme point    : %+.1f%%\n",
            100 * (res_B$rmse - res_A$rmse) / res_A$rmse))
cat(sprintf("  Maximin multiplet: %+.1f%%\n\n",
            100 * (res_C$rmse - res_A$rmse) / res_A$rmse))

cat("Relative NLPD vs. Default baseline:\n")
cat(sprintf("  Extreme point    : %+.4f\n", res_B$nlpd - res_A$nlpd))
cat(sprintf("  Maximin multiplet: %+.4f\n\n", res_C$nlpd - res_A$nlpd))

# Diagnostic plots — predicted vs actual for all three methods
lims <- range(c(y_test, res_A$mu_pred, res_B$mu_pred, res_C$mu_pred))

par(mfrow = c(1, 3), mar = c(4, 4, 4, 1))
for (res in results) {
  plot(y_test, res$mu_pred,
       pch = 19, cex = 0.6, col = "#2166ac88",
       xlim = lims, ylim = lims,
       xlab = "Actual cycle time (s)",
       ylab = "Predicted cycle time (s)",
       main = sprintf("%s\nRMSE=%.4f  NLPD=%.3f", res$label, res$rmse, res$nlpd))
  abline(0, 1, col = "red", lwd = 1.5)
  legend("topleft", bty = "n", cex = 0.75,
         legend = sprintf("R\u00b2=%.3f\nCov95=%.0f%%", res$r2, res$cov95 * 100))
}
par(mfrow = c(1, 1))
