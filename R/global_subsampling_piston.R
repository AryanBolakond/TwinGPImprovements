# =============================================================================
# Compare Global Subsampling Methods using FIRST package 
# =============================================================================
# Three ways to choose global points for GLGP:
#
#   (A) Use all variables with no screening
#   (B) Select important variables using FIRST
#   (C) Create weights by rescaling each column after FIRST is used (use sqrt of importance)
#
# The variable set is then the input to twinning that produces gIndices.
#
# For each method, the resulting gIndices are passed into GLGP
# pipeline (create_gp / estimate_gParams / estimate_sParams / gp_predict) and
# run time, RMSE, NLPD, etc. are recorded and a plot is created.
# =============================================================================


library(FNN)
library(nloptr)
library(first)
library(twinning)

# Change source as necessary
source("C:/Users/Aryan/Downloads/RoshanJosephResearch/glgp.R")

# Piston Simulation function implementation and bounds
# (credit to sfu.ca/~ssurjano/piston.html)
piston <- function(xx) {
  M  <- xx[1]; S  <- xx[2]; V0 <- xx[3]; k  <- xx[4]
  P0 <- xx[5]; Ta <- xx[6]; T0 <- xx[7]
  
  Aterm1 <- P0 * S
  Aterm2 <- 19.62 * M
  Aterm3 <- -k * V0 / S
  A <- Aterm1 + Aterm2 + Aterm3
  
  Vfact1 <- S / (2 * k)
  Vfact2 <- sqrt(A^2 + 4 * k * (P0 * V0 / T0) * Ta)
  V <- Vfact1 * (Vfact2 - A)
  
  fact1 <- M
  fact2 <- k + (S^2) * (P0 * V0 / T0) * (Ta / (V^2))
  
  2 * pi * sqrt(fact1 / fact2)
}

piston_bounds <- list(
  M  = c(30,    60),    # piston weight (kg)
  S  = c(0.005, 0.020), # piston surface area (m^2)
  V0 = c(0.002, 0.010), # initial gas volume (m^3)
  k  = c(1000,  5000),  # spring coefficient (N/m)
  P0 = c(90000, 110000),# atmospheric pressure (N/m^2)
  Ta = c(290,   296),   # ambient temperature (K)
  T0 = c(340,   360)    # filling gas temperature (K)
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
n_global <- 80     # target number of global points
lNum     <- 25
theta    <- 0.3
nugget   <- TRUE

# Generate data
X_train <- matrix(runif(n_train * dim_), ncol = dim_)
y_train <- apply(X_train, 1, piston_unit)

X_test  <- matrix(runif(n_test * dim_),  ncol = dim_)
y_test  <- apply(X_test,  1, piston_unit)

y_mean <- mean(y_train)
y_sd   <- sd(y_train)
y_train_sc <- (y_train - y_mean) / y_sd

xy     <- cbind(X_train, y_train_sc)
x_test <- X_test

predIndices <- NULL  # filled in per-method below, after gIndices is known

# Helper: run twin() on matrix, returning smaller twin's indices
get_twin_indices <- function(data_mat, n_target, leaf_size = 8L) {
  r <- max(2, round(nrow(data_mat) / n_target))
  idx <- twin(data_mat, r = r, leaf_size = leaf_size)
  idx
}

# -----------------------------------------------------------------------------
# Method A: All variables — Use twin() directly on raw inputs
# -----------------------------------------------------------------------------
gIndices_all <- get_twin_indices(X_train, n_global)
cat(sprintf("[A] All variables: %d global points selected\n", length(gIndices_all)))

# -----------------------------------------------------------------------------
# Method B: Variable Selection — FIRST picks important variables first,
#                                 run twin() on resulting subspace
# -----------------------------------------------------------------------------
imp_B <- first(X_train, y_train_sc, n.knn = 2, rescale = TRUE, n.forward = 2, verbose = FALSE)
cat("[B] FIRST importances:", round(imp_B, 4), "\n")

selected_vars <- which(imp_B > 0)
if (length(selected_vars) == 0L) selected_vars <- seq_len(dim_)  # safety fallback
cat(sprintf("[B] Selected variables: %s\n", paste(selected_vars, collapse = ", ")))

gIndices_varsel <- get_twin_indices(X_train[, selected_vars, drop = FALSE], n_global)
cat(sprintf("[B] Variable selection: %d global points selected\n", length(gIndices_varsel)))

# -----------------------------------------------------------------------------
# Method C: Weights — FIRST importances used as continuous weights instead of 0/1;
#           columns rescaled by sqrt(importance) so twin()'s energy-distance metric
#           weighs less towards unimportant variables instead of dropping
# -----------------------------------------------------------------------------
imp_C <- imp_B  # reuse importance vector to get weights instead of selection
weights_C <- sqrt(imp_C / sum(imp_C))           # normalize so weights sum to 1 pre-sqrt
X_weighted <- sweep(scale(X_train), 2, weights_C, `*`)  # standardize, then weight

gIndices_weighted <- get_twin_indices(X_weighted, n_global)
cat(sprintf("[C] Weighted variables: %d global points selected\n", length(gIndices_weighted)))

# Common prediction-index set for lambda/nugget tuning
make_predIndices <- function(gIndices, n_pred = 30) {
  pool <- setdiff(seq_len(min(100, n_train)), gIndices)
  pool[seq_len(min(n_pred, length(pool)))]
}


# Run GLGP using glgp.R for each method and collect metrics
run_method <- function(label, gIndices) {
  predIdx <- make_predIndices(gIndices)
  
  cat(sprintf("\n--- Fitting GLGP [%s] (%d global pts) ---\n", label, length(gIndices)))
  t0 <- proc.time()
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
  elapsed <- (proc.time() - t0)["elapsed"]
  
  mu_pred    <- result$mu    * y_sd + y_mean
  sigma_pred <- result$sigma * y_sd
  
  errors    <- y_test - mu_pred
  rmse      <- sqrt(mean(errors^2))
  mae       <- mean(abs(errors))
  r_squared <- 1 - sum(errors^2) / sum((y_test - mean(y_test))^2)
  nlpd      <- mean(0.5 * (errors / sigma_pred)^2 + log(sigma_pred) + 0.5 * log(2 * pi))
  
  lower95    <- mu_pred - 1.96 * sigma_pred
  upper95    <- mu_pred + 1.96 * sigma_pred
  coverage95 <- mean(y_test >= lower95 & y_test <= upper95)
  
  list(
    label      = label,
    n_global   = length(gIndices),
    elapsed    = elapsed,
    rmse       = rmse,
    mae        = mae,
    r_squared  = r_squared,
    nlpd       = nlpd,
    coverage95 = coverage95,
    mu_pred    = mu_pred,
    sigma_pred = sigma_pred
  )
}

res_A <- run_method("All variables",      gIndices_all)
res_B <- run_method("Variable selection", gIndices_varsel)
res_C <- run_method("Weighted variables", gIndices_weighted)

# Summary table
summary_df <- data.frame(
  Method     = c(res_A$label, res_B$label, res_C$label),
  N_Global   = c(res_A$n_global, res_B$n_global, res_C$n_global),
  RMSE       = c(res_A$rmse, res_B$rmse, res_C$rmse),
  MAE        = c(res_A$mae, res_B$mae, res_C$mae),
  R_squared  = c(res_A$r_squared, res_B$r_squared, res_C$r_squared),
  NLPD       = c(res_A$nlpd, res_B$nlpd, res_C$nlpd),
  Coverage95 = c(res_A$coverage95, res_B$coverage95, res_C$coverage95) * 100,
  Fit_time_s = c(res_A$elapsed, res_B$elapsed, res_C$elapsed)
)

cat("\n===================== METHOD COMPARISON =====================\n")
print(summary_df, row.names = FALSE, digits = 4)
cat("===============================================================\n\n")

# Relative comparison against Method A (all variables) as the baseline
cat("Relative RMSE change vs. All-Variables baseline:\n")
cat(sprintf("  Variable selection : %+.1f%%\n",
            100 * (res_B$rmse - res_A$rmse) / res_A$rmse))
cat(sprintf("  Weighted variables : %+.1f%%\n\n",
            100 * (res_C$rmse - res_A$rmse) / res_A$rmse))


# Diagnostic plots: predicted vs actual for the three methods side by side
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

plot_results <- list(res_A, res_B, res_C)
lims <- range(c(y_test, res_A$mu_pred, res_B$mu_pred, res_C$mu_pred))

for (res in plot_results) {
  plot(y_test, res$mu_pred,
       pch = 19, cex = 0.6, col = "#2166ac88",
       xlim = lims, ylim = lims,
       xlab = "Actual", ylab = "Predicted",
       main = sprintf("%s\nRMSE=%.2f, NLPD=%.2f", res$label, res$rmse, res$nlpd))
  abline(0, 1, col = "red", lwd = 1.5)
}

par(mfrow = c(1, 1))
