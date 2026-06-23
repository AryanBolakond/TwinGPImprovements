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

source("C:/Users/Aryan/Downloads/RoshanJosephResearch/glgp.R")


# Borehole function implementation and bounds
# (credit to sfu.ca/~ssurjano/borehole.html)
borehole <- function(xx) {
  rw <- xx[1]; r  <- xx[2]; Tu <- xx[3]; Hu <- xx[4]
  Tl <- xx[5]; Hl <- xx[6]; L  <- xx[7]; Kw <- xx[8]
  
  frac1  <- 2 * pi * Tu * (Hu - Hl)
  frac2a <- 2 * L * Tu / (log(r / rw) * rw^2 * Kw)
  frac2b <- Tu / Tl
  frac2  <- log(r / rw) * (1 + frac2a + frac2b)
  
  frac1 / frac2
}

# 
borehole_bounds <- list(
  rw = c(0.05,  0.15), r  = c(100,   50000), Tu = c(63070, 115600),
  Hu = c(990,   1110), Tl = c(63.1,  116),   Hl = c(700,   820),
  L  = c(1120,  1680), Kw = c(9855,  12045)
)

unit_to_physical <- function(u) {
  lo <- sapply(borehole_bounds, `[`, 1)
  hi <- sapply(borehole_bounds, `[`, 2)
  lo + u * (hi - lo)
}
borehole_unit <- function(u) borehole(unit_to_physical(u))


# Simulation settings
set.seed(42)

n_train  <- 300
n_test   <- 150
dim_     <- 8
n_global <- 80     # target number of global points
lNum     <- 25
theta    <- 0.3
nugget   <- TRUE

# Generate data
X_train <- matrix(runif(n_train * dim_), ncol = dim_)
y_train <- apply(X_train, 1, borehole_unit)

X_test  <- matrix(runif(n_test * dim_),  ncol = dim_)
y_test  <- apply(X_test,  1, borehole_unit)

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
make_predIndices <- function(gIndices, n_pred = 24) {
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