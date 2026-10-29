# =============================================================================
# Tutor competence as the pathway from learning service quality to 4C skills:
# Evidence from urban, coastal, and highland community learning centers in
# Indonesia  (F1000Research, manuscript 193003)
#
# Open-source PLS-SEM analysis script (R + seminr).
# The published estimates were produced with SmartPLS 4; this script
# re-implements the same model and procedures with the open-source R package
# seminr so that the analysis can be inspected and reproduced.
#
# Authors : Ahlun Ansar, Fauzan Hidayatullah, Savira Widya Puspitasari,
#           Latang, Fatmawati Gaffar (Universitas Negeri Makassar)
# License : MIT (see LICENSE)
#
# Usage   : Rscript R/analysis.R [path/to/data.csv]
#           Without an argument the script runs on the SIMULATED example data
#           in data/example_data_simulated.csv. The real dataset is available
#           from the corresponding author under a data-use agreement (see the
#           Data availability statement of the article).
# =============================================================================

# ---- 0. Setup ---------------------------------------------------------------
required <- c("seminr")
missing  <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) install.packages(missing, repos = "https://cloud.r-project.org")
library(seminr)

args      <- commandArgs(trailingOnly = TRUE)
data_path <- if (length(args) >= 1) args[1] else "data/example_data_simulated.csv"
out_dir   <- "output"
dir.create(out_dir, showWarnings = FALSE)

SEED   <- 2026      # random seed for bootstrapping and cross-validation
NBOOT  <- 5000      # bootstrap resamples (as in the article)
NFOLDS <- 10        # folds for PLSpredict
NREPS  <- 10        # repetitions for PLSpredict

dat <- read.csv(data_path, check.names = FALSE)
cat("Data:", data_path, "-", nrow(dat), "respondents,", ncol(dat), "items\n")

# Item codes in questionnaire order (1-10 semantic differential scale)
#   X1-X24 learning service quality (Instructional Skills Questionnaire, 8 dimensions x 3)
#   Z1-Z20 tutor competence (professional Z1-Z4, pedagogical Z5-Z10,
#          social Z11-Z14, personality Z15-Z20)
#   Y1-Y15 4C skills (collaboration Y1-Y4, critical thinking Y5-Y8,
#          creativity Y9-Y11, communication Y12-Y15)
all_items <- list(
  X = paste0("X", 1:24),
  Z = paste0("Z", 1:20),
  Y = paste0("Y", 1:15)
)

# Indicators retained after item purification (Table 1 of the article)
retained <- list(
  X = paste0("X", c(2, 4, 6, 7, 9, 12, 13, 14, 15, 16, 17, 19, 21, 22, 23)),
  Z = paste0("Z", c(2, 3, 4, 11, 12, 13, 15:20)),
  Y = paste0("Y", c(3, 5:8, 9, 10, 13))
)
stopifnot(all(unlist(all_items) %in% names(dat)))

# ---- 1. Model specification -------------------------------------------------
build_mm <- function(items) {
  constructs(
    composite("LSQ",   items$X, weights = mode_A),   # learning service quality (X)
    composite("TC",    items$Z, weights = mode_A),   # tutor competence (Z)
    composite("FourC", items$Y, weights = mode_A)    # 4C skills (Y)
  )
}

sm <- relationships(
  paths(from = "LSQ", to = c("TC", "FourC")),   # H2 (X -> Z), H1 (X -> Y)
  paths(from = "TC",  to = "FourC")             # H3 (Z -> Y)
)

# ---- 2. Initial model with all 59 items (basis for item purification) ------
pls_initial <- estimate_pls(data = dat, measurement_model = build_mm(all_items),
                            structural_model = sm)
sum_initial <- summary(pls_initial)
write.csv(round(sum_initial$loadings, 3), file.path(out_dir, "01_loadings_initial_59items.csv"))

# Indicators were removed iteratively from the weakest loading following
# Hair et al. (2019, 2022): loading >= 0.708 retained; < 0.40 removed;
# 0.40-0.708 removed only if deletion raised rho_C or AVE to the thresholds.
# The final retained set is given in `retained` above.

# ---- 3. Final measurement model (35 indicators) -----------------------------
pls_model <- estimate_pls(data = dat, measurement_model = build_mm(retained),
                          structural_model = sm)
s <- summary(pls_model)

# Outer loadings (only the loading of each item on its own construct)
load_mat <- s$loadings
own_load <- data.frame(
  construct = rep(c("LSQ", "TC", "FourC"), lengths(retained)),
  item      = unlist(retained),
  loading   = mapply(function(i, c) load_mat[i, c], unlist(retained),
                     rep(c("LSQ", "TC", "FourC"), lengths(retained)))
)
write.csv(own_load, file.path(out_dir, "02_outer_loadings.csv"), row.names = FALSE)

# Reliability and convergent validity (Table 3): alpha, rhoA, rhoC, AVE
reliab <- s$reliability
print(reliab)
write.csv(round(unclass(reliab), 3), file.path(out_dir, "03_reliability_validity.csv"))

# Discriminant validity (Table 4): Fornell-Larcker criterion and HTMT
write.csv(round(unclass(s$validity$fl_criteria), 3), file.path(out_dir, "04a_fornell_larcker.csv"))
write.csv(round(unclass(s$validity$htmt), 3),        file.path(out_dir, "04b_htmt.csv"))
print(s$validity$fl_criteria); print(s$validity$htmt)

# ---- 4. Common method bias: full collinearity VIF (Kock, 2015) ---------------
# Each construct score is regressed on all other construct scores;
# VIF = 1 / (1 - R^2). Values > 3.3 indicate potential contamination.
scores <- as.data.frame(pls_model$construct_scores)
full_vif <- sapply(names(scores), function(v) {
  r2 <- summary(lm(scores[[v]] ~ ., data = scores[setdiff(names(scores), v)]))$r.squared
  1 / (1 - r2)
})
print(round(full_vif, 3))
write.csv(data.frame(construct = names(full_vif), full_collinearity_VIF = round(full_vif, 3)),
          file.path(out_dir, "05_full_collinearity_vif.csv"), row.names = FALSE)

# ---- 5. Structural model: VIF, R2, f2 (Table 5 and Table 6) ------------------
inner_vif <- s$vif_antecedents
print(inner_vif)
capture.output(inner_vif, file = file.path(out_dir, "06a_inner_vif.txt"))
write.csv(round(unclass(s$paths), 3),   file.path(out_dir, "06b_paths_r2.csv"))
write.csv(round(unclass(s$fSquare), 3), file.path(out_dir, "06c_f_squared.csv"))

# ---- 6. Bootstrapping (5,000 resamples) -------------------------------------
boot_model <- bootstrap_model(seminr_model = pls_model, nboot = NBOOT,
                              cores = NULL, seed = SEED)
sb <- summary(boot_model, alpha = 0.05)

bp <- as.data.frame(unclass(sb$bootstrapped_paths))
path_tab <- data.frame(
  path = rownames(bp),
  beta = bp[, "Original Est."],
  SD   = bp[, "Bootstrap SD"],
  t    = bp[, "T Stat."]
)
path_tab$p        <- 2 * pnorm(-abs(path_tab$t))
path_tab$CI_lower <- path_tab$beta - 1.96 * path_tab$SD   # beta +/- 1.96 SE (as reported)
path_tab$CI_upper <- path_tab$beta + 1.96 * path_tab$SD
path_tab$CI_perc_2.5  <- bp[, "2.5% CI"]                  # percentile CI, for reference
path_tab$CI_perc_97.5 <- bp[, "97.5% CI"]
print(path_tab, digits = 3)
write.csv(path_tab, file.path(out_dir, "07_hypothesis_tests_direct.csv"), row.names = FALSE)

# ---- 7. Mediation: indirect effect X -> Z -> Y (H4) -------------------------
ind <- specific_effect_significance(boot_model, from = "LSQ", through = "TC",
                                    to = "FourC", alpha = 0.05)
ind_beta <- unname(ind["Original Est."]); ind_sd <- unname(ind["Bootstrap SD"])
ind_t    <- ind_beta / ind_sd
a <- pls_model$path_coef["LSQ", "TC"]; b <- pls_model$path_coef["TC", "FourC"]
upsilon <- (a^2) * (b^2)                                     # Lachowicz et al. (2018)
ind_tab <- data.frame(
  path = "LSQ -> TC -> FourC", beta = ind_beta, SD = ind_sd, t = ind_t,
  p = 2 * pnorm(-abs(ind_t)),
  CI_lower = ind_beta - 1.96 * ind_sd, CI_upper = ind_beta + 1.96 * ind_sd,
  CI_perc_2.5 = unname(ind["2.5% CI"]), CI_perc_97.5 = unname(ind["97.5% CI"]),
  upsilon = upsilon
)
print(ind_tab, digits = 3)
write.csv(ind_tab, file.path(out_dir, "08_indirect_effect.csv"), row.names = FALSE)

# Mediation type (Zhao et al., 2010; Nitzl et al., 2016)
direct_p <- path_tab$p[path_tab$path == "LSQ  ->  FourC" | grepl("LSQ.*FourC", path_tab$path)][1]
med_type <- if (ind_tab$p >= 0.05) "no mediation" else if (direct_p >= 0.05) "indirect-only (full) mediation" else "complementary/competitive (partial) mediation"
cat("Mediation type:", med_type, "\n")

# ---- 8. Out-of-sample predictive relevance: Q2_predict (PLSpredict) ---------
# k-fold cross-validation repeated NREPS times. For each endogenous construct,
# composite scores predicted from the exogenous indicators (direct antecedent
# approach) are compared with the observed composite scores of the holdout
# fold; the benchmark is the training-sample mean (Shmueli et al., 2019).
set.seed(SEED)
endog   <- c("TC", "FourC")
mm_fin  <- build_mm(retained)
sse_pls <- sse_naive <- setNames(numeric(length(endog)), endog)
abs_pls <- setNames(numeric(length(endog)), endog); n_obs <- 0
for (rep_i in seq_len(NREPS)) {
  folds <- sample(rep(seq_len(NFOLDS), length.out = nrow(dat)))
  for (k in seq_len(NFOLDS)) {
    train <- dat[folds != k, ]; test <- dat[folds == k, ]
    m_tr  <- estimate_pls(train, mm_fin, sm)
    pred  <- predict(m_tr, testData = test, technique = predict_DA)
    # observed composite scores of the test fold using training weights,
    # standardized with training means/SDs
    obs_scores <- function(cn) {
      items <- retained[[c(TC = "Z", FourC = "Y")[[cn]]]]
      z <- scale(test[, items], center = m_tr$meanData[items], scale = m_tr$sdData[items])
      as.vector(z %*% m_tr$outer_weights[items, cn])
    }
    for (cn in endog) {
      obs <- obs_scores(cn)
      prd <- pred$predicted_composite_scores[, cn]
      sse_pls[cn]   <- sse_pls[cn]   + sum((obs - prd)^2)
      sse_naive[cn] <- sse_naive[cn] + sum((obs - 0)^2)  # training mean of std. scores = 0
      abs_pls[cn]   <- abs_pls[cn]   + sum(abs(obs - prd))
    }
    n_obs <- n_obs + nrow(test)
  }
}
q2_tab <- data.frame(
  construct  = endog,
  R2         = s$paths["R^2", endog],
  adj_R2     = s$paths["AdjR^2", endog],
  Q2_predict = 1 - sse_pls / sse_naive,
  RMSE       = sqrt(sse_pls / n_obs),
  MAE        = abs_pls / n_obs
)
print(q2_tab, digits = 3)
write.csv(q2_tab, file.path(out_dir, "09_explanatory_predictive_power.csv"), row.names = FALSE)

# ---- 9. Session information -------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\nDone. Results written to the 'output' folder.\n")
