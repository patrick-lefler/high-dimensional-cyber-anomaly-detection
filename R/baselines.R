# Simple detectors scored on the same preprocessed matrices as the autoencoder.

fit_baselines <- function(X_train, k = 4, seed = 1) {
  list(
    pca = prcomp(X_train, center = TRUE, scale. = FALSE, rank. = k),
    mu = colMeans(X_train),
    S = cov(X_train),
    iso = isotree::isolation.forest(X_train, ntrees = 200, nthreads = 1, seed = seed)
  )
}

pca_recon_mse <- function(pca, X) {
  scores <- scale(X, center = pca$center, scale = FALSE) %*% pca$rotation
  recon <- scores %*% t(pca$rotation)
  rowMeans((scale(X, center = pca$center, scale = FALSE) - recon)^2)
}

score_baselines <- function(fit, X) {
  tibble::tibble(
    `Sum of squared z-scores` = rowSums(X^2),
    `Max |z|` = do.call(pmax, as.data.frame(abs(X))),
    `PCA reconstruction (k = 4)` = pca_recon_mse(fit$pca, X),
    `Mahalanobis distance` = mahalanobis(X, fit$mu, fit$S),
    `Isolation Forest` = predict(fit$iso, X)
  )
}

# Rank-based ROC AUC (Mann-Whitney); ties count one half
auc_rank <- function(score, is_attack) {
  r <- rank(score)
  n1 <- sum(is_attack)
  n0 <- sum(!is_attack)
  (sum(r[is_attack]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

pr_auc_vec <- function(score, is_attack) {
  truth <- factor(ifelse(is_attack, "Anomaly", "Normal"), levels = c("Anomaly", "Normal"))
  yardstick::pr_auc_vec(truth, score, event_level = "first")
}
