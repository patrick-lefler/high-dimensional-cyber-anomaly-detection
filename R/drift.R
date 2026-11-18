# Drift metrics on preprocessed features (units = training standard deviations).

psi <- function(ref, new, bins = 10) {
  vals <- sort(unique(ref))
  br <- if (length(vals) <= bins) {
    c(-Inf, head(vals, -1) + diff(vals) / 2, Inf)
  } else {
    b <- unique(quantile(ref, seq(0, 1, length.out = bins + 1)))
    c(-Inf, b[-c(1, length(b))], Inf)
  }
  p <- tabulate(cut(ref, br, labels = FALSE), length(br) - 1) / length(ref)
  q <- tabulate(cut(new, br, labels = FALSE), length(br) - 1) / length(new)
  p <- pmax(p, 1e-6)
  q <- pmax(q, 1e-6)
  sum((q - p) * log(q / p))
}

wasserstein1 <- function(ref, new, grid = (1:999) / 1000) {
  mean(abs(quantile(ref, grid, names = FALSE) - quantile(new, grid, names = FALSE)))
}

drift_metrics <- function(X_ref, X_new) {
  purrr::map(colnames(X_ref), \(f) {
    tibble::tibble(
      feature = f,
      psi = psi(X_ref[, f], X_new[, f]),
      wasserstein = wasserstein1(X_ref[, f], X_new[, f]),
      ks_p = suppressWarnings(ks.test(X_ref[, f], X_new[, f])$p.value)
    )
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(ks_p_holm = p.adjust(ks_p, method = "holm"))
}
