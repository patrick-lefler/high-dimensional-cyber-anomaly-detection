# Peaks-over-threshold helpers.

# GPD quantile threshold for exceedance probability p = 1 - q per window.
# Uses the exponential limit when xi is effectively zero.
gpd_tau <- function(u, sigma, xi, N, Nu, p) {
  if (abs(xi) < 1e-6) {
    u + sigma * log(Nu / (N * p))
  } else {
    u + (sigma / xi) * (((N / Nu) * p)^(-xi) - 1)
  }
}

fit_pot <- function(losses, u_prob = 0.95, std_err = TRUE) {
  u <- unname(quantile(losses, u_prob))
  fit <- evd::fpot(losses, threshold = u, model = "gpd", std.err = std_err)
  list(
    u = u,
    sigma = unname(fit$estimate["scale"]),
    xi = unname(fit$estimate["shape"]),
    se = if (std_err) fit$std.err else c(scale = NA_real_, shape = NA_real_),
    N = length(losses),
    Nu = sum(losses > u),
    fit = fit
  )
}

pot_taus <- function(pot, q) {
  purrr::map_dbl(q, \(qq) gpd_tau(pot$u, pot$sigma, pot$xi, pot$N, pot$Nu, 1 - qq))
}

# Block bootstrap: resample whole blocks (host-days) to respect autocorrelation
bootstrap_taus <- function(losses, block, q, B, u_prob) {
  blocks <- split(losses, block)
  purrr::map(seq_len(B), \(b) {
    x <- unlist(blocks[sample.int(length(blocks), replace = TRUE)], use.names = FALSE)
    pot <- tryCatch(fit_pot(x, u_prob, std_err = FALSE), error = \(e) NULL)
    if (is.null(pot)) return(NULL)
    tibble::tibble(b = b, q = q, tau = pot_taus(pot, q), xi = pot$xi, sigma = pot$sigma)
  }) |>
    dplyr::bind_rows()
}

empirical_quantile <- function(x, q, min_exceed = 10) {
  purrr::map_dbl(q, \(qq) {
    if (length(x) * (1 - qq) < min_exceed) NA_real_ else unname(quantile(x, qq))
  })
}
