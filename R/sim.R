# Synthetic telemetry: latent-factor Gaussian copula with diurnal / weekly
# seasonality and AR(1) persistence. One row = one host x 5-minute window.

# Loadings of each feature on the six latent factors:
#   A      business activity (diurnal + weekday cycle + AR(1))
#   W      host workload (AR(1))
#   R      user/host role (constant per host)
#   G_net  shared network-volume factor (AR(1))
#   G_auth shared authentication-volume factor (AR(1))
#   G_proc shared process-activity factor (AR(1))
feature_spec <- tibble::tribble(
  ~feature,                     ~plane,     ~A,    ~W,   ~R,   ~G_net, ~G_auth, ~G_proc,
  "auth_failed_logins",         "Identity",  0.30, 0.00, 0.00, 0.00,   0.55,    0.00,
  "auth_priv_escalations",      "Identity",  0.20, 0.00, 0.30, 0.00,   0.00,    0.00,
  "auth_unique_source_ips",     "Identity",  0.35, 0.00, 0.00, 0.00,   0.55,    0.00,
  "auth_off_hours_flag",        "Identity", -0.55, 0.00, 0.00, 0.00,   0.00,    0.00,
  "auth_token_refresh_rate",    "Identity",  0.40, 0.00, 0.00, 0.00,   0.15,    0.00,
  "auth_new_dest_hosts",        "Identity",  0.30, 0.00, 0.30, 0.00,   0.00,    0.00,
  "net_bytes_out_kb",           "Network",   0.50, 0.00, 0.00, 0.55,   0.00,    0.00,
  "net_packet_velocity",        "Network",   0.50, 0.00, 0.00, 0.60,   0.00,    0.00,
  "net_port_entropy",           "Network",   0.20, 0.00, 0.00, 0.20,   0.00,    0.00,
  "net_dns_query_rate",         "Network",   0.50, 0.20, 0.00, 0.00,   0.00,    0.00,
  "net_tls_duration_sec",       "Network",   0.30, 0.00, 0.00, 0.60,   0.00,    0.00,
  "net_syn_ack_ratio",          "Network",   0.00, 0.00, 0.00, 0.15,   0.00,    0.00,
  "net_smb_rdp_winrm_conns",    "Network",   0.40, 0.20, 0.30, 0.00,   0.00,    0.00,
  "host_cpu_utilization",       "Host",      0.30, 0.60, 0.00, 0.00,   0.00,    0.40,
  "host_memory_pressure",       "Host",      0.20, 0.50, 0.00, 0.00,   0.00,    0.00,
  "host_child_process_cnt",     "Host",      0.20, 0.30, 0.00, 0.00,   0.00,    0.50,
  "host_fs_write_bursts",       "Host",      0.20, 0.30, 0.00, 0.00,   0.00,    0.30,
  "host_active_threads",        "Host",      0.30, 0.60, 0.00, 0.00,   0.00,    0.40,
  "host_new_service_creations", "Host",      0.00, 0.20, 0.20, 0.00,   0.00,    0.00
)

# Target marginal distributions (quantile functions applied to copula uniforms)
feature_marginals <- list(
  auth_failed_logins         = \(u) qpois(u, 1.2),
  auth_priv_escalations      = \(u) qpois(u, 0.1),
  auth_unique_source_ips     = \(u) qpois(u, 2.0) + 1,
  auth_off_hours_flag        = \(u) as.numeric(u > 0.92),
  auth_token_refresh_rate    = \(u) qnorm(u, 12, 2.5),
  auth_new_dest_hosts        = \(u) qpois(u, 0.5),
  net_bytes_out_kb           = \(u) qlnorm(u, 7.5, 0.8),
  net_packet_velocity        = \(u) qnorm(u, 450, 60),
  net_port_entropy           = \(u) qunif(u, 0.15, 0.45),
  net_dns_query_rate         = \(u) qpois(u, 35),
  net_tls_duration_sec       = \(u) qexp(u, 1 / 45) + 5,
  net_syn_ack_ratio          = \(u) qnorm(u, 1.02, 0.05),
  net_smb_rdp_winrm_conns    = \(u) qpois(u, 1.5),
  host_cpu_utilization       = \(u) qbeta(u, 4, 6) * 100,
  host_memory_pressure       = \(u) qnorm(u, 62, 8),
  host_child_process_cnt     = \(u) qpois(u, 8),
  host_fs_write_bursts       = \(u) qpois(u, 14),
  host_active_threads        = \(u) qnorm(u, 120, 15),
  host_new_service_creations = \(u) qpois(u, 0.05)
)

feature_cols <- feature_spec$feature
latent_names <- c("A", "W", "R", "G_net", "G_auth", "G_proc")

# Stationary AR(1) with unit variance
ar1 <- function(n, phi) {
  e <- rnorm(n, sd = sqrt(1 - phi^2))
  as.numeric(stats::filter(e, phi, method = "recursive", init = rnorm(1)))
}

simulate_host <- function(host_id, sim, n_days = sim$n_days) {
  n <- n_days * 1440 / sim$window_minutes
  t0 <- as.POSIXct(sim$start_date, tz = "UTC")
  ts <- t0 + (seq_len(n) - 1) * sim$window_minutes * 60
  lt <- as.POSIXlt(ts)
  hour <- lt$hour + lt$min / 60
  weekday <- lt$wday %in% 1:5

  diurnal <- exp(-(hour - sim$diurnal_peak_hour)^2 / (2 * sim$diurnal_width_hours^2)) *
    ifelse(weekday, 1, sim$weekend_activity)
  a_raw <- sim$diurnal_amplitude * diurnal + ar1(n, sim$ar_phi)

  latent <- cbind(
    A      = (a_raw - mean(a_raw)) / sd(a_raw),
    W      = ar1(n, sim$ar_phi),
    R      = rep(rnorm(1), n),
    G_net  = ar1(n, sim$ar_phi),
    G_auth = ar1(n, sim$ar_phi),
    G_proc = ar1(n, sim$ar_phi)
  )
  L <- as.matrix(feature_spec[, latent_names])
  d <- nrow(L)
  idio_sd <- sqrt(1 - rowSums(L^2))
  Z <- latent %*% t(L) + matrix(rnorm(n * d), n, d) * rep(idio_sd, each = n)
  U <- pnorm(Z)
  colnames(U) <- feature_cols

  feats <- purrr::map(feature_cols, \(f) feature_marginals[[f]](U[, f])) |>
    rlang::set_names(feature_cols)

  tibble::tibble(
    host_id = host_id,
    window = seq_len(n),
    timestamp = ts,
    activity_index = latent[, "A"],
    !!!feats,
    label = "Normal",
    attack_type = "Normal (benign)",
    tier = NA_character_,
    campaign_id = NA_integer_
  )
}

simulate_hosts <- function(host_ids, sim, n_days = sim$n_days) {
  purrr::map(host_ids, \(h) simulate_host(h, sim, n_days)) |> dplyr::bind_rows()
}

# ---- Attack injection --------------------------------------------------------

attack_levels <- c(
  "Data Exfiltration", "Low-and-Slow Exfiltration", "Credential Stuffing",
  "Privilege Escalation", "Ransomware Staging", "Lateral Movement"
)
tier_levels <- c("low", "medium", "high")

# Features each scenario modifies (used to score attribution)
injected_features <- list(
  "Data Exfiltration"         = c("net_bytes_out_kb", "net_tls_duration_sec", "net_port_entropy"),
  "Low-and-Slow Exfiltration" = c("net_bytes_out_kb", "net_tls_duration_sec", "net_port_entropy"),
  "Credential Stuffing"       = c("auth_unique_source_ips", "auth_failed_logins"),
  "Privilege Escalation"      = c("auth_priv_escalations"),
  "Ransomware Staging"        = c("host_child_process_cnt", "host_fs_write_bursts", "net_dns_query_rate"),
  "Lateral Movement"          = c("auth_new_dest_hosts", "net_smb_rdp_winrm_conns", "host_new_service_creations")
)

runif_range <- function(n, r) runif(n, r[[1]], r[[2]])
add_counts <- function(n, add, pois) add + if (pois > 0) rpois(n, pois) else 0

inject_attack <- function(df, type, tier, spec) {
  n <- nrow(df)
  s <- spec[[type]][[tier]]
  out <- switch(type,
    data_exfiltration = dplyr::mutate(df,
      net_bytes_out_kb     = net_bytes_out_kb * runif_range(n, s$bytes_mult),
      net_tls_duration_sec = net_tls_duration_sec * runif_range(n, s$tls_mult),
      net_port_entropy     = net_port_entropy * runif_range(n, s$entropy_mult)
    ),
    credential_stuffing = dplyr::mutate(df,
      auth_unique_source_ips = auth_unique_source_ips + add_counts(n, s$ips_add, s$ips_pois),
      auth_failed_logins     = auth_failed_logins + add_counts(n, s$fails_add, s$fails_pois)
    ),
    privilege_escalation = dplyr::mutate(df,
      auth_priv_escalations = auth_priv_escalations + s$priv_add
    ),
    ransomware_staging = dplyr::mutate(df,
      host_child_process_cnt = host_child_process_cnt + add_counts(n, s$child_add, s$child_pois),
      host_fs_write_bursts   = host_fs_write_bursts + add_counts(n, s$fs_add, s$fs_pois),
      net_dns_query_rate     = net_dns_query_rate * runif_range(n, s$dns_mult)
    ),
    lateral_movement = dplyr::mutate(df,
      auth_new_dest_hosts        = auth_new_dest_hosts + add_counts(n, s$dest_add, s$dest_pois),
      net_smb_rdp_winrm_conns    = net_smb_rdp_winrm_conns + add_counts(n, s$smb_add, s$smb_pois),
      host_new_service_creations = host_new_service_creations + s$svc_add + rbinom(n, 1, s$svc_prob)
    )
  )
  dplyr::mutate(out, label = "Anomaly", attack_type = spec[[type]]$label, tier = .env$tier)
}

# Single-window attacks: every type x tier cell, injected into sampled normal windows
simulate_window_attacks <- function(pool, spec, n_per_cell) {
  cells <- tidyr::expand_grid(type = names(spec), tier = tier_levels)
  idx <- sample.int(nrow(pool), nrow(cells) * n_per_cell)
  purrr::map(seq_len(nrow(cells)), \(i) {
    rows <- idx[((i - 1) * n_per_cell + 1):(i * n_per_cell)]
    inject_attack(pool[rows, ], cells$type[i], cells$tier[i], spec)
  }) |>
    dplyr::bind_rows()
}

# Low-and-slow exfiltration: 1-4 hour campaigns of consecutive modest windows
simulate_low_and_slow <- function(pool, ls) {
  max_len <- ls$length[[2]]
  slots <- pool |>
    dplyr::mutate(slot = (window - 1) %/% max_len) |>
    dplyr::distinct(host_id, slot) |>
    dplyr::slice_sample(n = ls$n_campaigns)
  lens <- sample(ls$length[[1]]:ls$length[[2]], ls$n_campaigns, replace = TRUE)

  purrr::map(seq_len(ls$n_campaigns), \(i) {
    start <- slots$slot[i] * max_len + 1
    df <- dplyr::filter(pool, host_id == slots$host_id[i], window >= start, window < start + lens[i])
    n <- nrow(df)
    dplyr::mutate(df,
      net_bytes_out_kb     = net_bytes_out_kb * runif_range(n, ls$bytes_mult),
      net_tls_duration_sec = net_tls_duration_sec * runif_range(n, ls$tls_mult),
      net_port_entropy     = net_port_entropy * runif_range(n, ls$entropy_mult),
      label = "Anomaly", attack_type = ls$label, tier = "low", campaign_id = i
    )
  }) |>
    dplyr::bind_rows()
}

# Full synthetic study for one seed
simulate_study <- function(cfg) {
  sim <- cfg$simulation
  ids <- cumsum(c(sim$hosts_train, sim$hosts_val, sim$hosts_test, sim$hosts_attack_pool))
  train <- simulate_hosts(1:ids[1], sim)
  val   <- simulate_hosts((ids[1] + 1):ids[2], sim)
  test_normal <- simulate_hosts((ids[2] + 1):ids[3], sim)
  pool_ids <- (ids[3] + 1):ids[4]
  pool <- simulate_hosts(pool_ids, sim)
  half <- seq_len(length(pool_ids) / 2)

  window_attacks <- simulate_window_attacks(
    dplyr::filter(pool, host_id %in% pool_ids[half]), cfg$attacks, sim$n_attacks_per_cell
  )
  low_slow <- simulate_low_and_slow(
    dplyr::filter(pool, host_id %in% pool_ids[-half]), cfg$low_and_slow
  )

  test <- dplyr::bind_rows(test_normal, window_attacks, low_slow) |>
    dplyr::mutate(
      attack_type = factor(attack_type, levels = c("Normal (benign)", attack_levels)),
      tier = factor(tier, levels = tier_levels)
    )
  list(train = train, val = val, test = test)
}
