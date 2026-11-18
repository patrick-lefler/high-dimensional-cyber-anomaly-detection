# Regenerates every result the pages read. Run from the project root:
#   Rscript run_all.R
# Takes a few minutes; results are written to results/*.rds.

suppressPackageStartupMessages(library(tidyverse))
purrr::walk(list.files("R", pattern = "\\.R$", full.names = TRUE), source)
cfg <- yaml::read_yaml("params.yaml")
dir.create("results", showWarnings = FALSE)
seed <- cfg$seeds$main
d <- length(feature_cols)
q <- unlist(cfg$evt$q)
timer <- Sys.time()
say <- \(...) message(sprintf("[%5.0fs] ", as.numeric(Sys.time() - timer, units = "secs")), ...)

# ---- Main run ----------------------------------------------------------------
say("main run")
main <- run_pipeline(seed, cfg, keep = TRUE)

# ---- Gradient check ----------------------------------------------------------
say("gradient check")
set.seed(seed)
batch <- main$pp$X_train[sample.int(nrow(main$pp$X_train), 32), ]
grad_check <- bind_rows(
  gradient_check(ae_init(c(d, cfg$autoencoder$bottleneck, d)), batch) |>
    mutate(model = "Shallow (d-4-d)"),
  gradient_check(ae_init(c(d, cfg$autoencoder$deep_hidden, cfg$autoencoder$bottleneck,
                           cfg$autoencoder$deep_hidden, d)), batch) |>
    mutate(model = "Deeper (d-8-4-8-d)")
)

# ---- Architecture sensitivity --------------------------------------------------
say("architecture sweep")
a <- cfg$autoencoder
is_att <- main$test$label == "Anomaly"
archs <- c(
  purrr::map(unlist(a$sweep_bottleneck), \(k) list(label = sprintf("Shallow d-%d-d", k), dims = c(d, k, d))),
  list(list(label = sprintf("Deeper d-%d-%d-%d-d", a$deep_hidden, a$bottleneck, a$deep_hidden),
            dims = c(d, a$deep_hidden, a$bottleneck, a$deep_hidden, d)))
)
arch_sweep <- purrr::imap(archs, \(arch, i) {
  set.seed(seed + i)
  m <- train_ae(main$pp$X_train, main$pp$X_val, arch$dims, a$max_epochs, a$patience, a$batch_size, a$lr)
  tibble(
    model = arch$label,
    train_mse = m$train_losses[m$best_epoch],
    val_mse = m$best_val_loss,
    best_epoch = m$best_epoch,
    roc_auc = auc_rank(ae_reconstruct(m, main$pp$X_test)$sample_mse, is_att)
  )
}) |>
  bind_rows() |>
  bind_rows(purrr::map(unlist(a$sweep_bottleneck), \(k) {
    pca <- prcomp(main$pp$X_train, rank. = k)
    tibble(
      model = sprintf("PCA, k = %d (linear reference)", k),
      train_mse = mean(pca_recon_mse(pca, main$pp$X_train)),
      val_mse = mean(pca_recon_mse(pca, main$pp$X_val)),
      best_epoch = NA_real_,
      roc_auc = auc_rank(pca_recon_mse(pca, main$pp$X_test), is_att)
    )
  }) |> bind_rows())

# ---- EVT: bootstrap, calibration, declustering, tail composition --------------
say("EVT bootstrap")
val_block <- main$study$val |> mutate(block = paste(host_id, (window - 1) %/% 288))
set.seed(seed)
boot <- bootstrap_taus(main$loss_val, val_block$block, q, cfg$evt$bootstrap_B, cfg$evt$u_prob)

test_normal_loss <- main$test$recon_loss[main$test$label == "Normal"]
calibration <- tibble(q = q, tau = main$taus) |>
  mutate(
    n_normal = length(test_normal_loss),
    expected_fp = (1 - q) * n_normal,
    observed_fp = map_int(tau, \(t) sum(test_normal_loss > t)),
    p_value = map2_dbl(observed_fp, q, \(x, qq) binom.test(x, length(test_normal_loss), 1 - qq)$p.value)
  )

say("declustering")
r <- cfg$evt$decluster_run_length
decluster <- list(
  theta = evd::exi(main$loss_val, main$pot$u, r = r),
  fit = evd::fpot(main$loss_val, main$pot$u, cmax = TRUE, r = r),
  n_exceed = main$pot$Nu
)
decluster$n_clusters <- decluster$fit$nhigh

tail_rows <- which(main$loss_val > main$taus[["q0.99"]])
val_resid <- ae_reconstruct(main$model, main$pp$X_val)$residuals
tail_composition <- tibble(feature = feature_cols[max.col(val_resid[tail_rows, , drop = FALSE], ties.method = "first")]) |>
  count(feature, name = "n") |>
  mutate(share = n / sum(n)) |>
  arrange(desc(n))

# ---- Attribution: top-3 hit rate --------------------------------------------
say("attribution")
top3_hit <- function(M, types) {
  top3 <- t(apply(M, 1, \(r) feature_cols[order(r, decreasing = TRUE)[1:3]]))
  map_lgl(seq_along(types), \(i) any(top3[i, ] %in% injected_features[[types[i]]]))
}
detected <- which(main$test$label == "Anomaly" & main$test$flag_999)
attribution <- tibble(
  attack_type = as.character(main$test$attack_type[detected]),
  ae_hit = top3_hit(main$residuals[detected, , drop = FALSE], attack_type),
  z_hit = top3_hit(abs(main$pp$X_test[detected, , drop = FALSE]), attack_type)
) |>
  summarise(n = n(), ae_top3 = mean(ae_hit), z_top3 = mean(z_hit), .by = attack_type)

# ---- Drift demonstration -----------------------------------------------------
say("drift demo")
dcfg <- cfg$drift
set.seed(seed + 100)
week <- simulate_hosts(1000 + seq_len(dcfg$week_hosts), cfg$simulation, n_days = dcfg$week_days)
X_week <- as.matrix(recipes::bake(main$pp$recipe, new_data = week[feature_cols]))
X_ref <- rbind(main$pp$X_train, main$pp$X_val)
shift <- \(X, s) { X[, unlist(dcfg$shifted_features)] <- X[, unlist(dcfg$shifted_features)] + s; X }
drift <- bind_rows(
  drift_metrics(X_ref, X_week) |> mutate(scenario = "No shift"),
  drift_metrics(X_ref, shift(X_week, dcfg$harmless_shift_sd)) |>
    mutate(scenario = sprintf("Small shift (%.2f SD)", dcfg$harmless_shift_sd)),
  drift_metrics(X_ref, shift(X_week, dcfg$material_shift_sd)) |>
    mutate(scenario = sprintf("Material shift (%.1f SD)", dcfg$material_shift_sd))
) |>
  mutate(shifted = feature %in% unlist(dcfg$shifted_features))
week_loss <- ae_reconstruct(main$model, X_week)$sample_mse
alert_rate <- tibble(
  scenario = unique(drift$scenario),
  alert_rate = c(
    mean(week_loss > main$taus[["q0.999"]]),
    mean(ae_reconstruct(main$model, shift(X_week, dcfg$harmless_shift_sd))$sample_mse > main$taus[["q0.999"]]),
    mean(ae_reconstruct(main$model, shift(X_week, dcfg$material_shift_sd))$sample_mse > main$taus[["q0.999"]])
  ),
  n = nrow(X_week)
)

# ---- Service graph ----------------------------------------------------------
say("service graph")
graph <- analyse_service_graph(build_service_graph(cfg$graph$seed), cfg$graph)

# Containment routing: map each alerted test window to the service its host supports
set.seed(seed + 200)
routing <- main$test |>
  mutate(service = sample(graph$table$name, n(), replace = TRUE)) |>
  filter(flag_999) |>
  left_join(select(graph$table, name, criticality_tier), by = join_by(service == name)) |>
  count(criticality_tier, label)

# ---- Multi-seed robustness --------------------------------------------------
say("multi-seed robustness")
future::plan(future::multisession, workers = min(6, future::availableCores() - 1))
robust <- furrr::future_map(
  seq_len(cfg$seeds$robustness),
  \(s) {
    purrr::walk(list.files("R", pattern = "\\.R$", full.names = TRUE), source)
    run_pipeline(s, cfg)
  },
  .options = furrr::furrr_options(seed = TRUE, packages = c("dplyr", "tibble", "purrr"))
) |>
  bind_rows()
future::plan(future::sequential)

# ---- Save ---------------------------------------------------------------------
say("saving")
main_slim <- main
main_slim$study$train <- main$study$train |> filter(host_id == 1)
saveRDS(main_slim, "results/main.rds")
saveRDS(main$study$train |> select(all_of(feature_cols)) |> cor(method = "spearman"), "results/train_cor.rds")
saveRDS(list(
  grad_check = grad_check, arch_sweep = arch_sweep, boot = boot, calibration = calibration,
  decluster = decluster, tail_composition = tail_composition, attribution = attribution,
  drift = drift, alert_rate = alert_rate, routing = routing
), "results/analysis.rds")
saveRDS(graph, "results/graph.rds")
saveRDS(robust, "results/robust.rds")
say("done")
