# One full run: simulate -> preprocess -> train -> EVT on held-out losses ->
# evaluate against the attack set and the baseline detectors.

preprocess <- function(study) {
  rec <- recipes::recipe(~ ., data = study$train[feature_cols]) |>
    recipes::step_YeoJohnson(recipes::all_numeric_predictors()) |>
    recipes::step_normalize(recipes::all_numeric_predictors()) |>
    recipes::prep()
  bake_m <- \(df) as.matrix(recipes::bake(rec, new_data = df[feature_cols]))
  list(recipe = rec, X_train = bake_m(study$train), X_val = bake_m(study$val), X_test = bake_m(study$test))
}

ae_label <- "Autoencoder (shallow, k = 4)"

# Detection share by class and tier for a logical flag vector
detection_by_cell <- function(test, flag) {
  test |>
    dplyr::mutate(flag = flag) |>
    dplyr::filter(label == "Anomaly") |>
    dplyr::summarise(n = dplyr::n(), detected = sum(flag), .by = c(attack_type, tier))
}

# Chance that a run of normal windows of the same lengths would contain an alert
null_campaign_rate <- function(test, flag, lengths, n_draws = 2000) {
  normal <- test |>
    dplyr::mutate(flag = flag) |>
    dplyr::filter(label == "Normal") |>
    dplyr::arrange(host_id, window)
  by_host <- split(normal$flag, normal$host_id)
  hits <- purrr::map_lgl(seq_len(n_draws), \(i) {
    h <- by_host[[sample.int(length(by_host), 1)]]
    L <- sample(lengths, 1)
    s <- sample.int(length(h) - L + 1, 1)
    any(h[s:(s + L - 1)])
  })
  mean(hits)
}

run_pipeline <- function(seed, cfg, keep = FALSE) {
  set.seed(seed)
  study <- simulate_study(cfg)
  pp <- preprocess(study)
  d <- length(feature_cols)
  a <- cfg$autoencoder
  model <- train_ae(pp$X_train, pp$X_val, c(d, a$bottleneck, d),
                    a$max_epochs, a$patience, a$batch_size, a$lr)

  loss_train <- ae_reconstruct(model, pp$X_train)$sample_mse
  loss_val <- ae_reconstruct(model, pp$X_val)$sample_mse
  recon_test <- ae_reconstruct(model, pp$X_test)

  q <- unlist(cfg$evt$q)
  pot <- fit_pot(loss_val, cfg$evt$u_prob)
  pot_train <- fit_pot(loss_train, cfg$evt$u_prob)
  taus <- rlang::set_names(pot_taus(pot, q), paste0("q", q))

  test <- study$test |>
    dplyr::mutate(
      recon_loss = recon_test$sample_mse,
      flag_990 = recon_loss > taus[["q0.99"]],
      flag_999 = recon_loss > taus[["q0.999"]]
    )
  is_att <- test$label == "Anomaly"
  is_norm <- !is_att

  # Baselines, all thresholded at their own empirical 99.9% held-out quantile
  bfit <- fit_baselines(pp$X_train, a$bottleneck, seed)
  scores_val <- dplyr::bind_cols(!!ae_label := loss_val, score_baselines(bfit, pp$X_val))
  scores_test <- dplyr::bind_cols(!!ae_label := test$recon_loss, score_baselines(bfit, pp$X_test))
  method_metrics <- purrr::map(names(scores_test), \(m) {
    flag <- scores_test[[m]] > quantile(scores_val[[m]], 0.999)
    det <- detection_by_cell(test, flag)
    tibble::tibble(
      method = m,
      roc_auc = auc_rank(scores_test[[m]], is_att),
      pr_auc = pr_auc_vec(scores_test[[m]], is_att),
      false_alarm_rate = mean(flag[is_norm]),
      det_low = with(dplyr::filter(det, tier == "low", attack_type != "Low-and-Slow Exfiltration"), sum(detected) / sum(n)),
      det_medium = with(dplyr::filter(det, tier == "medium"), sum(detected) / sum(n)),
      det_high = with(dplyr::filter(det, tier == "high"), sum(detected) / sum(n)),
      det_low_slow = with(dplyr::filter(det, attack_type == "Low-and-Slow Exfiltration"), sum(detected) / sum(n))
    )
  }) |>
    dplyr::bind_rows()

  det999 <- detection_by_cell(test, test$flag_999)
  campaigns <- test |>
    dplyr::filter(!is.na(campaign_id)) |>
    dplyr::summarise(len = dplyr::n(), hit = any(flag_999), .by = campaign_id)

  summary <- dplyr::bind_rows(
    tibble::tibble(group = "Model", metric = c(
      "Training MSE (best epoch)", "Validation MSE (best epoch)", "Best epoch",
      "ROC AUC (all attacks)", "GPD shape xi", "GPD scale sigma",
      "tau 0.99", "tau 0.999",
      "False alarms at 99.0%", "False alarms at 99.9%",
      "Low-and-slow campaigns detected", "Null campaign alert rate"
    ), value = c(
      model$train_losses[model$best_epoch], model$best_val_loss, model$best_epoch,
      auc_rank(test$recon_loss, is_att), pot$xi, pot$sigma,
      taus[["q0.99"]], taus[["q0.999"]],
      sum(test$flag_990[is_norm]), sum(test$flag_999[is_norm]),
      mean(campaigns$hit), null_campaign_rate(test, test$flag_999, campaigns$len)
    )),
    det999 |>
      dplyr::summarise(value = sum(detected) / sum(n), .by = attack_type) |>
      dplyr::transmute(group = "Detection at 99.9% by attack", metric = as.character(attack_type), value),
    det999 |>
      dplyr::filter(attack_type != "Low-and-Slow Exfiltration") |>
      dplyr::summarise(value = sum(detected) / sum(n), .by = tier) |>
      dplyr::transmute(group = "Detection at 99.9% by intensity", metric = as.character(tier), value),
    method_metrics |>
      dplyr::transmute(group = "ROC AUC by method", metric = method, value = roc_auc),
    method_metrics |>
      dplyr::transmute(group = "Medium-tier detection at 0.1% FPR by method", metric = method, value = det_medium)
  ) |>
    dplyr::mutate(seed = seed, .before = 1)

  if (!keep) return(summary)

  list(
    seed = seed, study = study, pp = pp, model = model,
    loss_train = loss_train, loss_val = loss_val, residuals = recon_test$residuals,
    test = test, pot = pot, pot_train = pot_train, q = q, taus = taus,
    method_metrics = method_metrics, scores_test = scores_test, scores_val = scores_val,
    det999 = det999, campaigns = campaigns, summary = summary
  )
}
