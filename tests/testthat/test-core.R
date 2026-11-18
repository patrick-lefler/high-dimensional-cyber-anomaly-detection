root <- testthat::test_path("..", "..")
purrr::walk(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE), source)
cfg <- yaml::read_yaml(file.path(root, "params.yaml"))

test_that("gpd_tau handles xi = 0 with the exponential limit", {
  tau0 <- gpd_tau(u = 1, sigma = 0.3, xi = 0, N = 16000, Nu = 800, p = 0.001)
  expect_true(is.finite(tau0))
  expect_equal(tau0, 1 + 0.3 * log(800 / (16000 * 0.001)))
  expect_equal(gpd_tau(1, 0.3, 1e-4, 16000, 800, 0.001), tau0, tolerance = 1e-3)
})

test_that("gpd_tau returns u at the threshold's own exceedance probability", {
  expect_equal(gpd_tau(2, 0.5, 0.2, 1000, 50, 50 / 1000), 2)
})

test_that("analytic gradients match finite differences", {
  set.seed(1)
  X <- matrix(rnorm(32 * 6), 32, 6)
  expect_true(all(gradient_check(ae_init(c(6, 3, 6)), X)$rel_error < 1e-6))
  expect_true(all(gradient_check(ae_init(c(6, 4, 2, 4, 6)), X)$rel_error < 1e-6))
})

test_that("transition matrix is row-stochastic and PageRank matches igraph", {
  res <- analyse_service_graph(build_service_graph(cfg$graph$seed), cfg$graph)
  expect_true(all(abs(rowSums(res$P) - 1) < 1e-12))
  expect_equal(sum(res$pagerank$pagerank), 1, tolerance = 1e-10)
  expect_lt(res$max_diff_vs_igraph, 1e-8)
})

test_that("removal impact of an entry point is zero and of a shared sink is positive", {
  imp <- removal_impact(build_service_graph(cfg$graph$seed))
  expect_equal(unname(imp["Web Portal"]), 0)
  expect_gt(unname(imp["Key Management"]), 0.5)
})

test_that("drift metrics are near zero for identical samples", {
  set.seed(2)
  x <- rnorm(5000)
  expect_lt(psi(x, x), 1e-10)
  expect_lt(wasserstein1(x, x), 1e-10)
})

test_that("rank AUC is 1 for perfect separation and 0.5 for ties", {
  expect_equal(auc_rank(c(1, 2, 3, 4), c(FALSE, FALSE, TRUE, TRUE)), 1)
  expect_equal(auc_rank(rep(1, 4), c(FALSE, FALSE, TRUE, TRUE)), 0.5)
})

test_that("simulated host has one row per 5-minute window and all features", {
  set.seed(3)
  h <- simulate_host(1, cfg$simulation, n_days = 2)
  expect_equal(nrow(h), 2 * 288)
  expect_true(all(feature_cols %in% names(h)))
  expect_false(anyNA(h[feature_cols]))
})
