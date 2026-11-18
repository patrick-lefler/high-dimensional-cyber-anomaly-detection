# Service dependency graph, Markov transition matrix, PageRank, and
# structural single-point-of-failure measures.
# Edge direction: caller (dependent) -> callee (dependency).

service_nodes <- tibble::tribble(
  ~name,               ~layer,
  "Web Portal",        "Entry point",
  "Mobile App",        "Entry point",
  "Partner API",       "Entry point",
  "API Gateway",       "Edge & identity",
  "Identity Provider", "Edge & identity",
  "Auth Service",      "Edge & identity",
  "MFA Service",       "Edge & identity",
  "Customer Service",  "Business service",
  "Payments Service",  "Business service",
  "Orders Service",    "Business service",
  "Reporting Service", "Business service",
  "Notification Svc",  "Business service",
  "Message Queue",     "Platform",
  "Cache",             "Platform",
  "DNS Resolver",      "Platform",
  "SIEM Logging",      "Platform",
  "Key Management",    "Platform",
  "Customer DB",       "Data store",
  "Payments DB",       "Data store",
  "Orders DB",         "Data store",
  "Data Warehouse",    "Data store",
  "File Store",        "Data store",
  "SaaS CRM",          "External SaaS"
)

service_edges <- tibble::tribble(
  ~from,               ~to,
  "Web Portal",        "API Gateway",
  "Web Portal",        "Identity Provider",
  "Web Portal",        "DNS Resolver",
  "Mobile App",        "API Gateway",
  "Mobile App",        "Identity Provider",
  "Partner API",       "API Gateway",
  "API Gateway",       "Auth Service",
  "API Gateway",       "Customer Service",
  "API Gateway",       "Payments Service",
  "API Gateway",       "Orders Service",
  "API Gateway",       "Reporting Service",
  "API Gateway",       "SIEM Logging",
  "API Gateway",       "DNS Resolver",
  "Identity Provider", "MFA Service",
  "Identity Provider", "Key Management",
  "Identity Provider", "SIEM Logging",
  "Auth Service",      "Identity Provider",
  "Auth Service",      "Cache",
  "MFA Service",       "Notification Svc",
  "Customer Service",  "Customer DB",
  "Customer Service",  "Cache",
  "Customer Service",  "SaaS CRM",
  "Customer Service",  "Message Queue",
  "Payments Service",  "Payments DB",
  "Payments Service",  "Key Management",
  "Payments Service",  "Message Queue",
  "Payments Service",  "Auth Service",
  "Orders Service",    "Orders DB",
  "Orders Service",    "Message Queue",
  "Orders Service",    "Payments Service",
  "Orders Service",    "Customer Service",
  "Reporting Service", "Data Warehouse",
  "Reporting Service", "File Store",
  "Message Queue",     "Notification Svc",
  "Message Queue",     "Data Warehouse",
  "Notification Svc",  "DNS Resolver",
  "Data Warehouse",    "File Store",
  "Customer DB",       "Key Management",
  "Payments DB",       "Key Management",
  "Orders DB",         "Key Management",
  "File Store",        "Key Management",
  "SIEM Logging",      "File Store"
)

build_service_graph <- function(seed) {
  set.seed(seed)
  edges <- dplyr::mutate(service_edges, calls = round(rlnorm(dplyr::n(), log(1000), 1)))
  igraph::graph_from_data_frame(edges, directed = TRUE, vertices = service_nodes)
}

# Row-stochastic P: P[i, j] = Pr(next call goes to j | current service i).
# Dangling services (no outgoing calls) redistribute uniformly.
transition_matrix <- function(g) {
  W <- as.matrix(igraph::as_adjacency_matrix(g, attr = "calls", sparse = FALSE))
  out <- rowSums(W)
  n <- nrow(W)
  P <- W / ifelse(out > 0, out, 1)
  P[out == 0, ] <- 1 / n
  P
}

pagerank_power <- function(P, alpha = 0.85, tol = 1e-12, max_iter = 1000) {
  n <- nrow(P)
  pr <- rep(1 / n, n)
  for (iter in seq_len(max_iter)) {
    new <- alpha * as.numeric(pr %*% P) + (1 - alpha) / n
    err <- sum(abs(new - pr))
    pr <- new
    if (err < tol) break
  }
  list(pagerank = rlang::set_names(pr, rownames(P)), iterations = iter, l1_change = err)
}

# Share of the other services that lose at least one dependency they could
# reach before (including the removed service itself) when a node is removed.
removal_impact <- function(g) {
  nm <- igraph::V(g)$name
  reach <- purrr::map(nm, \(u) setdiff(names(igraph::subcomponent(g, u, mode = "out")), u)) |>
    rlang::set_names(nm)
  purrr::map_dbl(nm, \(v) {
    g2 <- igraph::delete_vertices(g, v)
    others <- setdiff(nm, v)
    lost <- purrr::map_lgl(others, \(u) {
      before <- reach[[u]]
      after <- setdiff(names(igraph::subcomponent(g2, u, mode = "out")), u)
      length(setdiff(before, after)) > 0
    })
    mean(lost)
  }) |>
    rlang::set_names(nm)
}

analyse_service_graph <- function(g, gcfg) {
  P <- transition_matrix(g)
  stopifnot(all(abs(rowSums(P) - 1) < 1e-12))

  pr <- pagerank_power(P, gcfg$damping, gcfg$tol, gcfg$max_iter)
  pr_igraph <- igraph::page_rank(g, damping = gcfg$damping, weights = igraph::E(g)$calls)$vector

  sens <- purrr::map(unlist(gcfg$damping_sensitivity), \(a) {
    tibble::tibble(alpha = a, name = rownames(P), pagerank = pagerank_power(P, a, gcfg$tol, gcfg$max_iter)$pagerank)
  }) |>
    dplyr::bind_rows()
  base <- dplyr::filter(sens, alpha == gcfg$damping)
  kendall <- sens |>
    dplyr::summarise(
      kendall_tau = cor(pagerank, base$pagerank[match(name, base$name)], method = "kendall"),
      .by = alpha
    )

  art <- names(igraph::articulation_points(igraph::as_undirected(g)))
  impact <- removal_impact(g)

  table <- tibble::tibble(
    name = rownames(P),
    layer = igraph::V(g)$layer,
    pagerank = pr$pagerank,
    betweenness = igraph::betweenness(g, directed = TRUE, weights = NA),
    articulation_point = name %in% art,
    removal_impact = impact[name]
  ) |>
    dplyr::mutate(
      criticality_tier = dplyr::case_when(
        articulation_point | removal_impact >= gcfg$tier1_removal_impact ~ "Tier 1",
        pagerank >= median(pagerank) ~ "Tier 2",
        .default = "Tier 3"
      )
    ) |>
    dplyr::arrange(dplyr::desc(pagerank))

  list(
    graph = g, P = P, pagerank = pr,
    max_diff_vs_igraph = max(abs(pr$pagerank - pr_igraph[names(pr$pagerank)])),
    sensitivity = sens, kendall = kendall, table = table,
    spearman_pr_impact = cor(table$pagerank, table$removal_impact, method = "spearman"),
    spearman_btw_impact = cor(table$betweenness, table$removal_impact, method = "spearman")
  )
}
