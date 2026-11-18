# Sourced by every page: packages, project functions, cached results, styling.

suppressPackageStartupMessages({
  library(tidyverse)
  library(gt)
  library(patchwork)
  library(scales)
  library(ggraph)
})

purrr::walk(list.files("R", pattern = "\\.R$", full.names = TRUE), source)
cfg <- yaml::read_yaml("params.yaml")

if (!file.exists("results/robust.rds")) stop("Results missing: run `Rscript run_all.R` first.")
main <- readRDS("results/main.rds")
analysis <- readRDS("results/analysis.rds")
graph <- readRDS("results/graph.rds")
robust <- readRDS("results/robust.rds")
train_cor <- readRDS("results/train_cor.rds")

theme_set(
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "#52514e", size = 11),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
)

# Categorical palette in fixed order; normal traffic is neutral grey
attack_colors <- c(
  "Normal (benign)"           = "#9a9a96",
  "Data Exfiltration"         = "#2a78d6",
  "Low-and-Slow Exfiltration" = "#eb6834",
  "Credential Stuffing"       = "#1baf7a",
  "Privilege Escalation"      = "#eda100",
  "Ransomware Staging"        = "#e87ba4",
  "Lateral Movement"          = "#008300"
)
tier_colors <- c("Tier 1" = "#2a78d6", "Tier 2" = "#eb6834", "Tier 3" = "#1baf7a")

# ---- Formatting helpers for inline numbers ------------------------------------
f_pct <- \(x, d = 1) scales::percent(x, accuracy = 10^-d)
f_num <- \(x, d = 2) formatC(x, format = "f", digits = d, big.mark = ",")
f_int <- \(x) scales::comma(round(x))

rob <- robust |>
  summarise(
    mean = mean(value), lo = quantile(value, 0.025), hi = quantile(value, 0.975),
    share_zero = mean(value == 0),
    .by = c(group, metric)
  )
rob_get <- function(metric, group = "Model") {
  rob[rob$metric == metric & rob$group == group, ]
}
# "mean [2.5%-97.5%]" across seeds
f_rob <- function(metric, group = "Model", f = f_num) {
  r <- rob_get(metric, group)
  sprintf("%s [%s–%s]", f(r$mean), f(r$lo), f(r$hi))
}

n_seeds <- length(unique(robust$seed))
d_features <- length(feature_cols)
n_normal_test <- sum(main$test$label == "Normal")
n_attack_test <- sum(main$test$label == "Anomaly")
tau <- main$taus
