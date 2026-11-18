# Cyber Anomaly Detection & Service Dependency Risk

Markov transition matrices, PageRank, a shallow autoencoder, and Extreme Value Theory thresholds: a simulation study in R.

**Author:** Patrick Lefler
**Published:** 2026-08-12 · **Revised:** 2026-11-18
**Site:** https://patrick-lefler.com/high-dimensional-cyber-anomaly-detection
**Repository:** https://github.com/patrick-lefler/high-dimensional-cyber-anomaly-detection

## What this project does

Signature rules and fixed thresholds miss attacks that stay inside normal ranges. This project tests a two-layer alternative on simulated enterprise telemetry:

- **Detection layer.** A shallow autoencoder, written in base R with hand-derived gradients and the Adam optimizer, learns what normal host activity looks like across 19 telemetry features. It scores each 5-minute window by how badly it can reconstruct it. Alert thresholds come from Extreme Value Theory: a Generalized Pareto Distribution fitted to the tail of reconstruction losses on held-out hosts.
- **Structural layer.** A simulated graph of 23 services becomes a Markov transition matrix of call flows. The services are ranked by PageRank, articulation points and a node-removal simulation. The ranking sets criticality tiers that decide which alerts may trigger automated containment and which need a person.

The project is built to be honest about what the detector can and cannot do. It compares the autoencoder against simple baselines, repeats the full pipeline over 30 random seeds, and reports uncertainty on every threshold.

## Key findings

All results come from simulated data with injected attacks. They show how the methods behave under controlled conditions, not how they would perform on a real network. Figures are means across 30 seeds, with the 95% range in brackets.

- **The autoencoder adds nothing over simpler statistics on this data.** Its ROC AUC is 0.794 [0.734–0.824]. Mahalanobis distance scores 0.817 [0.765–0.843]. The network learns roughly what principal component analysis learns.
- **It catches loud attacks and misses quiet ones.** At the 99.9% threshold it detects 77% of high-intensity, 25% of medium-intensity and under 1% of low-intensity attack windows. It misses data exfiltration at every intensity, and catches low-and-slow exfiltration campaigns no more often than chance.
- **False alarms run above the nominal rate and do not transfer cleanly to new hosts.** At the 99.9% threshold the detector raised 14 [1–52] false alarms per 8,064 normal windows, against about 8 expected.
- **The loss tail is heavy and hard to model.** Rare benign count spikes (new service creations, privilege escalations) cause 96% of the most extreme normal losses. The diagnostics show no stable tail cutoff for a single Generalized Pareto fit.
- **Base rates dominate.** At one attack per 100,000 windows and a million windows a day, the 99.9% threshold would raise about 2,000 false alarms a day, and fewer than 1% of alerts would be real.
- **PageRank is a reasonable first screen for dependency risk, not a single-point-of-failure test.** It tracks node-removal impact (Spearman correlation 0.73), but it ranks one of the two articulation points (Customer Service) 16th of 23 services.

**Recommendations:**

- Route containment by service criticality, with human approval for the most critical services.
- Report and validate the features behind each alert.
- Monitor drift by effect size (Wasserstein distance) and live alert rate, rather than weekly significance tests.

## Site structure

| Page | Contents |
|---|---|
| Overview (`index.qmd`) | Introduction, headline results, recommendations, limitations, reproducibility |
| Dependency Graph (`dependency-graph.qmd`) | Transition matrix, PageRank, articulation points, removal impact, criticality tiers |
| Detector (`detector.qmd`) | Telemetry simulator, attack scenarios, autoencoder, gradient check, architecture sensitivity |
| Thresholds (`thresholds.qmd`) | GPD fit, block-bootstrap intervals, diagnostics, calibration, declustering |
| Evaluation (`evaluation.qmd`) | Scorecard, baselines, 30-seed robustness, base rates, attribution, containment routing, drift |
| Reproducibility (`reproducibility.qmd`) | Quick start, pipeline, result fingerprints, seeds, test results, caching, licenses, software environment |

## Reproducing the results

```bash
git clone https://github.com/patrick-lefler/high-dimensional-cyber-anomaly-detection.git
cd high-dimensional-cyber-anomaly-detection
Rscript -e 'renv::restore()'                         # install pinned package versions
Rscript run_all.R                                    # regenerate results/*.rds (under a minute)
Rscript -e 'testthat::test_dir("tests/testthat")'    # unit tests
rm -rf _freeze && quarto render                      # render the site into docs/
```

`results/` is not committed, so run `run_all.R` before the first render. Every simulation parameter and seed is in `params.yaml`, labelled as an assumption. The site's Reproducibility page lists fingerprints of the key results, so you can confirm that a rerun matches exactly.

## Repository structure

```
├── renv.lock            # pinned package versions (renv::restore())
├── params.yaml          # all simulation, model, EVT and policy parameters
├── run_all.R            # one-command pipeline: writes results/*.rds
├── R/
│   ├── sim.R            # telemetry simulator and attack injection
│   ├── ae.R             # autoencoder, Adam, early stopping, gradient check
│   ├── evt.R            # GPD thresholds (incl. xi = 0), block bootstrap
│   ├── baselines.R      # z-score, PCA, Mahalanobis, Isolation Forest, AUC
│   ├── graph.R          # service graph, transition matrix, PageRank, removal impact
│   ├── drift.R          # PSI, Wasserstein distance, KS with Holm correction
│   ├── pipeline.R       # run_pipeline(seed): simulate -> train -> threshold -> evaluate
│   └── pages/setup.R    # shared setup sourced by every page
├── tests/testthat/      # unit tests
├── *.qmd                # site pages
├── _quarto.yml, theme.scss, _brand.yml
└── docs/                # rendered site
```

## Tech stack

R and Quarto. Packages: tidyverse, recipes, yardstick, evd, igraph, ggraph, isotree, binom, gt, patchwork, scales, furrr, testthat and sessioninfo. The autoencoder uses no deep-learning framework.

## License

MIT License. See [LICENSE](LICENSE).

## Contact

Patrick Lefler · [LinkedIn](https://www.linkedin.com/in/patricklefler/) · [patrick-lefler.github.io](https://patrick-lefler.github.io) · [Substack](https://substack.com/@pflefler)
