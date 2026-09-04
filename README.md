# Gridless SR-momentLS

Code for reproducing the numerical results in **A Gridless Support Reduction Algorithm for Nonparametric Mixture Moment Problems**.

## Files

- `MWE.R`: Minimal working example for fitting a gridless momentLS object
- `codes/`: Implementation of the gridless SR-momentLS
- `experiments/`: Code for reproducing the numerical experiments in the paper.
  - `grid_runtime_comparison/`: Runtime comparison in Table 1
  - `ldp_approximation_accuracy/`: Low-degree polynomial approximation in Section 5.1
  - `gridless_vs_grid_comparison/`: Gridless versus grid-based comparisons for the synthetic and Bayesian LASSO examples in Section 5.1 and 5.3.
- `figures/`: Code for generating illustrative figures in the paper.

Scripts beginning with `make_table_` or `make_figure_` generate the tables and figures reported in the paper.