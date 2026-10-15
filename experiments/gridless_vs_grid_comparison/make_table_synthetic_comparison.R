remove(list=ls())
library(dplyr)
library(knitr)
library(kableExtra)

filenames<- list.files(path="experiments/gridless_vs_grid_comparison/results/",pattern="obj_comparison_setting")
for(i in 1:length(filenames)){
  load(paste0("experiments/gridless_vs_grid_comparison/results/",filenames[i]))
}

results_all = lapply(ls(pattern = "obj_comparison_setting"),function(x){
  out=eval(parse(text=x));
  setting_id <- as.numeric(sub(".*setting([0-9]+)_.*", "\\1", x))
  out$setting_id = setting_id
  return(out)})
results_all = do.call("rbind",results_all)

all(results_all$conv) # check convergence

# scale_choice = "equidist"
make_comparison_table <- function(results_all, scale_choice) {
  summary_df <- results_all %>%
    filter(setting_id <= 4, grid_scale == scale_choice) %>%
    group_by(setting_id, grid) %>%
    summarise(
      T_Mean = mean(time),
      T_SE = sd(time) / sqrt(n()),
      E_Mean = mean(val_diff),
      E_SE = sd(val_diff) / sqrt(n()),
      R_Mean = mean(rel_val_diff),
      R_SE = sd(rel_val_diff) / sqrt(n()),
      SubOpt_Mean = mean(subopt_bounds),
      SubOpt_SE = sd(subopt_bounds) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      Method = ifelse(is.na(grid), "Gridless", paste0("Grid ($s=", grid, "$)")),
      Sort_Order = ifelse(is.na(grid), Inf, grid)
    ) %>%
    arrange(setting_id, Sort_Order)

  scientific_cell <- function(avg, se) {
    sprintf("\\num{%.2e} (\\num{%.2e})", avg, se)
  }

  latex_df <- summary_df %>%
    transmute(
      Setting = setting_id,
      Method,
      `Time (s)` = sprintf("%.4f (%.4f)", T_Mean, T_SE),
      `Excess L2 Error` = ifelse(
        is.na(grid), "0 (Baseline)", scientific_cell(E_Mean, E_SE)
      ),
      `Relative Excess L2 Error` = ifelse(
        is.na(grid), "0 (Baseline)", scientific_cell(R_Mean, R_SE)
      ),
      `Suboptimality Bound` = scientific_cell(SubOpt_Mean, SubOpt_SE)
    )

  latex_code <- kbl(
    latex_df,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = c("c", "l", "c", "c", "c", "c"),
    linesep = ""
  ) %>%
    kable_styling(latex_options = "hold_position") %>%
    collapse_rows(columns = 1, latex_hline = "major", valign = "middle")
  return(list(summary_df=summary_df, latex_code = latex_code))
}

out_equidist <- make_comparison_table(results_all, "equidist")
out_log <- make_comparison_table(results_all, "log")

cat(as.character(out_equidist$latex_code))
cat(as.character(out_log$latex_code))
