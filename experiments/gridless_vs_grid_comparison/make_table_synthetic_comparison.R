remove(list=ls())
library(dplyr)
library(knitr)
library(kableExtra)

filenames<- list.files(path="experiments/gridless_vs_grid_comparison/results/",pattern="obj_comparison_setting")
filenames = filenames[c(1,6,7,8)]
for(i in 1:4){
  load(paste0("experiments/gridless_vs_grid_comparison/results/",filenames[i]))
}

library(dplyr)
library(ggplot2)
results_all = lapply(ls(pattern = "obj_comparison_setting"),function(x){
  out=eval(parse(text=x));
  setting_id <- as.numeric(sub(".*setting([0-9]+)_.*", "\\1", x))
  out$setting_id = setting_id
  return(out)})
results_all = do.call("rbind",results_all)

df = results_all %>% filter(setting_id <=4) %>% 
  group_by(setting_id,grid) %>% 
  summarise(
    # Calculate Stats
    T_Mean = mean(time),
    T_SE   = sd(time) / sqrt(n()),
    E_Mean = mean(equi_diff),
    E_SE   = sd(equi_diff) / sqrt(n()),
    dh_Mean = mean(dh_diff),
    dh_SE = sd(dh_diff)/sqrt(n()),
    R_Mean = mean(rel_diff),
    R_SE = sd(rel_diff)/sqrt(n()),
    P_positive = mean(equi_diff > -1e-10),
    .groups = "drop"
  ) %>%
  mutate(Method = ifelse(is.na(grid), "Gridless", paste0("Grid ($s=", grid, "$)")),
         Sort_Order = ifelse(is.na(grid), 999999, grid) ) %>% 
  arrange(setting_id,Sort_Order) 

latex_df = df %>%
  mutate(
    Setting = setting_id,
    `Time (s)` = sprintf("%.4f (%.4f)", T_Mean, T_SE),
    `Excess L2 Error` = ifelse(Method == "Gridless","0 (Baseline)",
                               paste0("\\num{", sprintf("%.2e", E_Mean), "} (\\num{", sprintf("%.2e", E_SE), "})")),
    `Relative L2 Error` = ifelse(Method == "Gridless","0 (Baseline)",
                               paste0("\\num{", sprintf("%.2e", R_Mean), "} (\\num{", sprintf("%.2e", R_SE), "})")),
    `Suboptimality Bound` = paste0("\\num{", sprintf("%.2e", dh_Mean), "} (\\num{", sprintf("%.2e", dh_SE), "})")
    )%>%
  select(Setting, Method, `Time (s)`, `Excess L2 Error`, `Relative L2 Error`, `Suboptimality Bound`)

kbl(latex_df, 
    format = "latex", 
    booktabs = TRUE,     
    escape = FALSE, 
    align = c("c", "l", "c", "c", "c", "c"), 
    linesep = "") %>%   
  kable_styling(latex_options = "hold_position") %>%
  collapse_rows(columns = 1, latex_hline = "major", valign = "middle")
