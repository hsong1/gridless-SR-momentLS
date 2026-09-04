remove(list=ls())
library(dplyr)
library(knitr)
library(kableExtra)

filenames<- list.files(path="experiments/gridless_vs_grid_comparison/results/",pattern="obj_comparison_setting")
for(i in 1:length(filenames)){
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
    E_Mean = mean(val_diff),
    E_SE   = sd(val_diff) / sqrt(n()),
    P_positive = mean(val_diff > -1e-10),
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
                               paste0("\\num{", sprintf("%.2e", E_Mean), "} (\\num{", sprintf("%.2e", E_SE), "})")
    )
  ) %>%
  select(Setting, Method, `Time (s)`, `Excess L2 Error`)

kbl(latex_df, 
    format = "latex", 
    booktabs = TRUE,     
    escape = FALSE, 
    align = c("c", "l", "c", "c"), 
    linesep = "") %>%   
  kable_styling(latex_options = "hold_position") %>%
  collapse_rows(columns = 1, latex_hline = "major", valign = "middle")
