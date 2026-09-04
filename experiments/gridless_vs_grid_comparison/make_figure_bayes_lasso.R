remove(list=ls())
library(dplyr)

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

############ BLASSO chain summary ###########

BLasso_all = results_all %>% filter(setting_id >=5) %>% mutate(chain_id = setting_id-4) %>% select(sim_index,grid,time,val_diff,chain_id)

#BLasso_all %>% filter(chain_id==1) %>% head()
#obj_comparison_setting5_B100 %>% head()

# calculate Mean Time per chain and per grid
chain_grid_stats <- BLasso_all %>%
  filter(!is.na(grid)) %>%
  group_by(chain_id, grid) %>%
  summarise(
    chain_grid_mean_time = mean(time, na.rm=TRUE), 
    .groups = "drop"
  )


plot_data <- BLasso_all %>%
  filter(!is.na(grid)) %>%
  left_join(chain_grid_stats, by = c("chain_id", "grid")) %>%
  mutate(
    Grid_Label = factor(paste0("s=", grid), 
                        levels = c("s=21", "s=201", "s=2001", "s=20001"))
  )

# Gridless Stats PER CHAIN
gridless_stats <- BLasso_all %>%
  filter(is.na(grid)) %>%
  group_by(chain_id) %>%
  summarise(
    mean_gl_time = mean(time, na.rm=TRUE), 
    .groups = "drop"
  )

ggplot(plot_data, aes(x = chain_grid_mean_time, y = val_diff)) +
  facet_wrap(~chain_id, labeller = label_both, scales="free") +
  geom_vline(data = gridless_stats, aes(xintercept = mean_gl_time), 
             color = "firebrick", linetype = "dashed", linewidth = 0.8) +
  geom_boxplot(aes(group = grid, fill = Grid_Label), 
               width = 0.2, outlier.size=1, outlier.shape = 4,alpha = 0.8) +
  stat_summary(fun = median, geom = "line", aes(group = 1), 
               color = "grey30", linetype = "dotted") +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_manual(values = c("white", "grey85", "grey50", "grey20"), 
                    name = "Grid Size") +
  labs(
    x = "Mean Computation Time (s) [Log Scale]",
    y = "Excess L2 Error [Log Scale]"
  ) +
  theme_test() +
  theme(legend.position = "bottom")


ggsave("experiments/gridless_vs_grid_comparison/results/BayesLasso_plots.pdf", width = 10, height = 6)
