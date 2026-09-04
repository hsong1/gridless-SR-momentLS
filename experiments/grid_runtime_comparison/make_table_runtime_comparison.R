remove(list=ls())
library(dplyr)
library(knitr)

load("experiments/grid_runtime_comparison/results/grid_time.Rdata")

df_summary= comp1 %>% data.frame() %>% 
  group_by(method) %>% 
  mutate(time_ms = time_s*10^3) %>% 
  select(method,time_ms,value) %>% 
  summarise(
    mean_time_ms = mean(time_ms),
    se_time_ms = sd(time_ms) / sqrt(n()),
    mean_value = mean(value),
    se_value = sd(value) / sqrt(n())
  )

# latex tble
df_summary%>%
  mutate(
    Method = ifelse(method=="method1","Low-degree approximation", "Exact polynomial"),
    `Time (ms)` = sprintf("%.3f (%.3f)", mean_time_ms, se_time_ms),
    Value = sprintf("%.3f (%.3f)", mean_value, se_value)
  ) %>% select(Method, `Time (ms)`, Value) %>%
  kable(
    format = "latex",
    booktabs = F,  
    linesep = "",
    escape = FALSE,   
    align = "lcc"
  ) 


