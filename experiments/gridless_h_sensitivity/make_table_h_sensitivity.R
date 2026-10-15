remove(list=ls())
library(dplyr)
library(knitr)
library(kableExtra)

load("experiments/gridless_h_sensitivity/results/h_sensitivity_results.Rdata")
all(out$target_reached)
df = 
  out %>% group_by(setting_id, weight_option, h_power) %>%
  summarise(
    m_iter = mean(iterations),
    se_iter = sd(iterations)/sqrt(n()),
    m_time = mean(time),
    se_time = sd(time)/sqrt(n()))

df_formatted <- df %>%
  ungroup()%>%
  transmute(
    Setting = setting_id,
    Weight = ifelse(
      weight_option == "unweighted",
      "$1$",
      paste0("$1-|\\theta|^{", h_power, "}$")
    ),
    Iterations = sprintf("%.1f (%.1f)", m_iter, se_iter),
    `Time (s)` = sprintf("%.4f (%.4f)", m_time, se_time)
  )


latex_table <- kbl(
  df_formatted,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = c("c", "l", "c", "c"),
  linesep=""
) %>%
  kable_styling(latex_options = "hold_position") %>%
  collapse_rows(columns = 1,latex_hline="major",valign = "middle")

cat(as.character(latex_table))
