library(dplyr)
library(knitr)
library(kableExtra)

set.seed(123)
iters = matrix(rep(0,600), nrow = 100, ncol = 6)
errors = matrix(rep(0,600), nrow = 100, ncol = 6)
delta = 1e-6
rho = 0.9
V = "Z"
for(i in 1:100){
  y = autocov(generateChain(list(type = "AR", M = 10000, rho = 0.9))$x)
  for(j in 1:6){
    fit   <- SR1_V(y, V, delta =delta , h_power = (j-1))
    errors[i,j] = L2diff_L2Moment_V(y, fit$support, fit$weights, V)
    iters[i,j] = fit$opt$iter
  }
}
errs = errors - errors[,1]
itrs = colMeans(iters)
E_Mean = colMeans(errs)
Sitrs = apply(iters, 2, sd)/10
E_SE = apply( errs,2,sd)/10
p = seq(0,5,1)/3
df = data.frame(itrs, E_Mean, Sitrs, E_SE, p)
latex_df = df %>%
  mutate(
    `Weight Power (p)` = p,
    `Iterations` = sprintf("%.4f (%.4f)", itrs, Sitrs),
    `Excess L2 Error` = ifelse(p == 0,"0 (Baseline)",
                               paste0("\\num{", sprintf("%.2e", E_Mean), "} (\\num{", sprintf("%.2e", E_SE), "})")
    )
  ) %>%
  select(`Weight Power (p)`, `Iterations`, `Excess L2 Error`)

kbl(latex_df, 
    format = "latex", 
    booktabs = TRUE,     
    escape = FALSE, 
    align = c("c", "l", "c"), 
    linesep = "") %>%   
  kable_styling(latex_options = "hold_position") %>%
  collapse_rows(columns = 1, latex_hline = "major", valign = "middle")


