remove(list=ls())
library(JuliaCall)
library(foreach)
julia_setup()
julia_command('using Pkg; Pkg.add(["Clarabel", "SumOfSquares", "DynamicPolynomials", "RCall", "Random", "Distributions", "LinearAlgebra", "LinearSolve", "StatsBase"])')
julia_source("codes/SR1_gridless.jl")
julia_source("experiments/ldp_approximation_accuracy/ldp_approximation.jl")

eval_poly = function(cf,xs){
  p_xs = sapply(xs, function(x) sum(cf * x^{0:(length(cf)-1)}))
  return(p_xs)
}


#cf = cf_random[1:2]
#cf_list = cf_random[1:2]
run_ldp_comp <- function(cf, tol){
  
    
    #cf = cf_list[[i]]
    
    # obtain output
    out_r = list()
    pb = progress::progress_bar$new(total = length(tol))
    for(b in 1:length(tol)){
      pb$tick()
      out_r[[b]] = julia_call("ldp_test", cf, tol[b], F)
    }
    
    # summarize output
    d1 = sapply(1:length(tol), function(b){
      c(n = length(cf)-1,
        diff = max(out_r[[b]]$diff),
        normalized_diff = max(out_r[[b]]$normalized_diff),
        maxdegree= max(out_r[[b]]$maxdegree))
    })
    
    data.frame(cbind(tol, t(d1)))
}

  



d = 10000
seed0 = 1234
tol = c(1e-4, 1e-6, 1e-8, 1e-10, 1e-12)
# xs = seq(-1,1,length.out=10000)
# plot(xs,eval_poly(cf_random,xs),type="l")

cf_random = julia_call("coeffs_random_poly",d,seed0)
res1 = run_ldp_comp(cf_random, tol)

cf_truncated_rational = julia_call("coeffs_truncated_rational", d, 0.95)
res2 = run_ldp_comp(cf_truncated_rational, tol)
d_random = res1
d_rational = res2

save(d_random,d_rational,file = "experiments/ldp_approximation_accuracy/results/ldp_approximation_results.Rdata")


