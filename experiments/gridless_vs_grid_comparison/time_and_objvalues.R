remove(list=ls())
library(foreach)
library(doParallel)
library(momentLS)
############### simulation setting #####################
grid_sizes = c(20,200,2000,20000)+1
# grid_sizes=c(11,21,31,41)
# Read arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
sim_setting <- if (length(args) > 0) as.numeric(args[1]) else 13
B <- if (length(args) >= 2) as.numeric(args[2]) else 5

# set up simulation
filename = paste0("obj_comparison_setting",sim_setting,"_B",B)
filepath = paste0("experiments/gridless_vs_grid_comparison/results/",filename,".Rdata")

source("codes/SR1_V.R")
source("experiments/gridless_vs_grid_comparison/simulate_examples.R")
#############################################

#b=2
#pb <- progress::progress_bar$new(total=B)
results = foreach(b = 1:B,.combine="rbind")%do%{
  #pb$tick()
  cat("b:",b,"\n")
  example_b  = simulate_examples(sim_setting,seed0 = 1367+b)
  
  nsum = if(example_b$weight_option == "unweighted") sqrt(sum(example_b$y^2)) else 4*sum(abs(example_b$y))
  # fit gridless verison
  start_time <- proc.time()
  m_gl <- with(example_b, SR1_gridless_R(y, delta = delta, weight_option = weight_option, V = V, maxit = 1000, tol=1e-10))
  # m_gl <- with(example_b, SR1_gridless_R(y, delta = delta, weight_option = weight_option, V = V, maxit = 500, tol=1e-10, verbose =T))
  end_time <- proc.time()
  t_gl <- (end_time - start_time)["elapsed"]
  
  lb_gl <- abs(m_gl$min_grad) *nsum
  
  obj_gl <- with(m_gl,L2diff_L2Moment_V(y = example_b$y, support, weights, example_b$V) )
  
  # fit grid SR
  m_equi = list()
  m_log = list()
  grid_times <- c()
  log_times = c()
  dh_grid = c()
  for(g in 1:length(grid_sizes)){
    start_time <- proc.time()
    alphaGrid = momentLS::makeGrid(upper_threshold = 1-example_b$delta_grid,
                                   nX = grid_sizes[g],cm = ifelse(example_b$V=="Z",F,T),scale = "equidist")
    m_equi[[g]] =
      tryCatch(
        {with(example_b, SR1_V(y = y,V = V,
                               alphaGrid = alphaGrid,tol = 1e-6))},
      error = function(e){return(e$message)}
    )
    end_time <- proc.time()
    grid_times[g] <- (end_time - start_time)["elapsed"]
    start_time = proc.time()
    alphaGrid = momentLS::makeGrid(upper_threshold = 1-example_b$delta_grid,
                                   nX = grid_sizes[g],cm = ifelse(example_b$V=="Z",F,T),scale = "log")
    m_log[[g]] =
      tryCatch(
        {with(example_b, SR1_V(y = y,V = V,
                               alphaGrid = alphaGrid,tol = 1e-6))},
        error = function(e){return(e$message)}
      )
    end_time = proc.time()
    log_times[[g]] = (end_time-start_time)["elapsed"]
  }
  
  obj_equi <- sapply(m_equi, function(m){
    if(is.character(m)){return(NA)}else{
      return(with(m, L2diff_L2Moment_V(y = example_b$y,support = support,weights = weights,V = example_b$V)))
    }
  })
  obj_log <- sapply(m_log, function(m){
    if(is.character(m)){return(NA)}else{
      return(with(m, L2diff_L2Moment_V(y = example_b$y,support = support,weights = weights,V = example_b$V)))
    }
  })
  dh_grid <- sapply(m_equi, function(m){
    if(is.character(m)){return(NA)}else{
      return(with(m, SR1_gridless_R(example_b$y, delta = example_b$delta, weight_option = example_b$
                                      weight_option, V = example_b$V, maxit = 1,init_supp = support, init_weight = weights, tol=1)$min_grad))
    }
  })
  
  (equi_diff = (obj_equi - obj_gl))
  (log_diff = obj_log - obj_gl)
  (rel_diff = (obj_equi - obj_gl)/obj_gl)
  (lb_g = abs(dh_grid)* nsum)
  
  data.frame(sim_index= b, grid = c(NA,grid_sizes), time = c(t_gl,grid_times), equi_diff = c(NA,equi_diff), log_diff = c(NA, log_diff), rel_diff = c(NA, rel_diff), dh_diff = c(lb_gl,lb_g))
}

assign(filename,results)
save(list=filename,file = filepath)


