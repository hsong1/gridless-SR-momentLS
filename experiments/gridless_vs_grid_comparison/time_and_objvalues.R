remove(list=ls())
library(foreach)
library(doParallel)
library(momentLS)
############### simulation setting #####################
grid_sizes = c(20,200,2000,20000)+1
# grid_sizes=c(11,21,31,41)
# Read arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
sim_setting <- if (length(args) > 0) as.numeric(args[1]) else 2
B <- if (length(args) >= 2) as.numeric(args[2]) else 10

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
  
  # fit gridless verison
  start_time <- proc.time()
  m_gl <- with(example_b, SR1_gridless_R(y, delta = delta, weight_option = weight_option, V = V, maxit = 1000, tol=1e-10))
  # m_gl <- with(example_b, SR1_gridless_R(y, delta = delta, weight_option = weight_option, V = V, maxit = 500, tol=1e-10, verbose =T))
  end_time <- proc.time()
  t_gl <- (end_time - start_time)["elapsed"]
  
  obj_gl <- with(m_gl,L2diff_L2Moment_V(y = example_b$y, support, weights, example_b$V) )
  
  # fit grid SR
  m_grid <-list()
  grid_times <- c()
  for(g in 1:length(grid_sizes)){
    start_time <- proc.time()
    alphaGrid = momentLS::makeGrid(upper_threshold = 1-example_b$delta_grid,
                                   nX = grid_sizes[g],cm = ifelse(example_b$V=="Z",F,T),scale = "equidist")
    m_grid[[g]] = 
      tryCatch(
        {with(example_b, SR1_V(y = y,V = V,
                               alphaGrid = alphaGrid,tol = 1e-6))},
      error = function(e){return(e$message)}
    )
    end_time <- proc.time()
    grid_times[g] <- (end_time - start_time)["elapsed"]
  }
  
  obj_grid <- sapply(m_grid, function(m){
    if(is.character(m)){return(NA)}else{
      return(with(m, L2diff_L2Moment_V(y = example_b$y,support = support,weights = weights,V = example_b$V)))
    }
  })
  
  (val_diff = obj_grid - obj_gl)
  
  data.frame(sim_index= b, grid = c(NA,grid_sizes), time = c(t_gl,grid_times), val_diff = c(NA,val_diff))
}

assign(filename,results)
save(list=filename,file = filepath)

