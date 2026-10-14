remove(list=ls())
library(foreach)
library(doParallel)
library(momentLS)
############### simulation setting #####################
grid_sizes = c(20,200,2000,20000)+1
# grid_sizes=c(11,21,31,41)
# Read arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
sim_setting <- if (length(args) > 0) as.numeric(args[1]) else 3
B <- if (length(args) >= 2) as.numeric(args[2]) else 5
grid_scale <- if (length(args) >= 3) args[3] else "equidist"

# set up simulation
filename = paste0("obj_comparison_setting",sim_setting,"_B",B,"_scale_",grid_scale)
filepath = paste0("experiments/gridless_vs_grid_comparison/results/",filename,".Rdata")

source("codes/SR1_V.R")
source("experiments/gridless_vs_grid_comparison/simulate_examples.R")
#############################################

#b=37
#pb <- progress::progress_bar$new(total=B) 
results = foreach(b = 1:B,.combine="rbind")%do%{
  #pb$tick()
  cat("b:",b,"\n")
  example_b  = simulate_examples(sim_setting,seed0 = 1367+b)
  
  if(b==1){ # to remove first-call compiliation overhead
    invisible(with(example_b, SR1_gridless_R(y, delta = delta, weight_option = weight_option, V = V, maxit = 10000, tol=1e-10)))
  }

  # fit gridless verison
  start_time <- proc.time()
  m_gl <- with(example_b, SR1_gridless_R(y, delta = delta, weight_option = weight_option, V = V, maxit = 10000, tol=1e-10))
  # m_gl <- with(example_b, SR1_gridless_R(y, delta = delta, weight_option = weight_option, V = V, maxit = 2000, tol=1e-10, verbose =T))
  end_time <- proc.time()
  t_gl <- (end_time - start_time)["elapsed"]
  
  obj_gl <- with(m_gl,L2diff_L2Moment_V(y = example_b$y, support, weights, example_b$V) )
  
  
  # fit grid SR
  m_grid <-list()
  grid_times <- c()
  for(g in 1:length(grid_sizes)){
    start_time <- proc.time()
    alphaGrid = momentLS::makeGrid(upper_threshold = 1-example_b$delta_grid,
                                   nX = grid_sizes[g],cm = ifelse(example_b$V=="Z",F,T),scale = grid_scale)
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

  # suboptimality bounds
  norm1_y <- with(example_b, {
    if(V=="Z") abs(y[1])+2*sum(abs(y[-1])) else sum(abs(y))
  })
  ldp_tol <- 1e-12 / (2*norm1_y) # delta_h(mu)<= delta_h_ldp + 1e-12
  bound_gl <- with(example_b, subopt_bound(y,delta,weight_option,V,supp=m_gl$support,weight = m_gl$weights,ldp_tol=ldp_tol))
  bound_grid <- sapply(m_grid, function(m){
    if(is.character(m)){return(NA)}else{
      return(with(example_b, subopt_bound(y,delta,weight_option,V,supp=m$support,weight = m$weights,ldp_tol=ldp_tol)))
    }
  })
  
  (val_diff = obj_grid - obj_gl)
  (rel_val_diff = val_diff / obj_gl)

  # convergence
  conv_gl = m_gl$converged
  conv_grid = sapply(m_grid, function(m){
    if(is.character(m)){return(NA)}else{return(m$opt$convergence)}
  })

  data.frame(
    sim_index= b, 
    grid_scale = grid_scale,
    grid = c(NA,grid_sizes), 
    time = c(t_gl,grid_times), 
    val = c(obj_gl,obj_grid), 
    val_diff = c(0,val_diff),
    rel_val_diff = c(0, rel_val_diff),
    conv = c(conv_gl, conv_grid),
    subopt_bounds = c(bound_gl, bound_grid)
  )
}

assign(filename,results)
save(list=filename,file = filepath)

