remove(list=ls())
library(foreach)
library(doParallel)
library(momentLS)

source("codes/SR1_V.R")
source("experiments/gridless_vs_grid_comparison/simulate_examples.R")


### simulation setting ###
h_powers=c(1,2,5,10)
settings=1:4
B=100
###########################

# setting_id = 3; b = 1; j=1; 
out<- foreach(setting_id = settings,.combine="rbind")%do%{
  cat("setting_id:",setting_id,"\n")
  
  # h functions to test
  wchoices = data.frame(weight_option="weighted", h_powers)
  if(setting_id >=3){
    wchoices = rbind(data.frame(weight_option="unweighted",h_powers=1),wchoices)
  }

  foreach(b = 1:B,.combine = "rbind")%do%{
    
    cat("b:",b,"\n")
    # simulate example
    example_b <- simulate_examples(setting_id, seed0 = 1367 + b)

    # fit gridless_ref 
      m_gl_ref <- SR1_gridless_R(y = example_b$y, delta = example_b$delta, 
        weight_option = wchoices[1,]$weight_option, V = example_b$V, 
        maxit = 10000, tol=1e-12, h_power = wchoices[1,]$h_powers)
      obj_ref <- L2diff_L2Moment_V(y = example_b$y, m_gl_ref$support, m_gl_ref$weights, example_b$V)

    # for each weight function h
    foreach(j = 1:nrow(wchoices),.combine = "rbind")%do%{

      # fit gridless for current configuration
      wchoice = wchoices[j,]
      m_gl <- SR1_gridless_R(y = example_b$y, delta = example_b$delta, 
        weight_option = wchoice$weight_option, V = example_b$V, 
        maxit = 10000, tol=1e-10, h_power = wchoice$h_powers,save_trajectory = TRUE)
      
      obj_all <- sapply(m_gl$trajectory, function(x) L2diff_L2Moment_V(example_b$y,support=x$support,weights=x$weights,V=example_b$V))

      rel_target <- 1e-6
      rel_error <- (obj_all - obj_ref)/obj_ref

      maxit_j <- which.max(rel_error <= rel_target) # find the first TRUE

      # timed
      start_time <- proc.time()
      m_gl <- SR1_gridless_R(y = example_b$y, delta = example_b$delta, 
        weight_option = wchoice$weight_option, V = example_b$V, 
        maxit = maxit_j, tol=1e-10, h_power = wchoice$h_powers)
      end_time <- proc.time()

      # summary
      t_gl <- (end_time - start_time)["elapsed"]
      obj_gl <- with(m_gl,L2diff_L2Moment_V(y = example_b$y, support, weights, example_b$V) )
      rel_val_diff <- (obj_gl - obj_ref)/obj_ref



      # # suboptimality bounds
      # norm1_y <- with(example_b, {
      #   if(V=="Z") abs(y[1])+2*sum(abs(y[-1])) else sum(abs(y))
      # })
      # ldp_tol <- 1e-12 / (2*norm1_y) # delta_h(mu)<= delta_h_ldp + 1e-12
      # bound_gl <- subopt_bound(example_b$y,example_b$delta,
      #   weight_option = wchoice$weight_option,
      #   V = example_b$V,
      #   supp=m_gl$support,weight = m_gl$weights,ldp_tol=ldp_tol, h_power = wchoice$h_powers)

      # output
      data.frame(
        setting_id = setting_id,
        sim_index= b, 
        weight_option = wchoice$weight_option,
        h_power = wchoice$h_powers,
        time = unname(t_gl), 
        rel_target=rel_target,
        val = obj_gl,
        val_ref = obj_ref,
        rel_val_diff = rel_val_diff,
        iterations = m_gl$iterations,
        target_reached = rel_val_diff<=rel_target
      )

    } # end of j
  } #end of b
} # end of setting



save(out, file = "experiments/gridless_h_sensitivity/results/h_sensitivity_results.Rdata")
