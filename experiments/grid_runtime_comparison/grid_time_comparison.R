library(microbenchmark)
remove(list=ls())
# polyApprox = T
# one execution of SR algorithm 
SR1_one <- function(x,n_alphas,polyApprox=T){
  
  dhat = tune_delta(x,nSplits = 5)$delta*0.8
  r = autocov(x)
  alphaGrid = makeGrid(nX = n_alphas, upper_threshold = 1 - dhat, cm = FALSE, scale = "equidist")
  
  ### timed part ###
  start = proc.time()
  s_alpha = sqrt((1 + alphaGrid^2)/(1 - alphaGrid^2))
  XtX =  makeXtX(alphaGrid, s_x = s_alpha)
  if(polyApprox){
    Xtr = Xtr_cpp(x = alphaGrid, a = r,epsilon = 1e-10 * length(x)) # function ldp tolerance = epsilon / n
  }else{
    Xtr = computeXtr(x = alphaGrid,r = r,s_x = s_alpha )
  }
  m = SR1(r,alphaGrid = alphaGrid,precomputed=list(s_alpha=s_alpha,XtX=XtX,Xtr=Xtr,alphaGrid=alphaGrid,input=r))  
  end = proc.time()
  ###
  
  tm = end["elapsed"] - start["elapsed"]
  return(list(tm=tm,m = m))
}



# i=1
# M=10000
# n_alphas = 1001

compare_runtime<- function(M, n_alphas, runs=100){
  
  results <- data.frame(method = character(), time_s = numeric(), value = numeric(),  run_id = integer())
  pb=progress::progress_bar$new(total=runs)
  
  set.seed(1)
  rho = runif(1) # fix one rho
  
  for (i in 1:runs) {
    pb$tick()
    set.seed(12345+i)
    
    # Simulate an AR(1) chain
    chainParams = list(type="AR",  M = M, rho = rho)
    ch = generateChain(chainParams)
    out.pa = SR1_one(ch$x,n_alphas = n_alphas,polyApprox = T)
    out.exact = SR1_one(ch$x,n_alphas = n_alphas,polyApprox = F)
    
    r= autocov(ch$x)
    
    # function values
    val1<-momentLS::L2diff_L2Moment(r = r,support = out.pa$m$support,weights = out.pa$m$weights)
    val2<-momentLS::L2diff_L2Moment(r = r,support = out.exact$m$support,weights = out.exact$m$weights)
    
    results <- rbind(results,
                     data.frame(method = "method1", time_s = out.pa$tm,    value = val1, run_id=i),
                     data.frame(method = "method2", time_s = out.exact$tm, value = val2, run_id=i)
    )
  }
  return(results)
}

comp1 = compare_runtime(M = 1e5,n_alphas = 1001, runs = 100)
save(comp1, file="experiments/grid_runtime_comparison/results/grid_time.Rdata")

