simulate_examples<- function(sim_setting, seed0){
  # fix MC dynamics if sim_setting == 4
  if(sim_setting == 4){
    K = 1000
    set.seed(2468)
    g = matrix(rnorm(K*1),ncol=1)
    discreteMC = simulate_discreteMC(nStates = K)
  }
  
  set.seed(seed0)
  if(sim_setting == 1){
    n=1000
    weights <- c(0.25, 0.50, 0.25)
    thetas  <- c(0.10, 0.30, 0.80)
    components <- sample(1:3, size = n, replace = TRUE, prob = weights)
    x <- rgeom(n, prob = thetas[components])
    y = table(factor(x,levels=0:max(x)))/length(x)
    delta_grid =1e-6
    return(list(y=y,V="N",weight_option="weighted",delta = 0, delta_grid=delta_grid))
  }
  if(sim_setting == 2){
    a=1; b=3; n=1000
    theta <- rbeta(n, shape1 = a, shape2 = b)
    x = rgeom(n, prob = theta)
    #max(x)
    # convert it into pmf
    y = table(factor(x,levels=0:max(x)))/length(x)
    delta_grid =1e-6
    return(list(y=y,V="N",weight_option="weighted",delta = 0, delta_grid=delta_grid))
  }
  if(sim_setting == 3){
    n = 10000
    ch = momentLS::generateAR1Chain(n,0.95)
    dhat = momentLS::tune_delta(ch$x,nSplits = 5)$delta * 0.8
    r = autocov(ch$x)
    return(list(y=r,V="Z",weight_option="unweighted",delta=dhat,delta_grid=dhat))
  }
  if(sim_setting == 4){
    n = 10000
    chainParams = list(type="MH",  M = n,  nStates = K, g = g, discreteMC = discreteMC,d=1)
    ch = generateChain(chainParams)
    dhat = momentLS::tune_delta(ch$x,nSplits = 5)$delta * 0.8
    r = autocov(ch$x)
    return(list(y=r,V="Z",weight_option="unweighted",delta=dhat,delta_grid=dhat))
  }
  if(sim_setting >=5){
    # load data
    load("experiment/BayesLasso/concreteStrength.Rdata")
    # load bayesLasso sampler
    source("experiment/BayesLasso/bayesLasso.R")
    
    # sample
    warmup = 1000; n = 10000
    ch = run.bl(X=X_concrete_scaled,Y=Y_concrete-mean(Y_concrete),lambda=50,fast=TRUE,M=n,burnin=warmup)
    ch = cbind(ch$beta,ch$sigma2)
    colnames(ch)=c(colnames(X_concrete_scaled),"sigma2")
    
    ch_idx = sim_setting-4 
    x = ch[,ch_idx]
    dhat = momentLS::tune_delta(x,nSplits = 5)$delta *0.8
    r = momentLS::autocov(x)
    return(list(y=r,V="Z",weight_option="unweighted",delta=dhat,delta_grid=dhat))
  }
  stopifnot(sim_setting<=13)
}