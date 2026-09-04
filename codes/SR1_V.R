#remove(list=ls())
library(JuliaCall)
library(momentLS)
julia_setup()
julia_command('using Pkg; Pkg.add(["Clarabel", "SumOfSquares", "DynamicPolynomials", "RCall", "Random", "Distributions", "LinearAlgebra", "LinearSolve", "StatsBase","Polynomials","PolynomialRoots","Quadmath"])')
julia_source("codes/SR1_gridless.jl")
source("codes/supportReduction.R")

compute_XtX_V = function(V, support){
  XtX = outer(support,support)
  if(V=="N"){XtX = 1/(1-XtX)}else{XtX = (1 + XtX)/(1 - XtX)}
  
  return(XtX)
}

compute_Xty_V = function(V, y, support){
  n = length(y)
  Xty = sapply(1:length(support), function(i){
    x = support[i]^(0:(n - 1))
    if(V=="N"){
      out = sum(x * y)
    }else{
      out = (2 * sum(x * y) - y[1]) 
    }
    return(out)
  })
  return(Xty=Xty)
}


SR1_V = function (y, V, delta = NULL, alphaGrid = NULL, init = NULL, tol = 10^-8, 
          maxit = 1000, gradTrace = FALSE, n_alphas = 501, 
          gridControl = list(cm = FALSE, scale = "log"), Xtr_approx = TRUE, precomputed = list(s_alpha = NULL, XtX = NULL, Xtr = NULL, alphaGrid = NULL, input = NULL)) 
{
  n = length(y)
  
  if (!is.null(alphaGrid)) {
    if (!is.null(delta)) {
      warning("provided delta not used")
    }
  }
  if (is.null(alphaGrid)) {
    if (is.null(delta)) {
      stop("at least one of delta or alphaGrid needs to be provided")
    }
    alphaGrid = makeGrid(nX = n_alphas, upper_threshold = 1 - 
                           delta, cm = ifelse(V=="N",T,F), scale = gridControl$scale)
  }
  if (is.null(precomputed$s_alpha)) {
    s_alpha = rep(1,length(alphaGrid))
  }
  else {
    s_alpha = precomputed$s_alpha
  }
  
  if (is.null(precomputed$XtX)) {
    XtX = compute_XtX_V(V, alphaGrid)
  }
  else {
    if (is.null(precomputed$alphaGrid)) {
      stop("when XtX is not NULL, the precomputed list should contain alphaGrid")
    }
    else {
      stopifnot(abs(precomputed$alphaGrid - alphaGrid) < 
                  1e-04)
    }
    XtX = precomputed$XtX
  }
  if (is.null(precomputed$Xty)) {
    Xty = compute_Xty_V(V,y, alphaGrid)
    
  }
  
  else {
    if (is.null(precomputed$input)) {
      stop("When Xty is not NULL, the precomputed list should contain input autocovariance estimate")
    }
    else {
      stopifnot(abs(precomputed$input - y) < 1e-04)
    }
    if (is.null(precomputed$alphaGrid)) {
      stop("When Xty is not NULL, the precomputed list should contain alphaGrid")
    }
    else {
      stopifnot(abs(precomputed$alphaGrid - alphaGrid) < 
                  1e-04)
    }
    stopifnot(length(precomputed$Xty) == length(alphaGrid))
    Xty = precomputed$Xty
  }
  norm2_r = ifelse(V=="N", sum(y^2), 2 * sum(y^2) - y[1]^2)
  tol = tol * sqrt(norm2_r)
  
  if (is.null(init)) {
    inds = integer(0)
    beta = numeric(0)
  }
  else {
    if (any(init$weights < 0)) {
      stop("all weights need to be positive")
    }
    inds = findIndices(support = init$support, alphaGrid = alphaGrid)
    beta = init$weights/s_alpha[inds]
  }
  SR_out = supportReduction(XtX = XtX, Xtr = Xty, s_alpha = s_alpha, 
                            beta = beta, inds = inds, gradTrace = gradTrace, tol = tol, 
                            maxit = maxit)
  inds = SR_out$active_inds
  beta = SR_out$beta
  support = alphaGrid[inds]
  weights = beta/s_alpha[inds]
  weights = weights[order(support)]
  inds = inds[order(support)]
  support = support[order(support)]
  fvals = L2diff_L2Moment_V(y = y, support = support, weights = weights, V)
  res = structure(list(weights = weights, support = support, 
                       active_inds = inds, alphaGrid = alphaGrid, 
                       SRinputList = list(XtX = XtX, Xty = Xty, s_alpha = s_alpha), 
                       opt = list(iter = SR_out$iter, fvals = fvals, convergence = SR_out$convergence, gradient = SR_out$gradient, gradientExt = SR_out$gradientExt), 
                       M = length(y)), class = "SRfit1_V")
  return(res)
}


L2diff_L2Moment_V = function(y,support,weights,V){
  aat = outer(support, support)
  if(V=="N"){
    norm2_m = sum(1/(1 - aat) * outer(weights, weights))
    norm2_y = sum(y^2)
  }else{
    norm2_m = sum((1 + aat)/(1 - aat) * outer(weights, weights))
    norm2_y = 2*sum(y^2)-y[1]^2
  }
  n = length(y)
  innerProduct = sum(sapply(1:length(weights), function(i) {
    x = support[i]^(0:(n - 1))
    if(V=="N"){
      out = sum(x * y) * weights[i]
    }else{
      out = (2 * sum(x * y) - y[1]) * weights[i]
    }
    return(out)
  }))
  f = norm2_y + norm2_m - 2 * innerProduct
  return(f)
}


SR1_gridless_R<-function(y, 
                       delta, 
                       weight_option, 
                       V,
                       maxit=1000,
                       init_supp=0,
                       init_weight=0,
                       tol = 1e-10,
                       opt_type="root",
                       verbose=F){
  
  weight_option = match.arg(weight_option, c("unweighted","weighted"))
  V = match.arg(V,c("Z","N"))
  opt_type = match.arg(opt_type,c("root","SOS") )
  # (r, delta, weight_option, V,  maxit, init_supp, init_weight, tol, opt_type="root", verbose=false)
  out = julia_call("SR1_gridless",y, delta, weight_option, V,maxit,init_supp,init_weight,tol, opt_type,verbose)
  return(out)
}

