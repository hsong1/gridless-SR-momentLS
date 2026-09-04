# Slightly modified version of rinvgauss function from statmod package
rinvgauss. <- function(n,mean=1,shape=NULL,dispersion=1){
  if(!is.null(shape)){dispersion <- 1/shape}
  mu <- rep_len(mean,n)
  phi <- rep_len(dispersion,n)
  r <- rep_len(0,n)
  i <- (mu>0 & phi>0)
  if(!all(i)){
    r[!i] <- NA
    n <- sum(i)
  }
  phi[i] <- phi[i]*mu[i]
  Y <- rchisq(n,df=1)
  X1 <- 1+phi[i]/2*(Y-sqrt(4*Y/phi[i]+Y^2))
  # Note: The line above should yield all X1>0, but it occasionally doesn’t due to
  # numerical precision issues.  The line below detects this and recomputes
  # the relevant elements of X1 using a 2nd-order Taylor expansion of the
  # sqrt function, which is a good approximation whenever the problem occurs.
  if(any(X1<=0)){X1[X1<=0] <- (1/(Y*phi[i]))[X1<=0]}
  firstroot <- as.logical(rbinom(n,size=1L,prob=1/(1+X1)))
  r[i][firstroot] <- X1[firstroot]
  r[i][!firstroot] <- 1/X1[!firstroot]
  mu*r
}

# Draw from inverse-Gaussian distribution while avoiding potential numerical problems
rinvgaussian <- function(n,m,l){
  m. <- m/sqrt(m*l)
  l. <- l/sqrt(m*l)
  sqrt(m*l)*rinvgauss.(n,m.,l.)
} 


# Gibbs iteration functions for both Bayesian lassos
# Note: The versions have separate functions, as opposed to being different
#       options of the same function, since the latter would require checking any such
#       options every time the function is called, i.e., in every MCMC iteration.
# Note: The values of XTY, n, and p can obviously be calculated from X and Y, but they
#       are included as inputs to avoid recalculating them every time the function is
#       called, i.e., in every MCMC iteration.
iter.bl.original <- function(beta,sigma2,X,Y,XTY,n,p,lambda){
  d.tau.inv <- rinvgaussian(p,sqrt(lambda^2*sigma2/beta^2),lambda^2)
  A.chol <- chol(t(X)%*%X+diag(d.tau.inv))
  beta.tilde <- backsolve(A.chol,backsolve(t(A.chol),XTY,upper.tri=F))
  Z <- rnorm(p)
  beta.new <- beta.tilde+sqrt(sigma2)*backsolve(A.chol,Z)
  sigma2.new <- (sum((Y-drop(X%*%beta.new))^2)+
                   sum(beta.new^2*d.tau.inv))/rchisq(1,n+p-1)
  return(list(beta=beta.new,sigma2=sigma2.new))
}
iter.bl.fast <- function(beta,sigma2,X,Y,XTY,n,p,lambda){
  d.tau.inv <- rinvgaussian(p,sqrt(lambda^2*sigma2/beta^2),lambda^2)
  A.chol <- chol(t(X)%*%X+diag(d.tau.inv))
  beta.tilde <- backsolve(A.chol,backsolve(t(A.chol),XTY,upper.tri=F))
  sigma2.new <- (sum(Y^2)-sum(XTY*beta.tilde))/rchisq(1,n-1)
  Z <- rnorm(p)
  beta.new <- beta.tilde+sqrt(sigma2.new)*backsolve(A.chol,Z)
  return(list(beta=beta.new,sigma2=sigma2.new))
}
# Run original and two-step Bayesian lassos
run.bl <- function(X,Y,lambda,M,burnin=1000,fast=TRUE){
  XTY <- drop(t(X)%*%Y)
  n <- dim(X)[1]
  p <- dim(X)[2]
  iter.bl <- get(paste("iter.bl.",ifelse(fast,"fast","original"),sep=""))
  # chaindigits <- max(1,ceiling(log(M,10)))  # digits needed for chain label strings
  # for(chain in 0:(M-1)){
  beta <- rep(1,p)
  sigma2 <- 1
  # chaintext <- substring(format(chain/(10^chaindigits),nsmall=chaindigits),3)
  # outfile.beta <- paste(outfile.stem,"-",chaintext,"-b.txt",sep="")
  # outfile.sigma2 <- paste(outfile.stem,"-",chaintext,"-s.txt",sep="")
  # if(write.each){
  #   for(k in 1:K){
  #     iter.result <- iter.bl(beta,sigma2,X,Y,XTY,n,p,lambda)
  #     beta <- iter.result$beta
  #     sigma2 <- iter.result$sigma2
  #     if(keep.beta){
  #       beta.row <- matrix(beta,nrow=1)
  #       write.table(beta.row,outfile.beta,append=T,row.names=F,col.names=F)
  #     }
  #     write.table(sigma2,outfile.sigma2,append=T,row.names=F,col.names=F)
  #   }
  # }else{
  beta.chain <- matrix(NA,nrow=M,ncol=p)
  sigma2.chain <- rep(NA,M)
  for(k in 1:(burnin+M)){
    iter.result <- iter.bl(beta,sigma2,X,Y,XTY,n,p,lambda)
    beta <- iter.result$beta
    sigma2 <- iter.result$sigma2
    if (k>burnin){
      beta.chain[k-burnin,] <- beta
      sigma2.chain[k-burnin] <- sigma2
    }
  }
  
  return(list(beta=beta.chain,sigma2=sigma2.chain))
  #     if(keep.beta){
  #       write.table(beta.chain,outfile.beta,row.names=F,col.names=F)
  #     }
  #     write.table(sigma2.chain,outfile.sigma2,row.names=F,col.names=F)
  #   }
  #   typetext <- ifelse(fast,"Fast","Orig")
  #   print(paste(typetext,"chain",chain+1,"of",M,"complete at",date()))
  #   flush.console()
  # }
}
