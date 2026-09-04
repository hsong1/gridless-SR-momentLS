remove(list=ls())

save_ldp_comp = function(d_random, d_rational, filename){
  
  pdf(filename, width =10, height = 3.5)
  tol = d_random$tol
  par(mfrow = c(1,3))
  
  # plot 1
  yr <- range(c(d_random[,"diff"], d_rational[,"diff"]))
  plot(tol, d_random[,"diff"], log="xy", type="b", pch=19, lty=1,
       ylim = yr,
       xlab=expression(epsilon), ylab="Max approximation error",
       main="(a)")
  lines(tol, d_rational[,"diff"], type="b", pch=17, lty=2)
  #abline(0,1, lty=5)
  grid()
  legend("topleft",
         legend=c("Random coeff", "Trunc rational"),
         pch=c(19,17), lty=c(1,2), bty="n")
  
  # plot 2
  yr <- range(c(d_random[,"normalized_diff"], d_rational[,"normalized_diff"], 1))
  plot(tol, d_random[,"normalized_diff"], log="xy", type="b", pch=19, lty=1,
       ylim = yr,
       xlab=expression(epsilon), ylab="Normalized max approximation error",
       main="(b)")
  lines(tol, d_rational[,"normalized_diff"], type="b", pch=17, lty=2)
  abline(h=1, lty=5)
  
  grid()
  legend("topright",
         legend=c("Random coeff", "Trunc rational", "y=1"),
         pch=c(19,17,NA), lty=c(1,2,5), bty="n")
  
  # plot 3
  yr <- range(c(d_random[,"maxdegree"], d_rational[,"maxdegree"]))
  plot(tol, d_random[,"maxdegree"], log="x", type="b", pch=19, lty=1,
       ylim = yr,
       xlab=expression(epsilon), ylab=expression(paste("maxdegree (",epsilon, ")")),
       main="(c)")
  lines(tol, d_rational[,"maxdegree"], type="b", pch=17, lty=2)
  grid()
  legend("topright",
         legend=c("Random coeff", "Trunc rational"),
         pch=c(19,17), lty=c(1,2), bty="n")
  
  par(mfrow = c(1,1))
  dev.off()
}

load("experiments/ldp_approximation_accuracy/results/ldp_approximation_results.Rdata")
save_ldp_comp(d_random, d_rational, "experiments/ldp_approximation_accuracy/results/ldp_comp.pdf")
