source("codes/SR1_V.R")

grid <- momentLS::makeGrid(upper_threshold = 1-1e-6,nX=2001,scale="equidist")
rho0 <- mean(c(grid[(length(grid)-1)], grid[length(grid)]))
y <- rho0^{0:9999}


m_gl <- SR1_gridless_R(y,delta = 1e-6,weight_option = "unweighted",V = "Z")
m_grid <- SR1_V(y = y,V = "Z",alphaGrid = grid)

# objective values
(val_gl <- with(m_gl, L2diff_L2Moment_V(y=y, support, weights, "Z"))) # 0.08822682
(val_grid <- with(m_grid, L2diff_L2Moment_V(y=y, support, weights, "Z"))) # 219.5035
rel_error <- (val_grid-val_gl)/val_gl # 2486.945


# suboptimality bounds
norm1_y <- abs(y[1])+2*sum(abs(y[-1]))
ldp_tol_check = 1e-12/(2*norm1_y)

delta_h_gl = with(m_gl, julia_call("delta_h",y, 1e-6, "unweighted", "Z", support, weights, ldp_tol_check, 1, "root", F))
delta_h_grid = with(m_grid, julia_call("delta_h",y, 1e-6, "unweighted", "Z", support, weights, ldp_tol_check, 1, "root", F))

# estimated delta^h(mu) at mu = m_gl and m_grid
delta_h_gl$delta_h_upper # 7.033705e-11
delta_h_grid$delta_h_upper # 2179.924

-2*min(m_grid$opt$gradient) # max_{theta in grid} (-D_theta^h(mu))

# subopt_bounds
with(m_gl, subopt_bound(y,1e-6,"unweighted","Z",support,weights, ldp_tol = 1e-12/(2*norm1_y))) # 3.141966e-09
with(m_grid, subopt_bound(y,1e-6,"unweighted","Z",support,weights, ldp_tol = 1e-12/(2*norm1_y))) # 97377.53

# plot illustration
# dd(theta) = -2*<x_theta, y - f_mu>

dd_at_theta = function(theta, m){
  sapply(theta, function(t){
    term1 = innerProduct_L2Moment(r = y,support = t,weights =1)
    term2 = innerProduct_twoMoments(support1 = m$support,weights1 = m$weights, support2 = t, weights2 = 1)
  -2*(term1 - term2)
  })
}
# check
all.equal(dd_at_theta(grid,m_grid)/2, as.numeric(m_grid$opt$gradient))


# illustration
grid0 = grid[grid>=0.9965]
finer_grid = seq(min(grid0),max(grid0),length.out=10000)

pdf("experiments/grid_misspecification/directional_derivative.pdf",width = 7,height = 5)
plot(
  finer_grid, dd_at_theta(finer_grid, m_grid),
  type = "l",
  xlab = expression(theta),
  ylab = expression(D[theta]*psi(hat(mu)[grid]))
)
points(grid0, dd_at_theta(grid0, m_grid),pch = 19, col = "blue", cex=1.2)
abline(h = 0, lty = 2, col = "grey50")

theta_min <- delta_h_grid$theta_tilde # based on gridless search

points(theta_min,dd_at_theta(theta_min, m_grid),pch = 19,col = "red", cex=1.2)

legend(
  "bottomleft",
  legend = c("Directional derivative", "Grid points", "Continuum minimizer"),
  col = c("black", "blue", "red"),
  lty = c(1, NA, NA),
  pch = c(NA, 19, 19),
  bty = "n"
)
dev.off()
