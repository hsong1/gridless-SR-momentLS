# a minimal fitting example

# simulate a mixture of geometric distribution
n=1000
set.seed(12345)
weights <- c(0.25, 0.50, 0.25)
thetas  <- c(0.10, 0.30, 0.80)
components <- sample(1:3, size = n, replace = TRUE, prob = weights)
x <- rgeom(n, prob = thetas[components])
y = table(factor(x,levels=0:max(x)))/length(x)

# fit SR gridless
source("codes/SR1_V.R")
m_gl = SR1_gridless_R(y = y,delta = 0,weight_option = "weighted",V = "N",maxit = 1000,tol = 1e-10)

# fit SR grid
m_grid = SR1_V(y = y,V = "N",delta = 1e-6,n_alphas = 1001,tol = 1e-7)

# objective value comparison psi(m_grid) - psi(m_gl)
L2diff_L2Moment_V(y,m_grid$support,m_grid$weights,"N") - L2diff_L2Moment_V(y,m_gl$support,m_gl$weights,"N") 
# 1.359312e-08

# suboptimality bound 
subopt_bound(y,delta = 0,weight_option = "weighted",V = "N",supp = m_gl$support,weight = m_gl$weights,ldp_tol=1e-12)
# psi(m_gl) - psi(mu_hat) <= 2.504134e-10
subopt_bound(y,delta = 0,weight_option = "weighted",V = "N",supp = m_grid$support,weight = m_grid$weights,ldp_tol=1e-12)
# psi(m_grid) - psi(mu_hat) <= 3.780282e-07