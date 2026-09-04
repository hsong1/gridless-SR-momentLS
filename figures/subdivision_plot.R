library(latex2exp)
library(momentLS)
remove(list=ls())
set.seed(157911)
M =10000
y = rnorm(M,sd = 1/sqrt(M))

alphaGrid=makeGrid(upper_threshold=1-1e-6,nX=1001,scale = "equidist")

v=sapply(alphaGrid, function(a) sum(y*a^{0:(length(y)-1)}))

#v=Xtr_cpp(x=alphaGrid,a=x,standardization = FALSE)
pdf("figures/results/x_theta_y_innerproduct.pdf",width = 6,height = 4)


par(mfrow=c(1,1))
par(mar = c(4, 6, 3, 1)) # Bottom, left, top, right margins

plot(alphaGrid, v, type = "l",
     xlab = TeX("$\\theta$"),
     ylab = TeX("$\\sum_k c_k\\, \\theta^{|k|}$"))
     # main =TeX("$\\sum_k c_k\\, \\theta^{|k|}$") ,cex.lab = 1, cex.main = 1)
dev.off()
tail(v)

N=log(10000,base=2)
df=data.frame(k=1:4,Left=1-2^(-((1:4)-1)),Right=1-2^(-(1:4)))
knitr::kable(df)
pdf("figures/results/poly_approx.pdf",width = 8,height = 3)
par(mfrow=c(2,7))
par(mar = c(2.5, 1, 2, 1)) # Bottom, left, top, right margins
for (i in 1:N){
  L=1-2^(-(i-1))
  U=1-2^(-i)
  plot(alphaGrid,v,xlim=c(L,U),type="l",main=paste0("k=",i),xlab="")
}

L=1-2^(-N)
U=1
plot(alphaGrid,v,xlim=c(L,U),type="l",main="[1-1/M,  1]")
dev.off()
