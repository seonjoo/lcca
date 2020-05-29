## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(lcca)
library(MASS)
library(gplots)
library(tidyverse)

sim_lcca <- function(I, r){
  ccor.list = rep(0, 100)
  ccor.dim.list = rep(0, 100)
  xcv_x0.array= array(0, dim=c(144,1,100))
  xcv_x1.array = array(0, dim=c(144,1,100))
  xcv_y0.array = array(0, dim=c(81,1,100))
  xcv_y1.array = array(0, dim=c(81,1,100))

  cv_x1s = rep(0,100)
  cv_x1s = rep(0,100)
  cv_y0s = rep(0,100)
  cv_y1s = rep(0,100)

  set.seed(12345678)

  for(i in (1:100)){
    ## -----------------------------------------------------------------------------
    mu = c(0,0,0,0,0,0)
    stddev = rep(c(8,4,2),2)
    cormatx = diag(1,6,6)
    cormatx[1,5] <- r
    cormatx[5,1] <- r
    covmatx = stddev %*% t(stddev) * cormatx

    ## Generate scores
    xi = mvrnorm(n = I, mu = mu, Sigma = covmatx, empirical = FALSE)

    ## X
    visit.X =rpois(I,1)+3
    time.X = unlist(lapply(visit.X, function(x) scale(c(0,cumsum(rpois(x-1,1)+1)))))
    J.X = sum(visit.X)
    xi.X = xi[,1:3]
    V.x=144
    phix0 = matrix(0,V.x,3); phix0[1:12, 1]<-.1; phix0[1:12 + 12, 2]<-.1; phix0[1:12 + 12*2, 3]<-.1
    phix1 = matrix(0,V.x,3); phix1[1:12 + 12*3, 1]<-.1; phix1[1:12 + 12*4, 2]<-.1; phix1[1:12 + 12*5, 3]<-.1
    phixw = matrix(0,V.x,3); phixw[1:12 + 12*6, 1]<-.1; phixw[1:12 + 12*7, 2]<-.1; phixw[1:12 + 12*8, 3]<-.1
    zeta.X = t(matrix(rnorm(J.X*3), ncol=J.X)*c(8,4,2))*2
    X = phix0 %*% t(xi.X[rep(1:I, visit.X),]) + phix1 %*% t(time.X * xi.X[rep(1:I, visit.X),]) + phixw %*% t(zeta.X) + matrix(rnorm(V.x*J.X, 0, .1), V.x, J.X)

    ## Y
    visit.Y=rpois(I,1)+3
    time.Y = unlist(lapply(visit.Y, function(x) scale(c(0,cumsum(rpois(x-1,1)+1)))))
    K.Y = sum(visit.Y)

    V.y=81
    phiy0 = matrix(0,V.y,3); phiy0[1:9, 1]<-.1; phiy0[1:9 + 9, 2]<-.1; phiy0[1:9 + 9*2, 3]<-.1
    phiy1 = matrix(0,V.y,3); phiy1[1:9 + 9*3, 1]<-.1; phiy1[1:9 + 9*4, 2]<-.1; phiy1[1:9 + 9*5, 3]<-.1
    phiyw = matrix(0,V.y,3); phiyw[1:9 + 9*6, 1]<-.1; phiyw[1:9 + 9*7, 2]<-.1; phiyw[1:9 + 9*8, 3]<-.1
    zeta.Y = t(matrix(rnorm(K.Y*3), ncol=K.Y)*c(8,4,2))*2
    xi.Y = xi[,4:6]
    Y = phiy0 %*% t(xi.Y[rep(1:I, visit.Y),]) + phiy1 %*% t(time.Y * xi.Y[rep(1:I, visit.Y),]) + phiyw %*% t(zeta.Y) + matrix(rnorm(V.y*K.Y ,0, .1), V.y, K.Y)

    ## -----------------------------------------------------------------------------

    x = list(X=X, time=time.X, I=I, J=sum(visit.X),visit=visit.X)
    y = list(X=Y, time=time.Y, I=I, J=sum(visit.Y),visit=visit.Y)

    re = lcca.linear(x=x, y=y)

    ccor.dim.list[i] <- re$ccor.dim
    ccor.list[i] <- re$ccor
    xcv_x0.array[,,i] <- matrix(re$xcv_x0, nrow=144, ncol=1)
    xcv_x1.array[,,i] <- matrix(re$xcv_x1, nrow=144, ncol=1)
    xcv_y0.array[,,i] <- matrix(re$xcv_y0, nrow=81, ncol=1)
    xcv_y1.array[,,i] <- matrix(re$xcv_y1, nrow=81, ncol=1)

    cv_x1s[i] <- abs(cor(xcv_x0.array[,,i], phix0)[,1])
    cv_x1s[i] <- abs(cor(xcv_x1.array[,,i], phix1)[,1])

    cv_y0s[i] <- abs(cor(xcv_y0.array[,,i], phiy0)[,2])
    cv_y1s[i] <- abs(cor(xcv_y1.array[,,i], phiy1)[,2])
  }

  out=list(ccor.dim.list=ccor.dim.list, ccor.list=ccor.list,
           cv_x1s=cv_x1s,
           cv_x1s=cv_x1s,
           cv_y0s=cv_y0s,
           cv_y1s=cv_y1s)
  return(out)
}

get_bias = function(estimate, truth) {
  mean(estimate) - truth
}

## I=50, r=0.8
sim11 <- sim_lcca(I=50, r=0.8)
get_bias(mean(sim11$ccor.list), 0.8)

## I=50, r=0.5
sim12 <- sim_lcca(I=50, r=0.5)
get_bias(mean(sim12$ccor.list), 0.5)

## I=50, r=0.3
sim13 <- sim_lcca(I=50, r=0.3)
get_bias(mean(sim13$ccor.list), 0.3)

## I=50, r=0.1
sim14 <- sim_lcca(I=50, r=0.1)
get_bias(mean(sim14$ccor.list), 0.1)


## I=100, r=0.8
sim15 <- sim_lcca(I=100, r=0.8)
get_bias(mean(sim15$ccor.list), 0.8)
## I=100, r=0.5
sim16 <- sim_lcca(I=100, r=0.5)
get_bias(mean(sim16$ccor.list), 0.5)
## I=100, r=0.3
sim17 <- sim_lcca(I=100, r=0.3)
get_bias(mean(sim17$ccor.list), 0.3)
## I=100, r=0.1
sim18 <- sim_lcca(I=100, r=0.1)
get_bias(mean(sim18$ccor.list), 0.1)


## I=200, r=0.8
sim19 <- sim_lcca(I=200, r=0.8)
get_bias(mean(sim19$ccor.list), 0.8)
## I=200, r=0.5
sim110 <- sim_lcca(I=200, r=0.51) #error when r=.5
get_bias(mean(sim110$ccor.list), 0.51)
## I=200, r=0.3
sim111 <- sim_lcca(I=200, r=0.31) #error when r=.3
get_bias(mean(sim111$ccor.list), 0.31)
## I=200, r=0.1
sim112 <- sim_lcca(I=200, r=0.1)
get_bias(mean(sim112$ccor.list), 0.1)


## I=400, r=0.8
sim113 <- sim_lcca(I=400, r=0.8)
get_bias(mean(sim113$ccor.list), 0.8)
## I=400, r=0.5
sim114 <- sim_lcca(I=400, r=0.5)
get_bias(mean(sim114$ccor.list), 0.5)
## I=400, r=0.3
sim115 <- sim_lcca(I=400, r=0.3)
get_bias(mean(sim115$ccor.list), 0.3)
## I=400, r=0.1
sim116 <- sim_lcca(I=400, r=0.1)
get_bias(mean(sim116$ccor.list), 0.1)



dd <- data.frame(ccors = c(sim11$ccor.list, sim12$ccor.list, sim13$ccor.list, sim14$ccor.list,
                           sim15$ccor.list, sim16$ccor.list, sim17$ccor.list, sim18$ccor.list,
                           sim19$ccor.list, sim110$ccor.list, sim111$ccor.list, sim112$ccor.list,
                           sim113$ccor.list, sim114$ccor.list, sim115$ccor.list, sim116$ccor.list),
                 r = factor(rep(c(0.8,0.5,0.3,0.1), each=100)),
                 i = factor(rep(c(50,100,200,400), each=400)))

ggplot(dd , aes(x=r, y=ccors)) +
  geom_boxplot() +
  scale_y_continuous(name = "Estimated Canonical Correlations",
                     breaks = seq(0,1,.2),
                     limits=c(0, 1)) +
  scale_x_discrete(name = "Canonical Correlations") +
  ggtitle("Boxplot of estimated canonical correlations") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) +
  facet_grid(. ~ i)


dd <- data.frame(cv_x0 = c(sim11$cv_x0s, sim12$cv_x0s, sim13$cv_x0s, sim14$cv_x0s,
                           sim15$cv_x0s, sim16$cv_x0s, sim17$cv_x0s, sim18$cv_x0s,
                           sim19$cv_x0s, sim110$cv_x0s, sim111$cv_x0s, sim112$cv_x0s,
                           sim113$cv_x0s, sim114$cv_x0s, sim115$cv_x0s, sim116$cv_x0s),
                 r = factor(rep(c(0.8,0.5,0.3,0.1), each=100)),
                 i = factor(rep(c(50,100,200,400), each=400)))

ggplot(dd , aes(x=r, y=cv_x0)) +
  geom_boxplot() +
  scale_y_continuous(name = "Estimated CVs * True CVs X intercept",
                     breaks = seq(0,1,.2),
                     limits=c(0, 1)) +
  scale_x_discrete(name = "Canonical Correlations") +
  ggtitle("Boxplot of Inner product of Estimated CVs and True CVs X intercept") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) +
  facet_grid(. ~ i)

dd <- data.frame(cv_x1 = c(sim11$cv_x1s, sim12$cv_x1s, sim13$cv_x1s, sim14$cv_x1s,
                           sim15$cv_x1s, sim16$cv_x1s, sim17$cv_x1s, sim18$cv_x1s,
                           sim19$cv_x1s, sim110$cv_x1s, sim111$cv_x1s, sim112$cv_x1s,
                           sim113$cv_x1s, sim114$cv_x1s, sim115$cv_x1s, sim116$cv_x1s),
                 r = factor(rep(c(0.8,0.5,0.3,0.1), each=100)),
                 i = factor(rep(c(50,100,200,400), each=400)))

ggplot(dd , aes(x=r, y=cv_x1)) +
  geom_boxplot() +
  scale_y_continuous(name = "Estimated CVs * True CVs X slope",
                     breaks = seq(0,1,.2),
                     limits=c(0, 1)) +
  scale_x_discrete(name = "Canonical Correlations") +
  ggtitle("Boxplot of Inner product of Estimated CVs and True CVs X slope") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) +
  facet_grid(. ~ i)



dd <- data.frame(cv_y0 = c(sim11$cv_y0s, sim12$cv_y0s, sim13$cv_y0s, sim14$cv_y0s,
                           sim15$cv_y0s, sim16$cv_y0s, sim17$cv_y0s, sim18$cv_y0s,
                           sim19$cv_y0s, sim110$cv_y0s, sim111$cv_y0s, sim112$cv_y0s,
                           sim113$cv_y0s, sim114$cv_y0s, sim115$cv_y0s, sim116$cv_y0s),
                 r = factor(rep(c(0.8,0.5,0.3,0.1), each=100)),
                 i = factor(rep(c(50,100,200,400), each=400)))

ggplot(dd , aes(x=r, y=cv_y0)) +
  geom_boxplot() +
  scale_y_continuous(name = "Estimated CVs * True CVs Y intercept",
                     breaks = seq(0,1,.2),
                     limits=c(0, 1)) +
  scale_x_discrete(name = "Canonical Correlations") +
  ggtitle("Boxplot of Inner product of Estimated CVs and True CVs Y intercept") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) +
  facet_grid(. ~ i)

dd <- data.frame(cv_y1 = c(sim11$cv_y1s, sim12$cv_y1s, sim13$cv_y1s, sim14$cv_y1s,
                           sim15$cv_y1s, sim16$cv_y1s, sim17$cv_y1s, sim18$cv_y1s,
                           sim19$cv_y1s, sim110$cv_y1s, sim111$cv_y1s, sim112$cv_y1s,
                           sim113$cv_y1s, sim114$cv_y1s, sim115$cv_y1s, sim116$cv_y1s),
                 r = factor(rep(c(0.8,0.5,0.3,0.1), each=100)),
                 i = factor(rep(c(50,100,200,400), each=400)))

ggplot(dd , aes(x=r, y=cv_y1)) +
  geom_boxplot() +
  scale_y_continuous(name = "Estimated CVs * True CVs Y slope",
                     breaks = seq(0,1,.2),
                     limits=c(0, 1)) +
  scale_x_discrete(name = "Canonical Correlations") +
  ggtitle("Boxplot of Inner product of Estimated CVs and True CVs Y slope") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) +
  facet_grid(. ~ i)

save.image()
