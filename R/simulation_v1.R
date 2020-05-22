## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(lcca)
library(MASS)
library(gplots)

sim_lcca <- function(I, r, varthresh){
  ccor.list = rep(0, 100)
  ccor.dim.list = rep(0, 100)
  xcoef_x0.array= array(0, dim=c(144,1,100))
  xcoef_x1.array = array(0, dim=c(144,1,100))
  xcoef_y0.array = array(0, dim=c(81,1,100))
  xcoef_y1.array = array(0, dim=c(81,1,100))
  
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
    
    re = lcca.linear(x=x, y=y, varthresh=varthresh)
    
    ccor.dim.list[i] <- re$ccor.dim
    ccor.list[i] <- re$ccor
    xcoef_x0.array[,,i] <- matrix(re$xcoef_x0, nrow=144, ncol=1)
    xcoef_x1.array[,,i] <- matrix(re$xcoef_x1, nrow=144, ncol=1)
    xcoef_y0.array[,,i] <- matrix(re$xcoef_y0, nrow=81, ncol=1)
    xcoef_y1.array[,,i] <- matrix(re$xcoef_y1, nrow=81, ncol=1)
  }
  
  out=list(ccor.dim.list=ccor.dim.list, ccor.list=ccor.list, 
           xcoef_x0.array=xcoef_x0.array,
           xcoef_x1.array=xcoef_x1.array,
           xcoef_y0.array=xcoef_y0.array,
           xcoef_y1.array=xcoef_y1.array)
  return(out)
}


## I=100, r=0.8, varthresh=0.97
sim1 <- sim_lcca(I=100, r=0.8, varthresh=0.97)

sim1$ccor.dim.list

mean(sim1$ccor.dim.list)
sd(sim1$ccor.dim.list)

mean(sim1$ccor.list)
sd(sim1$ccor.list)

# X
xcoef_x0.array.avg <- apply(abs(sim1$xcoef_x0.array), c(1,2), mean)
xcoef_x1.array.avg <- apply(abs(sim1$xcoef_x1.array), c(1,2), mean)

cor(xcoef_x0.array.avg, phix0)
cor(xcoef_x1.array.avg, phix1)

par(mfrow=c(2,2), mar=rep(0.5,4), bg="gray")
bs=c(-100:100)/1000*1.5
image(phix0, axes=F, col=bluered(200), breaks=bs)
image(xcoef_x0.array.avg, axes=F, col=bluered(200), breaks=bs)
image(phix1, axes=F, col=bluered(200), breaks=bs)
image(xcoef_x1.array.avg, axes=F, col=bluered(200), breaks=bs)

# Y 
xcoef_y0.array.avg <- apply(abs(sim1$xcoef_y0.array), c(1,2), mean)
xcoef_y1.array.avg <- apply(abs(sim1$xcoef_y1.array), c(1,2), mean)

cor(xcoef_y0.array.avg, phiy0)
cor(xcoef_y1.array.avg, phiy1)

par(mfrow=c(2,2), mar=rep(0.5,4), bg="gray")
bs=c(-100:100)/1000*1.5
image(phiy0, axes=F, col=bluered(200), breaks=bs)
image(xcoef_y0.array.avg, axes=F, col=bluered(200), breaks=bs)
image(phiy1, axes=F, col=bluered(200), breaks=bs)
image(xcoef_y1.array.avg, axes=F, col=bluered(200), breaks=bs)




## I=100, r=0.5, varthresh=0.97
sim2 <- sim_lcca(I=100, r=0.5, varthresh=0.97)

## I=100, r=0.2, varthresh=0.97
sim3 <- sim_lcca(I=100, r=0.2, varthresh=0.97)

## I=100, r=0.8, varthresh=0.95
sim4 <- sim_lcca(I=100, r=0.8, varthresh=0.95)

## I=100, r=0.5, varthresh=0.95
sim5 <- sim_lcca(I=100, r=0.5, varthresh=0.95)

## I=100, r=0.2, varthresh=0.95
sim6 <- sim_lcca(I=100, r=0.2, varthresh=0.95)


## I=200, r=0.8, varthresh=0.97
sim2 <- sim_lcca(I=200, r=0.8, varthresh=0.97)

## I=200, r=0.5, varthresh=0.97
sim2 <- sim_lcca(I=200, r=0.5, varthresh=0.97)

## I=200, r=0.2, varthresh=0.97
sim3 <- sim_lcca(I=200, r=0.2, varthresh=0.97)

## I=200, r=0.8, varthresh=0.95
sim4 <- sim_lcca(I=200, r=0.8, varthresh=0.95)

## I=200, r=0.5, varthresh=0.95
sim5 <- sim_lcca(I=200, r=0.5, varthresh=0.95)

## I=200, r=0.2, varthresh=0.95
sim6 <- sim_lcca(I=200, r=0.2, varthresh=0.95)


## I=400, r=0.8, varthresh=0.97
sim2 <- sim_lcca(I=400, r=0.8, varthresh=0.97)

## I=400, r=0.5, varthresh=0.97
sim2 <- sim_lcca(I=400, r=0.5, varthresh=0.97)

## I=400, r=0.2, varthresh=0.97
sim3 <- sim_lcca(I=400, r=0.2, varthresh=0.97)

## I=400, r=0.8, varthresh=0.95
sim4 <- sim_lcca(I=400, r=0.8, varthresh=0.95)

## I=400, r=0.5, varthresh=0.95
sim5 <- sim_lcca(I=400, r=0.5, varthresh=0.95)

## I=400, r=0.2, varthresh=0.95
sim6 <- sim_lcca(I=400, r=0.2, varthresh=0.95)


