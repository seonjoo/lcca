
#' Longitudial canonical correlation analysis
#'
#' Current version is implemented for only the linear trajectory.
#'
#' @param x object for input list(x, time, J, I, visit)
#' @param varthreshold (default=0.95) threshold to detemined the number of components of lpcs.
#' @param projectthresh (default=1) threshold for the dimension reduction projection in lfpca.
#' @param timeadjust (defult=FALSE)
#' @param method (defualt='Wilks') test statistic to be used. "Wilks","Hotelling", "Pillai", or "Roy".
#' @return ccor
#' @return xcoef_x0: Canonical variable for the intercept for x
#' @return xcoef_x1: Canonical variable for the slope for x
#' @return xcoef_y0: Canonical variable for the intercept for y
#' @return xcoef_y1: Canonical variable for the slope for y
#'
#' @export
#'
#' @examples
#' set.seed(12345678)
#' r=0.8
#' mu = c(0,0,0,0,0,0)
#' stddev = rep(c(8,4,2),2)
#' cormatx = diag(1,6,6)
#' cormatx[1,5] <- r
#' cormatx[5,1] <- r
#' covmatx = stddev %*% t(stddev) * cormatx
#'
#' ## Generate scores
#' xi = mvrnorm(n = 100, mu = mu, Sigma = covmatx, empirical = FALSE)
#' I=100
#'
#' ## X
#' visit.X =rpois(I,1)+3
#' time.X = unlist(lapply(visit.X, function(x) scale(c(0,cumsum(rpois(x-1,1)+1)))))
#' J.X = sum(visit.X)
#' xi.X = xi[,1:3]
#' V.x=144
#' phix0 = matrix(0,V.x,3); phix0[1:12, 1]<-.1; phix0[1:12 + 12, 2]<-.1; phix0[1:12 + 12*2, 3]<-.1
#' phix1 = matrix(0,V.x,3); phix1[1:12 + 12*3, 1]<-.1; phix1[1:12 + 12*4, 2]<-.1; phix1[1:12 + 12*5, 3]<-.1
#' phixw = matrix(0,V.x,3); phixw[1:12 + 12*6, 1]<-.1; phixw[1:12 + 12*7, 2]<-.1; phixw[1:12 + 12*8, 3]<-.1
#' zeta.X = t(matrix(rnorm(J.X*3), ncol=J.X)*c(8,4,2))*2
#' X = phix0 %*% t(xi.X[rep(1:I, visit.X),]) + phix1 %*% t(time.X * xi.X[rep(1:I, visit.X),]) + phixw %*% t(zeta.X) + matrix(rnorm(V.x*J.X, 0, .1), V.x, J.X)
#'
#' ## Y
#' visit.Y=rpois(I,1)+3
#' time.Y = unlist(lapply(visit.Y, function(x) scale(c(0,cumsum(rpois(x-1,1)+1)))))
#' K.Y = sum(visit.Y)
#'
#' V.y=81
#' phiy0 = matrix(0,V.y,3); phiy0[1:9, 1]<-.1; phiy0[1:9 + 9, 2]<-.1; phiy0[1:9 + 9*2, 3]<-.1
#' phiy1 = matrix(0,V.y,3); phiy1[1:9 + 9*3, 1]<-.1; phiy1[1:9 + 9*4, 2]<-.1; phiy1[1:9 + 9*5, 3]<-.1
#' phiyw = matrix(0,V.y,3); phiyw[1:9 + 9*6, 1]<-.1; phiyw[1:9 + 9*7, 2]<-.1; phiyw[1:9 + 9*8, 3]<-.1
#' zeta.Y = t(matrix(rnorm(K.Y*3), ncol=K.Y)*c(8,4,2))*2
#' xi.Y = xi[,4:6]
#' Y = phiy0 %*% t(xi.Y[rep(1:I, visit.Y),]) + phiy1 %*% t(time.Y * xi.Y[rep(1:I, visit.Y),]) + phiyw %*% t(zeta.Y) + matrix(rnorm(V.y*K.Y ,0, .1), V.y, K.Y)
#' x = list(X=X, time=time.X, I=I, J=sum(visit.X),visit=visit.X)
#' y = list(X=Y, time=time.Y, I=I, J=sum(visit.Y),visit=visit.Y)
#' lcca.linear(x=x,y=y)
#' @import MASS, CCP, CCA

lcca.linear<-function(x,y,
               varthresh=0.95,
               projectthresh=1,
               method='Wilks'){

  re.X <- hd_lfpca(x$X, T=x$time, J=sum(x$visit), I=x$I, visit=x$visit, varthresh=varthresh, projectthresh=projectthresh, timeadjust=FALSE, figure=TRUE)

  re.Y <- hd_lfpca(y$X, T=y$time, J=y$J, I=y$I, visit=y$visit, varthresh=varthresh, projectthresh=projectthresh, timeadjust=FALSE, figure=TRUE)

  cca = cc(t(re.X$xi), t(re.Y$xi))
  tests = p.asym(rho = cca$cor, N = x$I,p = re.X$Nx,q = re.Y$Nx , tstat = method)
  ccor.dim=sum(tests$p.value<0.05)
  ccor=cca$cor[1:ccor.dim]
    xcv_x0 = re.X$phix0 %*% cca$scores$corr.X.xscores[,1:ccor.dim]
    xcv_x1 = re.X$phix1 %*% cca$scores$corr.X.xscores[,1:ccor.dim]
    xcv_y0 = re.Y$phix0 %*% cca$scores$corr.Y.yscores[,1:ccor.dim]
    xcv_y1 = re.Y$phix1 %*% cca$scores$corr.Y.yscores[,1:ccor.dim]
    scores=list(x=cca$scores$xscores[,1:ccor.dim],y= cca$scores$yscores[,1:ccor.dim])

    out=list(tests=tests, ccor.dim=ccor.dim, ccor=ccor,
             xcv_x0=xcv_x0,
             xcv_x1=xcv_x1,
             xcv_y0=xcv_y0,
             xcv_y1=xcv_y1,
             scores=scores)
  return(out)
}
