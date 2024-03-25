
#' Generalized Longitudial canonical correlation analysis
#'
#' GLCCA handles more than three sets of data
#' Current version is implemented for only the linear trajectory.
#'
#' @param x object for input list(x, time, J, I, visit)
#' @param varthreshold (default=0.95) threshold to detemined the number of components of lpcs.
#' @param projectthresh (default=1) threshold for the dimension reduction projection in lfpca.
#' @param timeadjust (defult=FALSE)
#' @param method (defualt='Wilks') test statistic to be used. "Wilks","Hotelling", "Pillai", or "Roy".
#' @param verbose (default=FALSE) print all details
#' @return ccor
#' @return xcv_x0: Longitudinal Canonical vector for the intercept for x
#' @return xcv_x1: Longitudinal Canonical vector for the slope for x
#' @return xcv_y0: Longitudinal Canonical vector for the intercept for y
#' @return xcv_y1: Longitudinal Canonical vector for the slope for y
#'
#' @export
#'
#' @examples
#' library(lcca)
#'
#' @import MASS
#' @import CCP
#' @import CCA

glcca.linear<-function(X,
               varthresh=0.95,
               projectthresh=1,
               nboot=1000,
               ncores=12,
               seed.offset=1234,
               verbose=FALSE){

  re<-lapply(X, function(obj){
    hd_lfpca(obj$X, T=obj$time, J=sum(obj$visit),
             I=obj$I, visit=obj$visit, varthresh=varthresh, projectthresh=projectthresh,
             timeadjust=FALSE,verbose=verbose)
  })

  xilist = lapply(re, function(obj)t(obj$xi) )
  cca = gCCA(xilist ,0)
  tests = gCCA_select(xilist, nboot=nboot,ncores=ncores,seed.offset=seed.offset )
  tests$r = cca$cor
  ccor.dim=sum(tests$p.value<0.05)
  ccor=cca$cor[1:ccor.dim]
    xcv_x0 = re.X$phix0 %*% cca$xcoef[,1:ccor.dim]
    xcv_x1 = re.X$phix1 %*% cca$xcoef[,1:ccor.dim]
    xcv_y0 = re.Y$phix0 %*% cca$ycoef[,1:ccor.dim]
    xcv_y1 = re.Y$phix1 %*% cca$ycoef[,1:ccor.dim]
    scores=list(x=cca$scores$xscores[,1:ccor.dim],y= cca$scores$yscores[,1:ccor.dim])

    out=list(tests=tests, ccor.dim=ccor.dim, ccor=ccor,
             xcv_x0=xcv_x0,
             xcv_x1=xcv_x1,
             xcv_y0=xcv_y0,
             xcv_y1=xcv_y1,
             scores=scores)
  return(out)
}
