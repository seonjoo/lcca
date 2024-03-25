
#' Generalized Longitudial canonical correlation analysis
#'
#' GLCCA handles more than three sets of data
#' Current version is implemented for only the linear trajectory.
#'
#' @param X object for input list(x, time, J, I, visit)
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
#' @import MASS
#' @import CCP
#' @import CCA

glcca_drcca.linear<-function(X,
               varthresh=0.95,
               projectthresh=1,
               nboot=50,
#               ncores=12,
#               seed.offset=1234,
               verbose=FALSE,
               reg=0,
                nfold=3){

  K = length(X)
  n=nrow(X[[1]])
  re<-lapply(X, function(obj){
    hd_lfpca(obj$X, T=obj$time, J=sum(obj$visit),
             I=obj$I, visit=obj$visit, varthresh=varthresh, projectthresh=projectthresh,
             timeadjust=FALSE,verbose=verbose)
  })

  xilist = lapply(re, function(obj)t(obj$xi) )

  cca = drCCAcombine(xilist,reg=reg,nfold=nfold,nrand=nboot)

  tmp.d = max(1, cca$n)
  loadings = lapply(1:K,
                    function(jj){
                      lcca.0 = re[[jj]]$phix0 %*% (cca$cca_data$white[[jj]]) %*% cca$cca_data$eigvecs[[jj]][,1:tmp.d]
                      lcca.1 = re[[jj]]$phix1 %*% (cca$cca_data$white[[jj]]) %*% cca$cca_data$eigvecs[[jj]][,1:tmp.d]
                      return(list(lcca.0=lcca.0,lcca.1=lcca.1))
                    }
                    )

  scores=cca$proj

  out=list(loadings = loadings, scores=scores, n=n, cca=cca)
  return(out)
}
