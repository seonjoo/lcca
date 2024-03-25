#' Dimension selection for gCCA using random permutation test
#'
#' @param X : input data list
#' @param nboot (default 1000) number of permutation
#' @param ncores (default 12) number of cores for permutation test
#' @param seed.offset (default 1234)
#'
#' @return glcca.dim : selected dimension
#' @return dat bootstrapped data
#' @import parallel
#' @import MASS
#' @export
#'
#' @examples
#' ndim=50
#' L1=L2=rep(0,ndim)
#' L1[1:5]<-1
#' L2[6:10]<-1
#' L = cbind(c(L1,L1,L1),c(L2,L2,rep(0, ndim)))#;L=L/sqrt(sum(L*L))
#' sigma =0.5*as.matrix(L) %*% t(as.matrix(L));diag(sigma)<-1
#' n=300
#' set.seed(123)
#' Xtmp = mvrnorm(n,mu=rep(0,ndim*3), Sigma=sigma)
#' X = list(Xtmp[,1:ndim],Xtmp[,1:ndim+ndim],Xtmp[,1:ndim+(2*ndim)])
#' fit.gcca <- gCCA(lapply(X, scale),0)
#' gCCA_select(X, nboot=300,ncores=12,seed.offset=1234)
#'
gCCA_select <-function(X,
                       nboot=1000,
                       ncores=12,
                       seed.offset=1234){



  N=nrow(X[[1]])
  ps = unlist(lapply(X, ncol))
  dims = min(N, min(ps))
  stat <- do.call(cbind,mclapply(1:nboot,
                                 function(seednum){
                                   set.seed(seednum + seed.offset)
                                   ind = 1:N
                                   ps = unlist(lapply(X, ncol))
                                   dims = min(N, min(ps))
                                   Xnew = lapply(1:length(X),
                                                 function(jj){
                                                   new<-X[[jj]]
                                                   if(jj>1){new<-X[[jj]][sample(ind, size = N, replace = FALSE),]}
                                                   return(new)})

                                   testnew =gCCA(lapply(Xnew, function(x)as.matrix(scale(x))),0)

                                   return(testnew$eigval[1:dims])
                                 }, mc.cores=ncores))


  pvalue=1-apply(stat < test$eigval[1:dims],1,sum)/nboot

  comput_fdr_alpha<-function(pvalue,q=0.05){
    m<-length(pvalue)
    pvalue.sort = sort(pvalue)
    R=max(which(pvalue.sort < (1:m)/m*q))
    return(1-R*q/m)
  }
  ufdr95 = apply(stat,1,  function(ss)quantile(ss, prob=comput_fdr_alpha(pvalue)))
  lfdr95 = apply(stat,1,  function(ss)quantile(ss, prob=1-comput_fdr_alpha(pvalue)))

  dat = data.frame(dim=1:length(test$eigval[1:dims]),
                   pvalue = pvalue,
                   gcorr = test$eigval[1:dims],
                   u95=ufdr95,
                   l95=lfdr95,
                   p.fdr = p.adjust(pvalue, method='BH'))

  glcca.dim <-0
  ii=1
  stop=0
  while(ii < nrow(dat) & stop==0 ){
    if (dat$p.fdr[ii]<0.05){
      glcca.dim<-ii
      stop=0
    }else{
      stop=1
    }
    ii <- ii+1
  }
  return(list(glcca.dim,dat))
  }
