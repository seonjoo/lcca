#' Conduct high-dimensional longitudinal principal component analysis
#'
#'
#' For V-by-J dimensioanl matrix of Y, HDLPCA first reduce the dimension of p to J if V>J:
#' Y = UDV. Then, longitudinal PCA is applied to the projected scores of Y on the p-dimensional subspace: Y*=DY.
#' The dimension of the subspace is selected based on the variance explained by first p-components.
#' The default value is projectthresh=1.
#'
#' @param Y p-by-J matrix,
#' @param T Time of theimage collection
#' @param J Total number of observations
#' @param I Total number of subjects
#' @param visit Vector of number of visits per subjects
#' @param verbose (default=FALSE)
#' @param Nx Dimension of the subject-specific components
#' @param Nw Dimension of the subject-visit specific components
#' @param varthresh (default=0.99) Threshold for variance explained for both subject-specific and subject-visit specific compoents for dimension selection
#' @param projectthresh Threshold for variance explain in the first step of SVD
#' @param timeadjust (default=FALSE) Scale time per subject
#' @return xi : Subject-specific principal component scores.
#' @return phix0: Subject-specific principal comonent loadings corresponding to intercept.
#' @return phix1: Subject-specific principal comonent loadings corresponding to slope.
#' @return zeta : Subject-visit-specific principal component scores.
#' @return phiw : Subject-visit-specific principal comonent loadings
#' @examples
#' set.seed(12345678)
#' I=100
#' visit=rpois(I,1)+3
#' time = unlist(lapply(visit, function(x) scale(c(0,cumsum(rpois(x-1,1)+1)))))
#' J = sum(visit)
#' V=2500
#' phix0 = matrix(0,V,3);phix0[1:50,1]<-.1;phix0[1:50 + 50,2]<-.1;phix0[1:50 + 100,3]<-.1
#' phix1 = matrix(0,V,3);phix1[1:50+150,1]<-.1;phix1[1:50 + 200,2]<-.1;phix1[1:50 + 250,3]<-.1
#' phiw = matrix(0,V,3);phiw[1:50+300,1]<-.1;phiw[1:50 + 350,2]<-.1;phiw[1:50 + 400,3]<-.1
#' xi = t(matrix(rnorm(I*3),ncol=I)*c(8,4,2))*3
#' zeta = t(matrix(rnorm(J*3),ncol=J)*c(8,4,2))*2
#' Y = phix0%*% t(xi[rep(1:I, visit),]) + phix1%*% t(time * xi[rep(1:I, visit),]) + phiw %*% t(zeta) + matrix(rnorm(V*J,0,.1),V,J)
#' library(Lpredict)
#' re<-hd_lfpca(Y,T=scale(time,center=TRUE,scale=TRUE),J=J,I=I,visit=visit, varthresh=0.95, projectthresh=1,timeadjust=FALSE,figure=TRUE)
#' cor(phix0, re$phix0)
#' cor(phix1, re$phix1)
#' library(gplots)
#' par(mfrow=c(2,2),mar=rep(0.5,4),bg="gray")
#' bs=c(-100:100)/1000*1.5
#' image(phix0, axes=F,col=bluered(200),breaks=bs)
#' image(re$phix0[,1:3], axes=F,col=bluered(200),breaks=bs)
#' image(phix1, axes=F,col=bluered(200),breaks=bs)
#' image(re$phix1[,1:3], axes=F,col=bluered(200),breaks=bs)
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references Zipunnikov, V., Greven, S., Shou, H., Caffo, B., Reich, D. S., & Crainiceanu, C. (2014). Longitudinal High-Dimensional Principal Components Analysis with Application to Diffusion Tensor Imaging of Multiple Sclerosis. The Annals of Applied Statistics, 8(4), 2175–2202.
#' @keywords hdlfpca glmnet
#' @import MASS
#' @export

hd_lfpca= function(Y,T,J,I,visit, verbose=FALSE, prefix=date(), Nx = NA,Nw = NA,
                   varthresh=0.95,projectthresh=1, timeadjust=FALSE,figure=FALSE){
 # library(MASS)
  	#J: Total number visits out of all subjects (sum of visit)
	#T must be centered+scaled time information
	# ylist : list of subjects.
	###
	### Y = [Y1 Y2... YI] = VDU' (p by J matrix)
	### Calcuate U using iterative calculation of sum(Yi'Yi)
#	Y_total_mean=apply(Y,1,mean)
	Y = round(t(scale(t(Y),center=TRUE,scale=FALSE)), digit=5) # Y is centered.

	if (dim(Y)[1]>=dim(Y)[2]){
	  svdy = svd(t(Y)%*%Y) # There is a sign ambiguity between SVD and eigen function.
	  if(verbose==TRUE){print("Reduce Dimension")}
	  N=J_projected=sum(cumsum(svdy$d^2)/sum(svdy$d^2)<projectthresh) #Dimensionality #
	  if(verbose==TRUE){print(J_projected)}
	  U = svdy$v[,1:J_projected]
	  D = diag(sqrt(svdy$d[1:J_projected]))
	  S = diag(svdy$d[1:J_projected])

	  V = Y %*% U %*% solve(D)
	  Ynew=t(V[,1:J_projected])%*% Y #D %*% t(U)
  	rm(list=c("svdy"))
	}else{
	  svdy = svd(t(Y)%*%Y) # There is a sign ambiguity between SVD and eigen function.
	  if(verbose==TRUE){print("Do not Reduce Dimension")}
	  N=J_projected=nrow(Y)-1#sum(cumsum(svdy$d^2)/sum(svdy$d^2)<projectthresh) #Dimensionality #
	  if(verbose==TRUE){print(J_projected)}
	  U = svdy$v[,1:J_projected]
	  D = diag(sqrt(svdy$d[1:J_projected]))
	  S = diag(svdy$d[1:J_projected])
	  V = Y %*% U %*% solve(D)
	  Ynew=t(V[,1:J_projected])%*% Y #D %*% t(U)
	  rm(list=c("svdy"))
	}


print("hi")
	Yvec=matrix(0,J_projected^2,sum(visit^2))
	X=matrix(0,sum(visit^2),5)
	k=0
	J_sq=0
	for (j in 1:I){
#	  print(j)
		s = 1;
        	Ji = visit[j];
		Ti = T[(1:Ji)+k];
		#Ti=(Ti-mean(Ti))
		if(timeadjust==TRUE){Ti=(Ti-mean(Ti))/sqrt(var(Ti))}
    Yvec_i = matrix(0,J_projected*J_projected,Ji^2)
		X_i=matrix(0,Ji*Ji,5)
		Ytmp = Ynew[,(1:Ji)+k]
    for (j1 in 1:Ji){
    	for (j2 in 1:Ji){
    		Yvec_i[,s] = as.vector(Ytmp[,j1]%*%t(Ytmp[,j2]))
		    X_i[s,] = c(1, Ti[j2], Ti[j1], Ti[j1]*Ti[j2],as.numeric(j1==j2))
    		s = s + 1
		  }
		}
    Yvec[,(1:(visit[j]^2))+J_sq] = Yvec_i;
		X[(1:visit[j]^2)+J_sq,] = X_i
    J_sq = J_sq + Ji^2;
		k = k + Ji;
	}
	if(verbose==TRUE){print("Estimate Covariance Functions")}

#system.time(beta<- Yvec %*% X %*% chol2inv(chol(t(X)%*%X)))
beta <- Yvec %*% X %*%solve(t(X)%*%X)
#sum((beta-beta2)^2)
K00 = matrix(beta[,1],J_projected,J_projected)
K01 = matrix(beta[,2],J_projected,J_projected)
K10 = matrix(beta[,3],J_projected,J_projected)
K11 = matrix(beta[,4],J_projected,J_projected)
Kw = matrix(beta[,5],J_projected,J_projected)
Kx = rbind(cbind(K00,K01),cbind(K10,K11))

Ax=eigen(Kx)
Aw=eigen(Kw)
#Estimate the components
if (is.na(Nx) | is.na(Nw)){
	lim=1;Nw0=Nx0=Nw=Nx=min(J_projected,sum(Ax$values>0),sum(Aw$values>0))

	while (lim>varthresh){
		if (min(Aw$values[1:Nw])<min(Ax$values[1:Nx])){Nw=Nw0-1}else{Nx=Nx0-1}
		lim=(sum(Ax$values[1:Nx])+sum(Aw$values[1:Nw]))/(sum(Ax$values[1:sum(Ax$values>0)])+sum(Aw$values[1:sum(Aw$values>0)]));
		Nx0<-Nx;Nw0<-Nw
#	if (lim>varthresh){Nx=Nx-1;Nw=Nw-1}
		if(verbose==TRUE){print(c(Nx,Nw,lim))}
}
}
Ax0 =  Ax$vectors[1:J_projected,1:Nx]
Ax1 = Ax$vectors[1:J_projected+J_projected,1:Nx]
S_x = Ax$values[1:Nx]
S_u = Aw$values[1:Nw]
Au =  Aw$vectors[,1:Nw]


phix0=V %*% Ax0
phix1=V %*% Ax1
phiw = V %*% Aw$vectors[,1:Nw]

#rm(list=c("Ax","Aw"))

#ramp = colorRamp(c("blue","gray","red"))
#colscheme = rgb( ramp(seq(0, 1, length = 240)), max = 255)#redgreed(255)#
b=rainbow(200)
phix=rbind(phix0,phix1)
#(diag(t(phix) %*% phix))
#(diag(t(phiw) %*% phiw))
#(diag(t(phix0) %*% phix0))
#(diag(t(phix1) %*% phix1))
total_lambda=sum(Ax$values[1:Nx])+sum(Aw$values[1:Nw])
a=(diag(t(phix0) %*% phix0))

if (figure==TRUE){
  x =cbind(Baseline=Ax$values[1:Nx]/total_lambda*a*100,Longitudinal=Ax$values[1:Nx]/total_lambda*(1-a)*100)
  maxval=max(apply(x,1,sum))

  if(Nx<15){h=barplot(t(x[1:Nx,]),border=b[c(100,170)], beside=FALSE,col=b[c(100,170)],names.arg=1:Nx,legend=TRUE,ylim=c(0,1.1*maxval))
    title("Variability Explained by subject-specific components")
    text(h[1:Nx],Ax$values[1:Nx]/total_lambda*100+0.06*maxval,paste(round(Ax$values[1:Nx]/total_lambda*10000)/100, '%',sep=""))
    text(h[1:Nx],Ax$values[1:Nx]/total_lambda*100+0.02*maxval,paste("(",round((1-a[1:Nx])*10000)/100, '%)',sep=""))
}
if(Nx>=15){h=barplot(t(x[1:15,]),border=b[c(100,170)], beside=FALSE,col=b[c(100,170)],names.arg=1:15,legend=TRUE,ylim=c(0,1.1*maxval))
    title("Variability Explained by First 15 Eigenvectors - subject-specific variation")
    text(h[1:15],Ax$values[1:15]/total_lambda*100+0.06*maxval,paste(round(Ax$values[1:15]/total_lambda*10000)/100, '%',sep=""))
    text(h[1:15],Ax$values[1:15]/total_lambda*100+0.02*maxval,paste("(",round((1-a[1:15])*10000)/100, '%)',sep=""))
  }
}

#par(mfrow=c(2,5))
#for (j in 1:5){image(matrix(phix0[,j],200,200),breaks=-120:120/2000,col=colscheme,)}
#for (j in 1:5){image(matrix(phix1[,j],200,200),breaks=-120:120/2000,col=colscheme,)}

## Step 4. EBLUPs
    xi_est = matrix(0,Nx,I);

    zeta_est = matrix(0,Nw,J);

    C00 = t(Ax0)%*%Ax0;
    C11 = t(Ax1)%*%Ax1;
    C01 = t(Ax0)%*%Ax1;
    C10 = t(C01)

    C0W = t(Ax0)%*%Au;
    C1W = t(Ax1)%*%Au;
    if(verbose==TRUE){print("Calcuate Scores")}
    Ut = t(U);
    k = 0;
     for (i in 1:I){
        Ji = visit[i];
        Ti = T[k+1:Ji];
	      if( timeadjust==TRUE ){Ti=(Ti-mean(Ti))/sqrt(var(Ti))}

        Ti_dot = sum(Ti); # it's 0
        Ti_sq_dot = sum(Ti^2); # it's 1
        Ui = Ut[,k+1:Ji];

        ones_Ji = rep(1,Ji);

#   %%%% BiBi' %%%%
        BBi_11 = Ji*C00 + Ti_sq_dot * C11 + Ti_dot * (C10+C01);
        BBi_12 = kronecker(t(ones_Ji),C0W) + kronecker(t(Ti),C1W);
        BBi_21 = t(BBi_12);
        BBi_22 = diag(rep(1,Nw*Ji))

 #       BBi = rbind(cbind(BBi_11, BBi_12),cbind(BBi_21, BBi_22))

       C1=ginv(BBi_11 -BBi_12%*%BBi_22%*%BBi_21)
       C2=BBi_22 + BBi_21%*%ginv(BBi_11-BBi_12 %*% BBi_21 )%*%BBi_12
       C3=ginv(BBi_11)%*%BBi_12%*%C2

#   %%%% Bi'Yi %%%%
        BYi_11 = as.vector(t(Ax0)%*%sqrt(S[1:N,])%*%Ui%*%ones_Ji+ t(Ax1)%*%sqrt(S[1:N,])%*%Ui%*%Ti)
        BYi_12 = as.vector(t(Au)%*%sqrt(S[1:N,])%*%Ui);

        BYi = c(BYi_11,BYi_12);
  #      system.time(wi<-chol2inv(chol(BBi))%*%BYi)
 	wi <- rbind(cbind(C1,-C3),cbind(-t(C3),C2)) %*%BYi
        xi_est[,i] = wi[1:Nx];
        zeta_i_vec = wi[Nx+1:(Ji*Nw)] # put the size
        for (j in 1:Ji){
            zeta_est[,k+j] = zeta_i_vec[((j-1)*Nw+1):(j*Nw)];
	}
        k = k + Ji;
    }


    if(verbose==TRUE){print("Compute Residual to the Demeaned Data.")}
	tmpx=matrix(0,nrow(Ax0),J)
	cumvisit=0
	for (j in 1:I){
		for (k in 1:visit[j]){tmpx[,cumvisit+k]=(Ax0+Ax1*T[cumvisit+k])%*%xi_est[,j]}
		cumvisit=sum(visit[1:j])
	}
	residual=sum((Ynew-tmpx-Aw$vectors[,1:Nw]%*%zeta_est)^2)
	if(verbose==TRUE){print(paste("Residual of LFPCA Model is:",residual))}
	result = list(
	  xi=xi_est,
	  zeta=zeta_est,
	  phix0=phix0,
	  phix1=phix1,
	  Nx=Nx,
	  Nw=Nw,
	  phiw=phiw,
	  sx=Ax$values[1:Nx],
	  sw=Aw$values[1:Nw],
	  residual=residual,
	  tij=T)
	return(result)
}

