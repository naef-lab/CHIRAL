#' Infer the circular phase of gene expression samples
#'
#' @param E Required matrix of gene expression. Samples should be on the columns, genes on the rows. Rows should be named to perform the inference on a subset of the matrix.
#' @param iterations Number of maximum iterations, default 500.
#' @param clockgenes Set of clock genes (or relevant genes) to be used for the inference. Should be a of the form c("gene1", "gene2") and subset of the rownames, default NULL.
#' @param tau2 Tau parameter for the prior on the gene coefficient, we suggest to leave the default. It is automatically calculated by the algorithm, can be changed if needed. Default NULL.
#' @param u u parameter for the prior on the gene means, we suggest to leave the default. It is automatically calculated by the algorithm, can be changed if needed. Default NULL.
#' @param sigma2 standard deviation of the data points fro the prediction, we suggest to leave the default. It is automatically calculated by the algorithm, can be changed if needed. default NULL.
#' @param TSM Switches the two state model in the EM on and off, default TRUE.
#' @param mean.centre.E parameter to decide if center around the empirical mean of the data, default TRUE.
#' @param pbar shows the progress bar on max number of iterations, default TRUE.
#' @param phi.start set of initial guess for the phases, default NULL.
#' @param standardize standardize the matrix used for the inference, default FALSE.
#' @param GTEX_names To set to true only if analyzing GTEx data, converts row names from ENSGXXX_GeneName|XX to GeneName. Default FALSE
#' @return A list containing: inferred phases in [0,2*pi], inferred standard deviation of the data, gene parameters gene weights (relevant with two state model),
#' @return iteration of stop, inferred standard deviation of data of rhythmic genes, the input matrix, history of the Q function of the EM, genes selected, matrix used for inference.
#' @examples
#' CHIRAL(data_exon)
#' CHIRAL(data_exon, iterations=50, clockgenes=gene_inf, standardize=TRUE)
#' @export



CHIRAL<- function(E, iterations=500, clockgenes=NULL,tau2=NULL, u=NULL,
                  sigma2=NULL, TSM=TRUE, mean.centre.E=TRUE, q=0.1, update.q=FALSE, pbar=TRUE,
                  phi.start=NULL, standardize=FALSE, GTEx_names=FALSE){

  # E is the data matrix
  #on the rows there are genes, rows should be named
  #on the columns there are samples
  #be sure that what you pass as "clockgenes" is a subset of the rownames

  id=as.vector(c(1,1,1,1))
  E=as.matrix(E)


  E=E[ , colSums(is.na(E)) == 0]
  if(standardize==TRUE){E=sweep(E,1,apply(E,1,sd),FUN="/")}
  E.full=E

  ### Clock gene selection ###

  rownames(E)=toupper(rownames(E))
  if(is.null(clockgenes)){
    if(is.null(rownames(E))){
      rownames(E)=1:nrow(E)
    }
    clockgenes=rownames(E)
  }
  clockgenes=toupper(clockgenes)
  gene.list=rownames(E)



  E.full=E

  if(GTEx_names){
    gene.list=gsub("^.*_", "",gene.list)
    gene.list=gsub("\\|.*$","", gene.list)
  }
  clock.coord=match(clockgenes,gene.list)
  clock.coord=clock.coord[!is.na(clock.coord)]
  E=E[clock.coord,]
  geni=gene.list[clock.coord]

  if(is.matrix(E)==FALSE){return()}
  if(mean.centre.E==TRUE){E=sweep(E,1,rowMeans(E),FUN="-")}

  ###  initialization of parameters of the probabilistic model ###

  if(is.null(sigma2)){sigma2 = mean(apply(E, 1, var))}
  if(is.null(u)){u=0.2}
  if(is.null(tau2)){tau2=4/(24+(ncol(E)))}
  if(nrow(E)>500){tau2=tau2/10}
  phi=phi.start

  ### Phase initialization using spin glass model ###
  if(is.null(phi.start)){
    beta=1000
    J=J.tilde(E)
    Zeta = Zeta.mf.ordered(J, beta, ncol(E))
    phi<-Zeta[,2]+runif(ncol(E),-0.5,0.5)
  }


  if(pbar){pb = txtProgressBar(min = 0, max = iterations, style=3)}


  Ng=length(E[,1])
  N=length(E[1,])
  Ns=N
  sigma2.0=sigma2
  T=diag(3)*tau2
  T[1,1]=u^2
  dTinv=1/det(T)
  sigma2.m0=apply(E,1,var)
  sigma2.m1=sigma2.m0


  phi.0=phi

  c=cos(phi)
  s=sin(phi)
  one=rep(1,Ns)


  X=cbind(one,c,s)

  curvature=matrix(0, nrow = 2, ncol = Ns)

  S=t(E)%*%E/Ng

  S<-list()
  for(l in 1:Ng){
    S[[l]]=E[l,]%o%E[l,]
  }


  W<-rep(1,nrow(E))
  W.0<-W

  ### Start of EM procedure ###

  for(i in 1:iterations){
    phiold=phi
    Xold=X
    sigma2old=sigma2
    W.old=W



    Nn=(t(X)%*%X)
    Nninv<-solve(Nn)
    Tinv=solve(T)
    M=Nn+sigma2*Tinv
    Minv=solve(M)


    alpha=Minv%*%t(X)%*%t(E)
    alpha<-t(alpha)
    colnames(alpha)=c("mu","a","b")
    alphat=as.matrix(alpha)
    alpha=as.data.frame(alpha)
    if(TSM==FALSE){W<-W.0}
    Tot<-sum(W)

    A=sum(W*(alpha$a*alpha$a/sigma2+Minv[2,2]))/Tot
    B=sum(W*(alpha$a*alpha$b/sigma2+Minv[2,3]))/Tot
    C=sum(W*(alpha$b*alpha$a/sigma2+Minv[3,2]))/Tot
    D=sum(W*(alpha$b*alpha$b/sigma2+Minv[3,3]))/Tot


    K=matrix(c(A,B,C,D),nrow = 2,ncol = 2,byrow = T)
    if (any(is.nan(K))){return(list(alpha=alpha,weights=W,iteration=i))}

    al=apply(E,2, function(x){
      return(sum(W*(alpha$a*(x-alpha$mu)/sigma2-Minv[1,2]))/Tot)
    })
    be=apply(E,2, function(x){
      return(sum(W*(alpha$b*(x-alpha$mu)/sigma2-Minv[1,3]))/Tot)
    })

    O=as.matrix(rbind(al,be))
    if (any(is.nan(O))){return(list(alpha=alpha,weights=W,iteration=i))}

    ### Find numerically roots of the polinomial obtained from Lagrange multipliers ###

    rooted=apply(O,2,function(x){
      zero=B*B*C*C+A*A*D*D-x[1]*x[1]*(D*D+C*C)-x[2]*x[2]*(A*A+B*B)+2*x[1]*x[2]*(A*B+C*D)-2*A*B*C*D
      one=2*((A+D)*(A*D-B*C)-x[1]*x[1]*D-x[2]*x[2]*A+x[1]*x[2]*(B+C))
      two=A*A+D*D+4*A*D-x[1]*x[1]-x[2]*x[2]-2*B*C
      three=2*(A+D)
      four=1
      roots=polyroot(c(zero,one,two,three,four))
      return(roots)
    })

    ### Test each of the 2 or 4 roots to find the minimum ###

    ze=lapply(1:Ns,function(j){
      possible=lapply(rooted[,j],function(y){
        if(abs(Im(y)) > 1.0e-8){return(c(0,0,100000,100000))}
        K.lambda=K+diag(2)*Re(y)
        zet=solve(K.lambda,O[,j])

        phit=atan2(zet[2],zet[1])
        Xt=c(1,cos(phit),sin(phit))
        Mt=Xt%o%Xt

        Xt.old=c(1,cos(phiold[j]),sin(phiold[j]))
        Mt.old=Xt.old%o%Xt.old

        Qs=sapply(1:Ng,function(k){
          return(alphat[k,] %*% Mt %*% alphat[k,] - 2*alphat[k,]%*%Xt*E[k,j]+sigma2*sum(diag(Minv%*%Mt)))
        })

        Qs.old=sapply(1:Ng,function(k){
          return(alphat[k,] %*% Mt.old %*% alphat[k,] - 2*alphat[k,]%*%Xt.old*E[k,j]+sigma2*sum(diag(Minv%*%Mt.old)))
        })

        Q=sum(Qs*W/sigma2)/Tot
        Q.old=sum(Qs.old*W/sigma2)/Tot

        return(c(zet, Q, Q.old))
      })

      pox=matrix(unlist(possible),4,4,byrow = T)
      mpox=which.min(pox[,3])
      if(pox[mpox,3]==100000){
        cat("\n",pox,"not any solution on the circle at iteration",i,"\n")
        stop()
        return(list(phi=phiold,Qhist=Qhist,sigma=sigma2old, alpha=alpha, weights=W))
      }
      return(c(pox[mpox,],pox[,3]))
    })

    ### Save a lot of metrics for possible diagnostics ###

    ze=matrix(unlist(ze),8,N)
    rownames(ze)=c("cos","sin","Q", "Q.old", "Q1", "Q2", "Q3", "Q4")
    phi=atan2(ze[2,],ze[1,])%%(2*pi)

    if(i==1){
      Qhist=as.data.frame(t(ze))
      colnames(Qhist)=c("cos","sin","Q", "Q.old", "Q1", "Q2", "Q3", "Q4")
      Qhist$iteration=i
      Qhist$sample=c(1:Ns)
    }
    else{
      Qtemp=data.frame(t(ze))
      colnames(Qtemp)=c("cos","sin","Q", "Q.old", "Q1", "Q2", "Q3", "Q4")
      Qtemp$iteration=i
      Qtemp$sample=c(1:Ns)
      Qhist=rbind(Qhist,Qtemp)}
    P1.0=sapply(1:Ng, function(p){
      gamm=sqrt(dTinv*sigma2^3/det(M))
      exponent=E[p,]%*%X%*%Minv%*%t(X)%*%E[p,]/(2*sigma2)
      return(exp(exponent)*gamm)
    })
    W<-q*P1.0/(1-q+q*P1.0)
    W[is.nan(W)]<-1

    ### Exit the loop if convergence is reached ###

    if(max(abs(phi-phiold))<(0.001)){
      if(pbar){cat("\n algorithm has converged \n")}
      return(list(phi=phi,sigma=sigma2, alpha=alpha, weights=W, iter=i, sigma.m1=sigma2.m1, E=E.full, Qhist=Qhist, geni=geni, clock=E))
    }

    ### Update parameters for next step of EM ###

    c=cos(phi)
    s=sin(phi)
    one=rep(1,length(phi))
    X=cbind(one,c,s)

    Mold=Nn+sigma2*Tinv
    Moldinv=solve(Mold)

    sigma2.m1=sapply(1:Ng, function(s){
      return(sum(diag(S[[s]]-S[[s]]%*%Xold%*%Moldinv%*%t(X)))/Ns+0.01)
    })
    sigma2.m0=apply(E,1,var)
    sigma2<-mean(sigma2.m1*W+sigma2.m0*(1-W))
    sigma2.m1=sum(sigma2.m1*W)/sum(W)

    if(update.q){
      q=mean(W)
      if(q<0.05) q=0.05
      if(q>0.3) q=0.3}
    if(pbar){Sys.sleep(0.1)
      setTxtProgressBar(pb, i)}
  }
  if(pbar){close(pb)
    cat("\n")}
  ### Close function after specified iterattions ###
  return(list(phi=phi,Qhist=Qhist,sigma=sigma2, alpha=alpha, weights=W, iter=i,sigma.m1=sigma2.m1, E=E.full))
}

#' Calculates the circular correlation
#'
#' According to A Correlation Coefficient for Circular Data N. I. Fisher and A. J. Lee, https://doi.org/10.2307/2335547
#' The correlation is calculated either between two vectors or between all the columns of a given matrix. Note that cor.c(matrix) != cor.c(t(matrix)).
#'
#' @param x Required vector or matrix. Values MUST be in radiants. If x is a matrix do not provide the parameter y.
#' @param y Vector to be provided only if x is a vector too. Must have same length as x. Default NULL
#' @return The value of the circular correlation if x and y are vectors.
#' @return A symmetric matrix of dimensions ncol(x) times ncol(x) with entry i,j equal to the circualr correlation between x[,i] and x[,j]
#' @examples
#' cor.c(matrix)
#' cor.c(vector_x, vector_y)
#' @export
cor.c<- function(x,y=NULL){
  if(is.null(y)){return(cor.c.g(x))}
  if(length(x)!=length(y)){
    cat("the two vectors have different length")
    return(NULL)
  }
  n=length(x)
  U1=0
  U2=0
  U3=0
  for(i in 1:(n-1)){
    for(j in ((i+1):n)){
      U1=U1+sin(x[i]-x[j])*sin(y[i]-y[j])
      U2=U2+sin(x[i]-x[j])^2
      U3=U3+sin(y[i]-y[j])^2
    }
  }
  corr=U1/sqrt(U2*U3)
  return(corr)
}

cor.c.g<-function(x){
  off=ncol(x)
  cr=matrix(0, ncol=off, nrow=off)
  dimnames(cr)=list(colnames(x), colnames(x))
  for(l in 1:off){
    for(m in l:off){
      cr[l,m]=cr[m,l]=cor.c(x[,l], x[,m])
    }
  }
  return(cr)
}

#' Calculates the best possible shift and direction of two phase distributions to minimize the median absolute error (MAD).
#'
#' @param phi.0 Required vector of true phases.
#' @param phi Required vector of inferred phases. Values MUST be matched to phi.0.
#' @param period period of the periodic objects. Default 2*pi
#' @param mode Character to decide what the function returns/prints.
#' If "say" the MAD is printed at screen and not returned, if "return" MAD is only returned, if "no_median" it nor returned not printed, in all other cases it is both returned and printed. Default "forgotten"
#' @param N the number of division of the circle where to look for the shift, default 200
#' @param median_scale scale the median absolute error when saving/reporting it. Useful in the circadian application to have a median absolute error in hours with phases in radiants.
#' @return Either the inferred phases shifted and or flipped such that the MAD is minimized or a list where "phi" are the best phases and "median" is the value of the MAD
#' @examples
#' delta.phi(true.phi, inf.phi, mode="say")
#' @export

delta.phi<-function(phi.0, phi,period=2*pi, mode="forgotten", N=200, median_scale=1){
  phi.0=phi.0%%(period)
  phi=phi%%(period)
  isP=abs(period-2*pi)>0.001
  if(!isP){
    obj=calc.delta(phi.0, phi, N)
    bestphi=obj$phi
    mad=obj$median*median_scale
  }
  else{
    obj=calc.delta(phi.0/period*2*pi, phi/period*2*pi, N)
    bestphi=obj$phi*period/2/pi
    mad=obj$median*period/2/pi*median_scale
  }
  if(mode=="say"){
    cat("median:", mad, "\n")
    return(bestphi)}
  else if(mode=="return"){
    return(list(phi=bestphi, median=mad))
  }
  else if(mode=="no_median"){
    return(bestphi)
  }
  else{
    cat("median:", mad, "\n")
    return(list(phi=bestphi, median=(mad)))
  }
}



calc.delta<-function(phi.0, phi, N=200){
  mad<-12
  sdel=NULL
  for(j in 1:N){
    offset<-j/N*2*pi
    theta<-(phi-offset)%%(2*pi)
    del<-abs(theta-phi.0)%%(2*pi)
    delta<-del
    for (i in 1:length(phi)) {
      delta[i]=min(del[i], 2*pi-del[i])
    }
    if(median(delta)<mad){
      mad<-median(delta)
      bestphi<-theta
      j_temp=j
      sdel=delta}
    phi<-(-phi)%%(2*pi)
    theta<-(phi-offset)%%(2*pi)
    del<-abs(theta-phi.0)%%(2*pi)
    delta<-del
    for (i in 1:length(phi)) {
      delta[i]=min(del[i], 2*pi-del[i])
    }
    if(median(delta)<mad){
      mad<-median(delta)
      bestphi<-theta
      j_temp=-j
      sdel=delta}
  }
  return(list(phi=bestphi, median=(mad)))
}

#' Changes the inferred phases to obtain visually nicer plots using the fact that the phases live on a circular space where period-a=-a.
#'
#' @param realphi Required vector of true phases.
#' @param infphi Required vector of inferred phases.
#' @param period Period of the space where the phases live. Default 2*pi
#' @return Set of shifted inferred phases.
#' @examples
#' adjust.phases(true.phi, inf.phi)
#' @export

adjust.phases<-function(realphi, infphi, period=2*pi){
  if(abs(period-2*pi)>0.001){return(adjust.phases(realphi*2*pi/period, infphi*2*pi/period)*period/pi/2)}
  for(i in 1:length(realphi)){
    if(realphi[i]-infphi[i]>pi){infphi[i]=2*pi+infphi[i]}
    if(realphi[i]-infphi[i]<(-pi)){infphi[i]=infphi[i]-2*pi}
  }
  return(infphi)
}





### Spin glass approximation to have initial condition for EM ###

Zeta.mf.ordered<-function(J, beta, samples, A.0=0.1){
  iterations<-1000
  time_symmetry<-0
  A<-rep(A.0, samples)
  Theta<-runif(samples, 0, 2*pi)
  for (time in 1: iterations) {
    A.c=A*cos(Theta)
    A.s=A*sin(Theta)
    for (k in 1: samples) {
      u=beta*sum(A.c*J[k,])
      v=beta*sum(A.s*J[k,])

      modulo<-sqrt(u*u+v*v)
      Zeta_k<-c(u/modulo, v/modulo)
      Theta[k]<-atan2(Zeta_k[2],Zeta_k[1])
      A[k]<-(besselI(modulo,1)/(besselI(modulo,0)))
      if(is.na(modulo)){A[k] <- 1}
      else if(modulo>20){A[k] <- 1}
      if(is.nan(u) || is.nan(v)){Theta[k]<-runif(1,0,2*pi)}
      else if(u==0){
        if(v==0){
          Theta[k]<-runif(1,0,2*pi)
        }
      }
    }
  }
  for (k in 1:(samples-1)) {
    if (Theta[k] < Theta[k+1]){
      time_symmetry<-time_symmetry+1
    }
    if (Theta[k] > Theta[k+1]){
      time_symmetry<-time_symmetry-1
    }
  }

  time_symmetry<-sign(time_symmetry)
  if(time_symmetry==-1){
    for (k in 1:samples) {
      Theta[k]<-2*pi-Theta[k]
    }
  }
  Theta<-Theta%%(2*pi)
  Zeta<-cbind(A, Theta)
  return(Zeta)
}


### Calculation of spin glass interaction matrix ###

J.tilde<-function(E,n.genes=0,n.samples=0){
  sda<-as.matrix(E)
  if(n.genes==0){n.genes=nrow(E)}
  if(n.samples==0){n.samples=ncol(E)}
  Jtilde<-matrix(0,ncol = n.samples, nrow = n.samples)
  for(i in 1:n.samples){
    for (j in 1:n.samples) {
      Jtilde[i,j]<-sum(sda[,i]*sda[,j])/(n.samples*n.genes)
    }
  }
  diag(Jtilde)<-0
  rownames(Jtilde)<-colnames(E)
  colnames(Jtilde)<-colnames(E)
  return(Jtilde)
}


