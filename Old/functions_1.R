#things to be careful about:
#reference file for ordering

#These are the gense I use

library(ggplot2)
library(tidyverse)
#library(ggpubr)
library(doParallel)
library(ggrepel)
library(viridis)
library(tidyr)
library(data.table)
library(reshape2)# data table and this fight, maybe we can not load this
library(gplots)
library(ggforce)
library(preprocessCore)
library(gridExtra)


#If you want this function does the inference.
#it takes as arguments:
#The matrix you want to re-order
#The number of max Iterations
#It returns a list with:
#phi: a vector containing the phases
#sigma2: a number estimatig the probabilistically correct error of the data
#alpha: a tibble with the genes parameter
#pesi: a vector of with the probaility that each gene is indeed rhythmic
#other useless things

reconstruct.phi.from.ref<-function(E, fix.var=FALSE, gene.set=c("DBP","PER3", "NPAS2","ARNTL","NR1D1")){
  sx=x
  

  load("data/riferimento_mice_truephi.RData")

  mc=match(c("CLOCK","NFIL3"),names(full.ref))
  mc=mc[!is.na(mc)]
  ig=match(gene.set,names(full.ref))
  ig=ig[!is.na(ig)]
  if (length(ig)<2) {ig=match(c("DBP","PER3", "NPAS2","ARNTL","NR1D1"),names(full.ref))}
  full.ref=full.ref[ig]
  shifts=c(1:1000)*pi/500
  ts=complex(real=cos(shifts), imaginary=sin(shifts))
  
  
  se=E
  gene.list=gsub("^.*_", "",rownames(E))
  gene.list=gsub("\\|.*$","", gene.list)
  common=intersect(names(full.ref), gene.list)
  ref.com=full.ref[common]
  clock.coord=match(common,gene.list)
  clock.coord=clock.coord[!is.na(clock.coord)]
  E=E[clock.coord,]
  E=sweep(E,1,rowMeans(E),FUN="-")
  vars=apply(E,1,var)
  ampl=sqrt(2*vars)
  prop=ampl/Mod(ref.com)
  phi=c(1:ncol(E))
  if(fix.var==FALSE){
    pred.genes=ref.com*prop
    pred.mat=Re(pred.genes%o%Conj(ts))
    for(k in 1:ncol(E)){
      delta=colSums((pred.mat-E[,k])^2)
      phi[k]=which.min(delta)/500*pi
    }
  }
  if(fix.var==TRUE){
    a=min(prop)
    b=max(prop)
    px=a+c(0:1000)/1000*b
    px.mat=Mod(ref.com)%o%px
    dt=colSums((px.mat-vars)^2)
    lambda=px[which.min(dt)]
    pred.genes=ref.com*lambda
    pred.mat=Re(pred.genes%o%Conj(ts))
    for(k in 1:ncol(E)){
      delta=colSums(pred.mat-E[,k])
      phi[k]=which.min(delta^2)/500*pi
    }
  }
  return(list(E=se,phi=phi))
}




Phase_reconstruction<- function(E, iterations, clockgenes=NULL,tau2=NULL, u=NULL, sigma2=NULL, TSM=TRUE, mean.centre.E=TRUE, q=0.1, update.q=FALSE, wt=0, pbar=TRUE, pca=FALSE, pca_comp=30, normalize_small=FALSE,blood=FALSE, phi.start=NULL, rem.1.svd=FALSE, supervised.start=FALSE, standardize=FALSE, multi_tix=FALSE){
  require(MASS)
  id=as.vector(c(1,1,1,1))
  E=as.matrix(E)
  if(pca==TRUE){
    tempo=princomp(E)
    E=tempo$loadings
    if(ncol(E)>pca_comp){
      E=E[,1:pca_comp]
    }
    E=t(E)
  }
  
  
  
  
  
  cgenes.caps<-c("DBP"   ,  "PER3"  ,  "TEF"   ,  "NR1D2" ,  "PER1" ,   "PER2"  ,  "NPAS2" ,  "ARNTL" ,  "NR1D1"  , "TSC22D3", "LONRF3" , "FMO2"  ,  "CRY1"  ,  "CRY2","CIART" )
  TS_genes=c('DDIT4', 'GZMB', 'CAMKK1', 'GHRL', 'CLEC10 A', 'DTYMK', 'PER1', 'PDK1', 'NPEPL1', 'EPHX2', 'GPCPD1', 'MS4A3', 'GNG2', 'MUM1', 'IL13RA1', 'IL1B', 'STIP1', 'ID3', 'DHRS13', 'CHSY1', 'MEGF6', 'NR1D1', 'AK5', 'TCN1', 'ZNF438', 'CYB561', 'NSUN3', 'NR1D2', 'SLPI', 'POLH', 'CD38', 'PARP2', 'SYT11', 'TIAM2', 'PGPEP1', 'SH2D1B', 'CD1C', 'C12orf75', 'REM2', 'LLGL2', 'FKBP4')
  blood.clock.1 = c('PER1', 'PER2', 'CRY1', 'CRY2', 'NR1D2','LGALS3','HSPH1','FKBP4','ELMO2','CRISPLD2')
  blood.clock.2= c('CRISPLD2','TSC22D3','RBM3','PER1','DBP','IRAK3','CD99','FUS','DPYSL2','FAM129A','CX3CR1','UBE2J1','SERPINB9','AGFG1','C7orf50','STK17B','MLKL')
  blood.clock.3=unique(c(blood.clock.1,blood.clock.2))
  blood.clock.full=unique(c(blood.clock.3,TS_genes))
  mid.clock=unique(c(TS_genes,cgenes.caps))
  E=E[ , colSums(is.na(E)) == 0]
  if(standardize==TRUE){E=sweep(E,1,apply(E,1,sd),FUN="/")}
  E.full=E
  
  
  
  rownames(E)=toupper(rownames(E))
  if(is.null(clockgenes)){clockgenes=mid.clock}
  if(blood==TRUE){clockgenes=TS_genes}
  clockgenes=toupper(clockgenes)
  
  if(multi_tix==TRUE){
    tix=gsub("_.*$", "",rownames(E))
    gene.list=gsub("^.*_", "",rownames(E))
    gene.list=gsub("\\|.*$","", gene.list)
    gene.list=paste(tix,gene.list)
    clock.tix=NULL
    for(tx in unique(tix)){
      clock.tix=c(clock.tix, paste(tx, clockgenes))
    }
    clock.coord=match(clock.tix,gene.list)
    clock.coord=clock.coord[!is.na(clock.coord)]
    E=E[clock.coord,]
    geni=gene.list[clock.coord]
  }
  else{
    gene.list=gsub("^.*_", "",rownames(E))
    gene.list=gsub("\\|.*$","", gene.list)
    clock.coord=match(clockgenes,gene.list)
    clock.coord=clock.coord[!is.na(clock.coord)]
    E=E[clock.coord,]
    geni=gene.list[clock.coord]
  }
  if(normalize_small==TRUE){
    require(preprocessCore)
    TMP=normalize.quantiles(E)
    dimnames(TMP)=dimnames(E)
    E=TMP
  }
  SVD=svd(E)
  if(is.matrix(E)==FALSE){return()}
  if(mean.centre.E==TRUE){E=sweep(E,1,rowMeans(E),FUN="-")}
  if(rem.1.svd==TRUE){E=E-(SVD$u[,1]%o%SVD$v[,1])*SVD$d[1]}
  if(is.null(sigma2)){sigma2 = mean(apply(E, 1, var))}
  if(is.null(u)){u=0.2}
  if(is.null(tau2)){tau2=4/(24+(ncol(E)))}
  if(nrow(E)>500){tau2=tau2/10}
  phi=phi.start
  if(is.null(phi.start)){
    beta=1000
    J=J.tilde(E)
    Zeta = Zeta.mf.ordered(J, beta, ncol(E))
    phi<-Zeta[,2]+runif(ncol(E),-0.5,0.5)}
  if(supervised.start==TRUE){
    rec=reconstruct.phi.from.ref(E)
    phi=rec$phi
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
  
  # S is size Ns x Ns
  S=t(E)%*%E/Ng
  
  S<-list()
  for(l in 1:Ng){
    S[[l]]=E[l,]%o%E[l,]
  }
  
  
  W<-rep(1,nrow(E))
  W.0<-W
  for(i in 1:iterations){
    phiold=phi
    Xold=X
    sigma2old=sigma2
    W.old=W
    
    
    
    Nn=(t(X)%*%X)
    Nninv<-ginv(Nn)
    Tinv=ginv(T)
    M=Nn+sigma2*Tinv
    Minv=ginv(M)
    
    
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
    
    
    rooted=apply(O,2,function(x){
      zero=B*B*C*C+A*A*D*D-x[1]*x[1]*(D*D+C*C)-x[2]*x[2]*(A*A+B*B)+2*x[1]*x[2]*(A*B+C*D)-2*A*B*C*D
      one=2*((A+D)*(A*D-B*C)-x[1]*x[1]*D-x[2]*x[2]*A+x[1]*x[2]*(B+C))
      two=A*A+D*D+4*A*D-x[1]*x[1]-x[2]*x[2]-2*B*C
      three=2*(A+D)
      four=1
      roots=polyroot(c(zero,one,two,three,four))
      return(roots)
    })
    
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
        return(list(phi=phiold,Qhist=Qhist,sigma=sigma2old, alpha=alpha, pesi=W))
      }
      return(c(pox[mpox,],pox[,3]))
    })
    
    ze=matrix(unlist(ze),8,N)
    rownames(ze)=c("cos","sin","Q", "Q.old", "Q1", "Q2", "Q3", "Q4")
    phi=atan2(ze[2,],ze[1,])%%(2*pi)
    
    #    plot(phi.0, phi)
    
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
    
    if(max(abs(phi-phiold))<(0.001)){
      if(pbar){cat("\n algorithm has converged \n")}
      return(list(phi=phi,sigma=sigma2, alpha=alpha, pesi=W, iter=i, sigma.m1=sigma2.m1, E=E.full, Qhist=Qhist, geni=geni, clock=E))
    }
    c=cos(phi)
    s=sin(phi)
    one=rep(1,length(phi))
    X=cbind(one,c,s)
    
    Mold=Nn+sigma2*Tinv
    Moldinv=ginv(Mold)
    
    sigma2.m1=sapply(1:Ng, function(s){
      return(sum(diag(S[[s]]-S[[s]]%*%Xold%*%Moldinv%*%t(X)))/Ns+0.01)
    })
    sigma2.m0=apply(E,1,var) 
    sigma2<-mean(sigma2.m1*W+sigma2.m0*(1-W))
    sigma2.m1=sum(sigma2.m1*W)/sum(W)
    #if(all(W==0)){W=W+1}
    
    if(update.q){
      q=mean(W)
      if(q<0.05) q=0.05
      if(q>0.3) q=0.3}
    if(pbar){Sys.sleep(0.1)
      setTxtProgressBar(pb, i)}
  } 
  if(pbar){close(pb)
    cat("\n")}
  return(list(phi=phi,Qhist=Qhist,sigma=sigma2, alpha=alpha, pesi=W, iter=i,sigma.m1=sigma2.m1, E=E.full))
}


#This functions imports the data
#the are 3 arguments: 
#which clock genes to use
#if to mean center (it does so if Mean==FALSE)
#ir only import clockgenes (done to speed up the process)
#It returns a list with:
#E: matrix of all genes
#clock: matrix of clock genes
#N: number of genes
#M: number of samples


#these functions are not important to understand or use, they are just called by other functions

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

#input is a matrix with first column (row) a, second column (row) b and rownames (colnames) the genes
#a is cos coefficient, b is sin one
time.order<-function(A){
  
  symm=1
  phi_A=atan2(A["ARNTL","b"],A["ARNTL","a"])%%(2*pi)
  phi_P=atan2(A["PER3","b"],A["PER3","a"])%%(2*pi)
  dentro=0
  fuori=0
  for (gene in c("DBP", "NR1D1", "NPAS2")){
    phi_temp=atan2(A[gene,"b"],A[gene,"a"])%%(2*pi)
    # cat(phi_temp,"\n")
    if(phi_temp> min(phi_A,phi_P) && phi_temp< max(phi_A,phi_P)){
      dentro=dentro+1
    }
    if(phi_temp> max(phi_A,phi_P) || phi_temp< min(phi_A,phi_P)){
      fuori=fuori+1
    }
  }
  if(dentro>fuori){
    if(phi_P<phi_A){symm=-1}
  }
  else if(dentro<fuori){
    if(phi_P>phi_A){symm=-1}
  }
  A[,"b"]=A[,"b"]*symm
  return(A)
}


svd.from.out<-function(OUT, inter.genes, lock.gene="PER3"){
  gene.c=matrix(0, nrow = length(inter.genes),ncol=length(names(OUT)))
  dimnames(gene.c)=list(inter.genes,names(OUT))
  for(i in names(OUT)){ 
    out=OUT[[i]]
    out=order.out.setgene(out, gene=lock.gene)
    fit=out$data.fit
    inf.phi=out$phi
    exprx=OUT[[i]]$E
    gene.list=gsub("\\|.*$","",gsub("^.*_", "",fit$genes))
    clock.coord=sapply(inter.genes,function(x){match(x, gene.list)})
    clock.coord=clock.coord[!is.na(clock.coord)]
    fit=fit[clock.coord,]
    A=fit[,c("mean","amp", "a","b","R2", "genes")]
    rownames(A)=gsub("\\|.*$","",gsub("^.*_", "",A$genes))
    B=A
    A=complex(real = A[,"a"], imaginary =  A[,"b"])
    names(A)=rownames(B)
    #A=order.from.ref(A)
    #nr=rep(-Conj(A["PER3"])/sqrt(A["PER3"]*Conj(A["PER3"])),length(A))
    #A=A*nr
    gene.c[,i]=A[inter.genes]
  }
  gene.c[is.na(gene.c)]=0
  SVD=svd(gene.c)
  if(!is.null(lock.gene)&&!is.na(match(lock.gene,inter.genes))){
    i=lock.gene
    for(k in 1:length(SVD$d)){
      rot=Conj(SVD$u[which(rownames(gene.c)==i),k])/Mod(SVD$u[which(rownames(gene.c)==i),k])
      SVD$u[,k]=-SVD$u[,k]*rot*max(Mod(SVD$v[,k]))*SVD$d[k]
      SVD$v[,k]=-SVD$v[,k]*rot/max(Mod(SVD$v[,k]))
    }
  }
  rownames(SVD$u)=rownames(gene.c)
  rownames(SVD$v)=names(OUT)
  return(SVD)
}
cor.c<- function(x,y){
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


order.from.ref<-function(x, ref="mice truephi"){
  sx=x
  
  if(ref=="mice"){load("/data/shared/FelixProjects/XYmodels/comparisons/GTEx/riferimento_mice.RData")}
  else if(ref=="mice b"){load("/data/shared/FelixProjects/XYmodels/comparisons/GTEx/riferimento_mice_b.RData")}
  else if(ref=="gtex b"){load("/data/shared/FelixProjects/XYmodels/comparisons/GTEx/riferimento_gtex_b.RData")}
  else if(ref=="full"){load("/data/shared/FelixProjects/XYmodels/comparisons/GTEx/riferimento_TCGA_GTEX.RData")}
  else if(ref=="male"){load("/data/shared/FelixProjects/XYmodels/comparisons/GTEx/riferimento_gtex_male.RData")}
  else if(ref=="female"){load("/data/shared/FelixProjects/XYmodels/comparisons/GTEx/riferimento_gtex_female.RData")}
  else if(ref=="mice truephi"){load("/data/shared/FelixProjects/XYmodels/comparisons/GTEx/riferimento_mice_truuephi.RData")}
  else if(ref=="mice truephi old"){load("/data/shared/FelixProjects/XYmodels/comparisons/GTEx/riferimento_mice_truephi.RData")}
  else if(ref=="gtex old"){load("/data/shared/FelixProjects/XYmodels/comparisons/GTEx/riferimento_gtex.RData")}
  else {load("/data/shared/FelixProjects/XYmodels/comparisons/GTEx/riferimento_gtex_infphi.RData")}
  mc=match(c("CLOCK","NFIL3"),names(full.ref))
  mc=mc[!is.na(mc)]
  full.ref=full.ref[-mc]
  shifts=c(1:1000)*pi/500
  ts=complex(real=cos(shifts), imaginary=sin(shifts))
  common=intersect(names(full.ref), names(x))
  ref.mat=full.ref[common]%o%ts
  x=x[common]
  x=x/(sqrt(sum(Re(x*Conj(x)))))
  gen.p.scal=max(Re(t(ref.mat)%*%Conj(x)))#remember the inversion 
  inv.p.scal=max(Re(t(ref.mat)%*%x))
  if(gen.p.scal>inv.p.scal){return(sx)}
  return(Conj(sx))
}


order.out.setgene<-function(out,ref="mice truephi",gene="PER3"){
  out$has.been.flipped=0
  d.fit=out[["data.fit"]]
  A=complex(real=d.fit$a, imaginary = d.fit$b)
  names(A)=gsub("\\|.*$","",gsub("^.*_", "",d.fit$genes))
  AA=order.from.ref(A, ref=ref)
  if(sum(Im(A)*Im(AA))<0){
    out[["phi"]]=(-out[["phi"]])%%(2*pi)
    out$has.been.flipped=1
  }
  nz=complex(argument=(pi*(1-1/4)-Arg(AA["PER3"])))
  nr=rep(nz,length(AA))
  AA=AA*nr
  out[["phi"]]=(out[["phi"]]+Arg(nz))%%(2*pi)
  d.fit$b=Im(AA)
  d.fit$a=Re(AA)
  d.fit$phase=(Arg(AA)%%(2*pi))*12/pi
  out[["data.fit"]]=d.fit
  return(out)
}


infer_l=function(k,clockgenes=NULL){
  
  tix=k$tissue
  tp=k$type
  v=k$E
  out=Phase_reconstruction(v,500,clockgenes = clockgenes,standardize = TRUE) 
  out$tissue=tix
  out$type=tp
  return(out)
}


data.preparation.all.healthy<-function(clockgenes, Mean=FALSE, onlyclock=FALSE){
  exprs=vroom("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct")
  names=colnames(exprs)[-c(1:2)]
  samp <- fread('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')
  clockgenes=intersect(exprs$Description,clockgenes)
  if (onlyclock){
    exprs.clock.tibble=exprs[exprs$Description %in% clockgenes,]
    rn=exprs.clock.tibble$Description
    exprs.clock=as.matrix(exprs.clock.tibble[,-c(1,2)])
    rownames(exprs.clock)=rn
    log.mat=log2(exprs.clock+1)}
  else{
    rn=exprs$Description
    exprs.mat=as.matrix(exprs[,-c(1,2)])
    rownames(exprs.mat)=rn
    log.mat=log2(exprs.mat+1)
  }
  E.mat=list()
  for(tissue in unique(samp$SMTSD)){
    tiss=samp$SAMPID[samp$SMTSD==tissue]
    tiss.id=intersect(tiss, colnames(log.mat))
    mat=log.mat[,tiss.id]
    if(ncol(mat)>24){
      if(!Mean){
        mat<-mat-(replicate(ncol(mat), rowMeans(mat)))
      }
      clock=mat[clockgenes,]
      tx=paste(tissue, "GTEX", sep=" ")
      E.mat[[tx]]=list(E=mat, clock=clock, N=nrow(mat), M=ncol(mat), gene=rownames(mat))
    }
  }
  return(E.mat)
}