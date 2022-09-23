library(class)
library(scales)
library(ggplot2)
library(tidyverse)
#library(foreach)
#library(ggpubr)
library(doParallel)
library(parallel)
library(corrplot)
library(ggExtra)
#library(png)
library(ggrepel)
library(foreach)
library(tidyr)
library(vroom)
library(data.table)
library(reshape2)# data table and this fight, maybe we can not load this
library(gplots)
#library(ggforce)
library(ggstar)
library(lmtest)
library(reshape2)
library(grid)
library(gridExtra)
library(Chiral)


spliti= function(x,sp,nb){
  
  v=sapply(strsplit(x,sp),"[[",nb)
  return(v)
}

remove_covariates=function(x, samp.all){
  tt=x
  tt=subset(tt,rowMeans(tt)>0)
  tt=tt[,names(tt)%in%samp.all$SAMPID]
  tt.info=samp.all[names(tt),]
  tt.norm=tt
  
  if(ncol(tt)>10){
    for(j in 1:nrow(tt)){
      for.fit=data.frame(y=as.numeric(tt[j,]), isch=tt.info$SMSISH.f, age=tt.info$AGE, sex=as.factor(tt.info$SEX),death.ttype=as.factor(tt.info$DTHHRDY))
      if(length(unique(for.fit$sex))==1){
        resid=summary(lm(formula = y ~ isch + age  + death.ttype, data = for.fit))$residuals
      }else{
        resid=summary(lm(formula = y ~ isch + age + sex + death.ttype, data = for.fit))$residuals
      }
      tt.norm[j,]=NA
      tt.norm[j,as.numeric(names(resid))]=resid
    }
    
  }else{
    tt.norm= NULL
  }
  tt.norm
}

Norm.CPM<- function(CPM.all, high_filter=T, ncores=18){
  samp <- fread('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')
  samp.2 <- fread('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
  samp=samp[,c('SAMPID','SMRIN',"SMTSISCH","SMMAPRT","SMUNMPRT","SMMPPD","SMNABTCH","SMATSSCR","SMTSPAX","SMPTHNTS")]
  
  if (high_filter) samp=subset(samp,SMRIN > 6 & SMMPPD > 40000000 & SMMAPRT > 0.8 &  SMTSISCH > 0 & ((SMATSSCR == 0) | (SMATSSCR == 1) | (SMATSSCR == 2) | (SMATSSCR == 3) | (SMATSSCR == 4)  ))
  samp=subset(samp,SMRIN > 4)
  
  samp$subj.id=paste(spliti(samp$SAMPID,"-",1),spliti(samp$SAMPID,"-",2),sep="-")
  samp.all= data.frame(samp,samp.2[match(samp$subj.id,samp.2$SUBJID),])
  samp.all$SAMPID=gsub('-',"\\.",samp.all$SAMPID)
  rownames(samp.all) =samp.all$SAMPID
  samp.all$SMSISH.f= cut(samp.all$SMTSISCH, c(0,200,500,800,2000))
  
  CPM.all.norm=mclapply(CPM.all,remove_covariates,samp.all,mc.cores = ncores)
  return(CPM.all.norm)
}

CPM_to_E<-function(CPM.all.norm, min.samp=24, sep=NULL, samp=NULL){
  E.matrix=list()
  if(length(which(startsWith(names(CPM.all.norm), "Cells")))>0) CPM.all.norm=CPM.all.norm[-which(startsWith(names(CPM.all.norm), "Cells"))]
  if(length(which(startsWith(names(CPM.all.norm), "Whole")))>0) CPM.all.norm=CPM.all.norm[-which(startsWith(names(CPM.all.norm), "Whole"))]
  if(!is.null(samp)){
    if( tolower(sep)=="sex"){
      male=samp$sub.id[samp$SEX==1]
      female=samp$sub.id[samp$SEX==2]
      for(name in names(CPM.all.norm)){
        E=list()
        E$E=CPM.all.norm[[name]]
        E$E=E$E[,!is.na(E$E[1,])]
        if(!is.null(E$E)){
          E$tissue=name
          ml=fm=list()
          #gene="PER3"
          #E$E=as.matrix(subset(E$E,rowMeans(E$E)>0))
          #cat(name, " ", dim(E$E), gene, "is present", any(gsub("^.*_", "",rownames(E$E))==gene),  "\n")
          colnames(E$E)=spliti(colnames(E$E), "\\.", 2)
          ml$E=E$E[,colnames(E$E) %in% male]
          ml$type="GTEX-male"
          ml$tissue=name
          if(ncol(ml$E)>=min.samp){E.matrix[[paste(name, "Male", sep="-")]]=ml}
          fm$E=E$E[,colnames(E$E) %in% female]
          fm$type="GTEX-female"
          fm$tissue=name
          if(ncol(fm$E)>=min.samp){E.matrix[[paste(name, "Female", sep="-")]]=fm}
        }
      }
    }
    if( tolower(sep)=="age"){
      young=samp$sub.id[samp$AGE>60]
      old=samp$sub.id[samp$AGE<50]
      for(name in names(CPM.all.norm)){
        E=list()
        E$E=CPM.all.norm[[name]]
        E$E=E$E[,!is.na(E$E[1,])]
        if(!is.null(E$E)){
          E$tissue=name
          ml=fm=list()
          #gene="PER3"
          #E$E=as.matrix(subset(E$E,rowMeans(E$E)>0))
          #cat(name, " ", dim(E$E), gene, "is present", any(gsub("^.*_", "",rownames(E$E))==gene),  "\n")
          ml$E=E$E[,colnames(E$E) %in% young]
          ml$type="GTEX-old"
          ml$tissue=name
          if(ncol(ml$E)>=min.samp){E.matrix[[paste(name, "Young", sep="-")]]=ml}
          fm$E=E$E[,colnames(E$E) %in% old]
          fm$type="GTEX-young"
          fm$tissue=name
          if(ncol(fm$E)>=min.samp){E.matrix[[paste(name, "Old", sep="-")]]=fm}
        }
      }
    }
  }
  if(length(names(E.matrix))==0){
    for(name in names(CPM.all.norm)){
      E=list()
      E$E=CPM.all.norm[[name]]
      E$E=E$E[,!is.na(E$E[1,])]
      if(!is.null(E$E)){
        if(ncol(E$E)>min.samp){
          E$tissue=name
          gene="PER3"
          #E$E=as.matrix(subset(E$E,rowMeans(E$E)>0))
          #cat(name, " ", dim(E$E), gene, "is present", any(gsub("^.*_", "",rownames(E$E))==gene),  "\n")
          E.matrix[[name]]=E}
      }
    }
  }
  
  return(E.matrix)
}


infer_l=function(k,clockgenes=NULL){
  
  tix=k$tissue
  tp=k$type
  v=k$E
  out=CHIRAL(v,500,clockgenes = clockgenes, standardize = TRUE, GTEx_names=TRUE) 
  out$tissue=tix
  out$type=tp
  return(out)
}

f24_R2_cycling=function(x, t=2*(0:(length(x)-1)), period=24, offset=0, na.lm=5){

  kk = which(!is.na(x)==TRUE)
  x = x[kk]
  t = t[kk]
  n=length(x)
  #mu=mean(x)
  nb.timepoints=length(x)
  if(n<4)
  { 
    if(n==0) c(nb.timepoints=nb.timepoints, mean=NA, amp=NA, relamp=NA,phase=NA,pval=NA) 
    else 
    {
      c(nb.timepoints=nb.timepoints, mean=mean(x), amp=NA, relamp=NA,phase=NA,pval=NA)
    }
  }
  else
  {
    sig2=var(x)
    c=cos(2*pi*t/period)
    s=sin(2*pi*t/period)
    A = mean(x*c)-mean(x)*mean(c)
    B = mean(x*s)-mean(x)*mean(s)
    c1 = mean(c^2)-mean(c)^2
    c2 = mean(c*s)-mean(c)*mean(s)
    c3 = mean(s^2)-mean(s)^2
    b = (A*c2-B*c1)/(c2^2-c1*c3)
    a = (A-b*c2)/c1
    mu = mean(x)-a*mean(c)-b*mean(s)
    #	b=2*mean(x*s)
    x.hat=mu+a*c+b*s 
    sig2.1=var(x-x.hat)
    if(is.na(a)||is.na(b)) {c(nb.timepoints=nb.timepoints, mean=mean(x), amp=NA, relamp=NA,phase=NA,pval=NA)}
    else
    {
      p=3
      R2=0
      if(sig2>0) R2=1-sig2.1/sig2
      # http://www.combustion-modeling.com/downloads/beta-distribution-for-testing-r-squared.pdf
      # I checked that it works
      amp=max(x)-min(x)
      phase=period/(2*pi)*atan2(b, a)
      if(phase<0) phase=phase+period
      if(phase>period) phase=phase-period
      phase=(phase+offset)%%period
      pval = pbeta(R2, (p-1)/2, (n-p)/2, lower.tail = FALSE, log.p = FALSE)
      
      c(nb.timepoints=nb.timepoints, mean =mean(x), amp=2*sqrt(a^2+b^2),relamp=sqrt(a^2+b^2)/(mu),phase=phase, pval=pval,tot_err=sum((x-x.hat)^2),a=a,b=b, R2=R2, var=var(x))
    }
  }
}

harm_reg<-function(x, t, period){
  n=length(x)
  fit0=lm(x~1)
  c=cos(2*pi*t/period)
  s=sin(2*pi*t/period)
  fit1=lm(x~c+s)
  a=coef(fit1)[2]
  b=coef(fit1)[3]
  p.val=lrtest(fit1, fit0)$Pr[2]
  amp=2*sqrt(a^2+b^2)
  phase=atan2(b,a)%%(2*pi)
  phase=period*phase/(2*pi)
  return(c(pval=p.val,phase=phase,amp=amp,a=a,b=b, period=period, R2=summary(fit1)$r.squared))
}

Fit_OUT<-function(OUT,period=24, NA5=T){
  for(i in names(OUT)){ 
    E=OUT[[i]]$E
    de=sweep(E,1,rowMeans(E),FUN="-")
    if(NA5) E[abs(de)>5]=NA
    phase=OUT[[i]]$phi
    if(length(phase)<24){
      OUT[[i]]=NULL
      next
    }
    #dat.fit=as.data.frame(t(apply(E,1,f24_R2_cycling,t=24*as.numeric(phase)/(2*pi))))
    dat.fit=as.data.frame(t(apply(E, 1, harm_reg, 12*as.numeric(phase)/pi, period=period)))
    dat.fit$qval=p.adjust(dat.fit$pval, "BH")
    OUT[[i]]$data.fit=data.frame(dat.fit,E,genes=rownames(E))
    colnames(OUT[[i]]$data.fit)=c("pval","phase", "amp","mu", "a","b", "period", "R2","qval", colnames(E), "genes")
    OUT[[i]]$E=OUT[[i]]$data.fit[,colnames(E)]
  }
  return(OUT)
}


order.from.hc.ref<-function(x){
  sx=x
  #hard coded refernce
  full.ref=c(-0.43043911+0.02646768i, -0.21520829-0.08481562i, -0.18976262-0.06323798i, -0.26402067+0.03233651i, 
             -0.15460540-0.17813898i,  0.34503079+0.12988785i,  0.39656252+0.03842312i, -0.30343360+0.29074942i,  0.06226694-0.12542420i, -0.06159846-0.04409744i, -0.18018038-0.01939843i)
  names(full.ref)=c("DBP"   ,  "PER3"  ,  "TEF"   ,  "NR1D2"  ,   "PER2"  ,  "NPAS2" ,  "ARNTL" ,  "NR1D1"  , "CRY1"  ,  "CRY2",  "PER1")
  
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

order.out.setgene.hc<-function(out,gene="PER3"){
  out$has.been.flipped=0
  d.fit=out[["data.fit"]]
  A=complex(real=d.fit$a, imaginary = d.fit$b)
  names(A)=gsub("\\|.*$","",gsub("^.*_", "",d.fit$genes))
  AA=order.from.hc.ref(A)
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

Set_OUT<-function(OUT){
  for(i in names(OUT)){
    out=OUT[[i]]
    out=order.out.setgene.hc(out) #to ser per3 at 9
    sampz=gsub("\\..*$","", gsub("GTEX.","", colnames(out$E)))
    colnames(out$E)=sampz
    names(out$phi)=sampz
    OUT[[i]]=out
  }
  return(OUT)
}

Unique_donors<-function(OUT){
  all_people=NULL
  for(i in names(OUT)){
    sampz=colnames(OUT[[i]]$E)
    all_people=c(all_people, sampz)
  }
  people=unique(all_people)
  return(people)
}

Create_phi_matrix<-function(OUT, people){
  tix.study=names(OUT)
  phi_mat=matrix(nrow = length(people), ncol= length(names(OUT)))
  dimnames(phi_mat)=list(people, names(OUT))
  for(i in names(OUT)){
    out=OUT[[i]]
    sampz=colnames(out$E)
    phi_mat[sampz,i]=out$phi
  }
  phi_study=phi_mat[,tix.study]
  return(phi_study)
}

Phi.from.phi_mat<-function(phi.study, return.all=TRUE, min.tix=1, ct=1.95){
  phi.study=phi.study[rowSums(!is.na(phi.study))>=min.tix,]
  phi.comp=phi.study
  for(i in 1:ncol(phi.study)){
    phi.comp[,i]=complex(argument=phi.study[,i])
  }
  phit=rowMeans(phi.comp,na.rm = TRUE)
  phis=phit
  for(i in 1:nrow(phi.comp)){
    pp=phi.comp[i,!is.na(phi.comp[i,])]
    pc=Re(pp)
    if(length(pp)==0){
      phis=phis[-i]
    }
    else{
      for(k in 1:length(pp)){
        pc[k]=sum(Mod(pp+pp[k])>ct)
      }
      kstar=which.max(pc)
      gp=sum(pp[Mod(pp+pp[kstar])>ct])
      phis[i]=gp/Mod(gp)}
  }
  phit=phit[!is.na(phit)]
  score=Mod(phit)
  phit=phis
  if(return.all)  return(Arg(phit)%%(2*pi))
  sum(is.na(phit))
  good.s=names(phit[score>0.25])
  length(good.s)
  phi.fin=phi_mat[good.s,]
  samp.tix=nrow(phi.fin)-colSums(is.na(phi.fin))
  #hist(samp.tix)
  sum(samp.tix>100)
  phi=Arg(phit[good.s])%%(2*pi)
  return(phi)
}

Make_big_OUT<-function(E,phi){
  donors=names(phi)
  OUT=list()
  out=list()
  for (i in names(E)) {
    out=E[[i]]
    colnames(out$E)=gsub("\\..*$","", gsub("GTEX.","", colnames(out$E)))
    int_donors=intersect(donors, colnames(out$E))
    if(length(int_donors)>0){
      out$E=out$E[,int_donors]
      out$phi=phi[int_donors]
      OUT[[i]]=out
    }
  }
  return(OUT)
}

split_E_sex<-function(E, samp){
  males=samp$sub.id[samp$SEX==1]
  females=samp$sub.id[samp$SEX==2]
  A=list()
  for(i in names(E)){
    e=E[[i]]
    ee=e$E
    cn=colnames(ee)
    cn=gsub("\\..*$","", gsub("GTEX.","", cn))
    midx=match(males, cn)
    midx=midx[!is.na(midx)]
    if (length(midx)>24){
    e$E=ee[,midx]
    A[[paste(i, "Male", sep="-")]]=e
    }
    fidx=match(females, cn)
    fidx=fidx[!is.na(fidx)]
    if (length(fidx)>24){
    e$E=ee[,fidx]
    A[[paste(i, "Female", sep="-")]]=e
    }
  }
  return(A)
}

split_E_age<-function(E, samp){
  males=samp$sub.id[samp$AGE>60]
  females=samp$sub.id[samp$AGE<50]
  A=list()
  for(i in names(E)){
    e=E[[i]]
    ee=e$E
    cn=colnames(ee)
    cn=gsub("\\..*$","", gsub("GTEX.","", cn))
    midx=match(males, cn)
    midx=midx[!is.na(midx)]
    if (length(midx)>24){
      e$E=ee[,midx]
      A[[paste(i, "Male", sep="-")]]=e
    }
    fidx=match(females, cn)
    fidx=fidx[!is.na(fidx)]
    if (length(fidx)>24){
      e$E=ee[,fidx]
      A[[paste(i, "Female", sep="-")]]=e
    }  }
  return(A)
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

svd.from.out<-function(OUT, inter.genes, ENSG=F){
  gene.c=matrix(0, nrow = length(inter.genes),ncol=length(names(OUT)))
  dimnames(gene.c)=list(inter.genes,names(OUT))
  gtot=NULL
  for(i in names(OUT)){ 
    out=OUT[[i]]
    #out=order.out.setgene(out, gene=lock.gene)
    fit=out$data.fit
    gene.list=gsub("\\|.*$","",gsub("^.*_", "",fit$genes))
    if(ENSG){
      gls=gene.list
      gene.list=gsub("\\..*$","",gsub("_.*$", "",fit$genes))
    }
    clock.coord=sapply(inter.genes,function(x){match(x, gene.list)})
    clock.coord=clock.coord[!is.na(clock.coord)]
    fit=fit[clock.coord,]
    if(ENSG){
      gtot=union(gtot,fit$genes)
    }
    A=fit[,c("amp", "a","b", "genes")]
    rownames(A)=gsub("\\|.*$","",gsub("^.*_", "",A$genes))
    if(ENSG){
      rownames(A)=gsub("\\..*$","",gsub("_.*$", "",A$genes))
    }
    B=A
    A=complex(real = A[,"a"], imaginary =  A[,"b"])
    names(A)=rownames(B)
    #A=order.from.ref(A)
    #nr=rep(-Conj(A["PER3"])/sqrt(A["PER3"]*Conj(A["PER3"])),length(A))
    #A=A*nr
    gene.c[,i]=A[inter.genes]
  }
  gene.c[is.na(gene.c)]=0
  if(ncol(gene.c)*nrow(gene.c)==0){return(NULL)}
  SVD=svd(gene.c)
  # if(!is.null(lock.gene)&&!is.na(match(lock.gene,inter.genes))){
  #   i=lock.gene
  for(k in 1:length(SVD$d)){
    #rot=Conj(SVD$u[which(rownames(gene.c)==i),k])/Mod(SVD$u[which(rownames(gene.c)==i),k])
    mn=sum(SVD$v[,k])
    rot=Conj(mn)/Mod(mn)
    SVD$u[,k]=SVD$u[,k]*rot*max(Mod(SVD$v[,k]))*SVD$d[k]
    SVD$v[,k]=Conj(SVD$v[,k]*rot/max(Mod(SVD$v[,k])))
  }
  #}
  rownames(SVD$u)=rownames(gene.c)
  if(ENSG){
    gtot=unique(gtot)
    pox=match(rownames(SVD$u), gsub("\\..*$","",gsub("_.*$", "",gtot)))
    rownames(SVD$u)=gsub("\\|.*$","",gsub("^.*_", "",gtot))[pox]
  }
  rownames(SVD$v)=names(OUT)
  return(SVD)
}

svd.from.gene.big<-function(gene.big, inter.genez, ENSG=F, per3set=FALSE){
  genz=rownames(gene.big)
  if(ENSG){
    rownames(gene.big)=gsub("\\..*$","",gsub("_.*$", "",rownames(gene.big)))
  }
  if(!ENSG){rownames(gene.big)=gsub("\\|.*$","",gsub("^.*_", "",rownames(gene.big)))}
  inter.genes=intersect(inter.genez,rownames(gene.big))
  if(length(inter.genes)==0){return(NULL)}
  gene.c=gene.big[inter.genes,]
  gene.c[is.na(gene.c)]=0
  SVD=svd(gene.c)
  # if(!is.null(lock.gene)&&!is.na(match(lock.gene,inter.genes))){
  #   i=lock.gene
  for(k in 1:length(SVD$d)){
    #rot=Conj(SVD$u[which(rownames(gene.c)==i),k])/Mod(SVD$u[which(rownames(gene.c)==i),k])
    mn=sum(SVD$v[,k])
    rot=Conj(mn)/Mod(mn)
    if(per3set){
      if(any(rownames(gene.c)=="PER3")){
        rot=complex(argument=(-Arg(SVD$u[which(rownames(gene.c)=="PER3"),k])+9/12*pi))
      }
    }
    SVD$u[,k]=SVD$u[,k]*rot*max(Mod(SVD$v[,k]))*SVD$d[k]
    SVD$v[,k]=Conj(SVD$v[,k]*rot/max(Mod(SVD$v[,k])))
  }
  #}
  rownames(SVD$u)=rownames(gene.c)
  if(ENSG){
    gtot=genz
    pox=match(rownames(SVD$u), gsub("\\..*$","",gsub("_.*$", "",gtot)))
    rownames(SVD$u)=gsub("\\|.*$","",gsub("^.*_", "",gtot))[pox]
  }
  rownames(SVD$v)=colnames(gene.c)
  return(SVD)
}

create.gene.matrix<-function(OUT, norm=F){
  gene.traduction=tibble(ens=character(), symbol=character(), full=character())
  for(i in names(OUT)){ 
    out=OUT[[i]]
    #out=order.out.setgene(out, gene=lock.gene)
    fit=out$data.fit
    gene.tr.temp=tibble(ens=gsub("\\..*$","",gsub("_.*$", "",fit$genes)), symbol=gsub("\\|.*$","",gsub("^.*_", "",fit$genes)),full=fit$genes)
    gene.traduction=rbind(gene.traduction,gene.tr.temp)
  }
  inter.genes=as.character(unique(gene.traduction$full))
  gene.c=matrix(0, nrow = length(inter.genes),ncol=length(names(OUT)))
  dimnames(gene.c)=list(inter.genes,names(OUT))
  for(i in names(OUT)){ 
    out=OUT[[i]]
    fit=out$data.fit
    dats=out$E
    gene.list=fit$genes
    clock.coord=sapply(inter.genes,function(x){match(x, gene.list)})
    clock.coord=clock.coord[!is.na(clock.coord)]
    fit=fit[clock.coord,]
    fit=fit[!is.na(fit[,"amp"]),]
    #dats=as.matrix(dats[clock.coord,])
    #vars=apply(dats, 1, var)
    if(nrow(fit)>0){
      if(norm){A=complex(modulus = fit[,"R2"], argument = (fit[,"phase"]*pi/12))}
      else{A=complex(real = fit[,"a"], imaginary =  fit[,"b"])}
      names(A)=fit$genes
      gene.c[,i]=A[inter.genes]}
  }
  gene.c[is.na(gene.c)]=0
  return(gene.c)
}

Plot_density<-function(OUT, phi, R.plot=FALSE, R.df=FALSE,cut=0.1, varz="pval", comp="small", title_param="", compet="", sz=20, th=1){
  phi.df=tibble(phi=rep((1:10000-1)/5000*pi,1))
  phi.df$hour=phi.df$phi/pi*12
  phi.df$count=phi.df$phi
  kapp=20
  for(s in 1:length(phi.df$phi)){
    #if(s<length(phi.df$phi)/2)
    phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(phi))))
  }
  phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2] #normalize in the linear sense
  phi.df$dens=phi.df$count/sqrt(sum(phi.df$count^2)*phi.df$phi[2]/2)#to romalize the radar plot integral to one
  #ggplot(phi.df)+geom_line(aes(x=hour,y=dens), colour="darkblue")+coord_polar()+scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ylim(0,NA)+theme_minimal()+labs(x="time of day", y="density", title="Donor circadian phase")
  phi.df$kind="donors"
  s.phi.df=phi.df
  pho=NULL
  for(name in names(OUT)){
    out=OUT[[name]]
    if(comp=="big") {phi=out$data.fit[which(out$data.fit[,varz]>cut), "phase"]*pi/12
    compet="bigger than"
    }
    else if(comp=="small") {phi=out$data.fit[which(out$data.fit[,varz]<cut), "phase"]*pi/12
    compet="smaller than"
    }
    else { cat("comp variable can be set ti either big or small")
      stop()
    }
    pho=c(pho,phi)
  }
  for(s in 1:length(phi.df$phi)){
    #if(s<length(phi.df$phi)/2)
    phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(pho))))
  }
  phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2] #normalize in the linear sense
  phi.df$dens=phi.df$count/sqrt(sum(phi.df$count^2)*phi.df$phi[2]/2)#to romalize the radar plot integral to one
  phi.df$kind="genes"
  phid=rbind(phi.df, s.phi.df)
  plt=ggplot(phid)+geom_line(aes(x=hour,y=dens, color=kind),size=th)+
    scale_color_manual(values= c("#008ecc", "#111e6c"))+coord_polar()+scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ylim(0,NA)+
    theme_minimal()+labs(x="time of day", y="density", title=paste("Density" , title_param))+
    theme(legend.position = NULL, text = element_text(size=sz))
  print(plt)
  if(R.plot & R.df) return(list(plot=plt, df=phi.df))
  else if(R.plot) return(plt)
  else if(R.df) return(phi.df)
}

PhiOUT_from_CPM<-function(CPM){
  E=CPM_to_E(CPM)
  if(length(which(startsWith(names(E), "Cells")))>0) OUT=OUT[-which(startsWith(names(E), "Cells"))]
  if(length(which(startsWith(names(E), "Whole")))>0) OUT=OUT[-which(startsWith(names(E), "Whole"))]
  
  gene_inf<-c("DBP"   ,  "PER3"  ,  "TEF"   ,  "NR1D2" ,  "PER1" ,   "PER2"  ,  "NPAS2" ,  "ARNTL" ,  "NR1D1",  "CRY1"  ,  "CRY2","CIART" )
  
  OUT= mclapply(E,infer_l, gene_inf, mc.cores=16)
  
  OUT=Fit_OUT(OUT)
  
  OUT=Set_OUT(OUT)
  
  donors=Unique_donors(OUT)
  
  phi_matrix=Create_phi_matrix(OUT, donors)
  
  phi=Phi.from.phi_mat(phi_matrix)
  
  return(list(phi=phi, OUT=OUT))
}

Tissue_density<-function(OUT, full_col, loc=NULL, cut=0.1, varz="pval"){
  full_col=full_col[full_col$Class!="Cells",]
  dec_names=full_col$`Short name`
  names(dec_names)=full_col$`Full name`
  colorandum=full_col$`# color`
  names(colorandum)=full_col$`Full name`
  nmz=unique(full_col$Class)
  if(is.null(loc)){loc=paste(getwd(),"/",sep="")}
  loco=paste(loc,"Tissue_density_",gsub("\\.", "-", cut),"_", varz, ".pdf", sep="")
  pdf(loco)
  phi.full=NULL
  pho=NULL
  for(name in names(OUT)){
    out=OUT[[name]]
    phi=out$data.fit[which(out$data.fit[,varz]<cut), "phase"]*pi/12
    pho=c(pho,phi)
    phi.df=tibble(phi=rep((1:1000-1)/500*pi,1))
    phi.df$hour=phi.df$phi/pi*12
    phi.df$count=phi.df$phi
    kapp=20
    for(s in 1:length(phi.df$phi)){
      #if(s<length(phi.df$phi)/2)
      phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(phi))))
    }
    phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2]
    ppt=ggplot(phi.df)+geom_line(aes(x=hour,y=dens), colour=colorandum[name])+coord_polar()+scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ylim(0,NA)+theme_minimal()+labs(x="time of day", y="density", title=name)
    print(ppt)
    phi.df$tissue=name
    phi.full=rbind(phi.full, phi.df)
  }
  ppt=ggplot(phi.full)+geom_line(aes(x=hour,y=dens,colour=tissue))+coord_polar()+scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ylim(0,NA)+theme_minimal()+labs(x="time of day", y="density", title="All Tissues")+scale_color_manual(values=colorandum[names(OUT)])+theme(legend.position = "none")
  print(ppt)
  dev.off()
}

Plot_profiles<-function(OUT, tissues, geni, plt.hist=FALSE, R.plot=FALSE, val="pval", text_size=12, title_add="", always_plot_fit=F, vcut=0.2, xl=NULL, period=24){
  p2=period/2
  nmz=names(OUT)
  pltlist=list()
  l=1
  tissuex= intersect(nmz, tissues)
  if(length(tissuex)==0){
    tissuex=NULL
    for (nm in tissues) {
      ttix=names(OUT)[startsWith(names(OUT),nm)]
      tissuex=c(tissuex,ttix)
    }
  }
  for(i in tissuex){ 
    sz=7
    ang=90
    hjst=1
    out=OUT[[i]]
    fit=out$data.fit
    fut=fit
    inf.phi=out$phi
    exprx=out$E
    gene.list=gsub("\\|.*$","",gsub("^.*_", "",fit$genes))
    full=fit[,c("amp","a","b","R2", "genes")]
    full$genz=gsub("\\|.*$","",gsub("^.*_", "",full$genes))
    full$R=sqrt(full$a^2+full$b^2)
    GE=fit[,c("a","b", "genes")]
    GL=GE$genes
    #mus=fit$mean
    GE=complex(real = GE[,"a"], imaginary =  GE[,"b"])
    names(GE)=gsub("\\|.*$","",gsub("^.*_", "",GL))
    #names(mus)=names(GE)
    if(is.null(xl)) xl="Donor internal phase [h]"
    for(k in 1:length(geni)){
      fg=geni[k]
      idx=match(fg, gsub("\\|.*$","",gsub("^.*_", "",rownames(out$E))))
      if(is.na(idx)){
        cat("gene ", fg, " not present in ", i, "\n")
        stop()
      }
      ts=(1:1000)*pi/500
      pred=GE[fg]*Conj(complex(argument = (ts*24/period)))#+mus[fg]
      fit=tibble(gene=Re(pred), phase=ts,type="fit", time=ts*12/pi, sd=0, w=0)#, conf=1)
      geneplt=tibble(gene=as.matrix(out$E)[idx,],phase=out$phi,type="data", time=out$phi*12/pi, sd=0, w=0)#, conf=out$score)
      aver=sdev=1:12
      phibinz=c(0:11)*pi/6
      for (bin in 1:(length(phibinz))) {
        phizb=which(out$phi>phibinz[bin])
        phizs=which(out$phi<phibinz[bin]+pi/6)
        phiz=intersect(phizb, phizs)
        aver[bin]=mean(as.matrix(out$E)[idx,phiz], na.rm = T)
        sdev[bin]=sd(as.matrix(out$E)[idx,phiz], na.rm=T)/sqrt(length(phiz))
      }
      tibaver=tibble(gene=aver, phase=(phibinz+pi/12), type="avg", time=(phibinz+pi/24)*12/pi, sd=sdev, w= 0.2)
      if (always_plot_fit) geneplt=rbind(geneplt, tibaver, fit)
      else{
        if(fut[[val]][idx]<vcut) geneplt=rbind(geneplt, tibaver, fit)
        else geneplt=rbind(geneplt, tibaver)
      }
      if(! (val%in%colnames(fut))){val="pval"}
      if(k==1) {
        phi.tb=tibble(phi=inf.phi*12/pi)
        brk=seq(0, 24, 1)
        
        if(!is.character(title_add)) ttl=NULL
        else ttl=paste(i, title_add)
        
        if(plt.hist) {pltlist[[l]]=ggplot(phi.tb, aes(x=phi))+geom_histogram(breaks=brk,fill="#008ECC")+
          theme_minimal()+ lims(x=c(0,24))+labs(x=xl,y="Frequency", title=ttl) +
          theme(text=element_text(size=text_size))
        l=l+1
        pltlist[[l]]=ggplot(geneplt)+geom_point(aes(x=time, y=gene, color=type, size=type,alpha=type))+
          labs(x=xl,y="Centred expression [log2]", subtitle=paste(fg, val, signif(fut[[val]][idx],2)))+theme_minimal()+ lims(x=c(0,24))+#c(-5,5))+
          theme(legend.position = "none",text=element_text(size=text_size))+scale_size_manual(values=c(3,1,0.00007))+scale_alpha_manual(values=c(1,0.8,1))+scale_color_manual(values=c("#008ECC","#B0DFE5","#111E6C","#0B49BD"))+geom_errorbar(aes(x=time, y=gene, ymin=gene-sd, ymax=gene+sd,colour="errorbar", width=w))
        }
        else{ pltlist[[l]]=ggplot(geneplt)+geom_point(aes(x=time, y=gene, color=type, size=type,alpha=type))+
          labs(x=xl,y="Centred expression [log2]", title=ttl, subtitle=paste(fg, val, signif(fut[[val]][idx],2)))+
          theme_minimal()+ lims(x=c(0,24))+#c(-5,5))+
          theme(legend.position = "none",text=element_text(size=text_size))+scale_size_manual(values=c(3,1,0.00007))+scale_alpha_manual(values=c(1,0.8,1))+scale_color_manual(values=c("#008ECC","#B0DFE5","#111E6C","#0B49BD"))+geom_errorbar(aes(x=time, y=gene, ymin=gene-sd, ymax=gene+sd,colour="errorbar", width=w))
        }
        l=l+1
      }
      else {pltlist[[l]]=ggplot(geneplt)+geom_point(aes(x=time, y=gene, color=type, size=type))+
        labs(x=xl ,y="Centred expression [log2]", subtitle=paste(fg, val, signif(fut[[val]][idx], 2)))+theme_minimal()+ lims(x=c(0,24))+#c(-5,5))+
        theme(legend.position = "none", text=element_text(size=text_size))+scale_size_manual(values=c(3, 1, 0.00007))+scale_color_manual(values=c("#008ECC","#B0DFE5","#111E6C","#0B49BD"))+geom_errorbar(aes(x=time, y=gene, ymin=gene-sd, ymax=gene+sd, color="errorbar",width=w))
      l=l+1}
    }
  }
  #if(is.null(top_title)) 
  if(length(pltlist)==1){f.plot=pltlist[[1]]}
  else f.plot=do.call(grid.arrange,c(grobs=pltlist, as.table=FALSE))
  #if(!is.null(top_title))  f.plot=do.call(grid.arrange,c(grobs=pltlist, as.table=FALSE))
  print(f.plot)
  if(R.plot) return(f.plot)
}

Plot_cSVD<-function(input, genes, full_col=NULL,loc=NULL, ENSG=F, ncomp=NULL, CT=15, ymax=1.5, dot_size=1.5, label_size=3.5, text_size=12, sectors=1, max_ov=Inf){
  if(!is.null(full_col)){
  full_col=full_col[full_col$Class!="Cells",]
  dec_names=full_col$`Short name`
  names(dec_names)=full_col$`Full name`
  colorandum=full_col$`# color`
  names(colorandum)=full_col$`Full name`}
  if(is.null(full_col)){
    nmtz=names(input)
    if(!is.list(input)) nmtz=colnames(input)
    colorandum=nmtz
    names(colorandum)=nmtz
    colorandum[nmtz]="#008ECC"
  }
  if(is.null(loc)){loc=paste(getwd(),"/",sep="")}
  loco= ifelse(endsWith(loc, ".pdf"), loc, paste(loc,"_cSVD", ".pdf", sep=""))
  pdf(loco)
  if (is.list(input))SVD=svd.from.out(input,genes,ENSG=ENSG)
  else SVD=svd.from.gene.big(input ,genes, ENSG =ENSG)
  
  if(is.null(ncomp)) ncomp=ncol(SVD$u)
  for(i in 1:ncomp){
    #i=1
    labs=rownames(SVD$u)
    for (j in 1:sectors) {
      in_sector=Arg(SVD$u[,i])%%(2*pi)>(j-1)*(2*pi/sectors) &  Arg(SVD$u[,i])%%(2*pi)<j*(2*pi/sectors)
      tempu=Mod(SVD$u[in_sector,i])
    if(length(tempu)>CT){
      cuts=tempu[base::order(-tempu)][CT]
      #labs[which(!(labs %in% gene_inf))]=""
      smz=Mod(SVD$u[,i])<cuts
      labs[as.logical(smz*in_sector)]=""
    }
    }
    colz=rownames(SVD$v)
    shepz=rep("all",nrow(SVD$v))
    deco_names=dec_names
    if(any(endsWith(rownames(SVD$v),"Male")) || any(endsWith(rownames(SVD$v),"young"))){
      colz=gsub("-old", "",gsub("-young", "",gsub("-Female", "",gsub("-Male", "", rownames(SVD$v)))))
      shepz=gsub("^.*-", "", rownames(SVD$v))
      dec1_names=paste(dec_names, unique(shepz)[1], sep="-")
      dec2_names=paste(dec_names,unique(shepz)[2], sep="-")
      names(dec1_names)=paste(names(dec_names), unique(shepz)[1], sep="-")
      names(dec2_names)=paste(names(dec_names), unique(shepz)[2], sep="-")
      deco_names=c(dec1_names,dec2_names)
    }
    varexp=round((SVD$d[i]^2)/sum(SVD$d^2)*100,0)
    clock.spread=round(1-Mod(sum(SVD$u[,i]^2/Mod(SVD$u[,i])))/sum(Mod(SVD$u[,i])),4)
    print(qplot(Arg(SVD$u[,i])%%(2*pi)*12/pi,Mod(SVD$u[,i]))+geom_point(size=dot_size/2)+
            geom_label_repel(aes(label=labs),size = label_size,max.overlaps = max_ov)+coord_polar()+ylim(0,ymax)+ theme_minimal()+
            scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ 
            labs(title=paste("Genes SVD, mode", i,"E.V.",varexp, "%"),x =paste("Argument of component", i), y = paste("Modulus of component", i))+
            theme(text=element_text(size=text_size)))
    
    print(qplot(Arg(SVD$v[,i])%%(2*pi)*12/pi,Mod(SVD$v[,i]),colour=colz, shape=shepz)+geom_point(size=dot_size/2)+ 
            theme_minimal()+scale_color_manual(values=colorandum[colz])+coord_polar()+aes(ymin=0, xmin=0,xmax=2*pi)+
            geom_label_repel(aes(label=deco_names[rownames(SVD$v)]),size = label_size, max.overlaps = max_ov)+theme(legend.position = "none")+
            scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+scale_y_continuous(breaks=seq(0, 1, by=0.25),expand=c(0,0), lim=c(0, 1.1))+ 
            labs(title=paste("Tissues SVD, mode", i, "E.V.",varexp, "%"),x =paste("Argument of component", i), y = paste("Modulus of component", i))+
            theme(text=element_text(size=text_size)))
    
    print(qplot(Arg(SVD$u[,i])%%(2*pi)*12/pi,Mod(SVD$u[,i]))+geom_point(size=dot_size)+coord_polar()+
            aes(ymin=0, xmin=0,xmax=2*pi)+ylim(0,ymax)+theme_minimal()+
            scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ 
            labs(title=paste("Tissues SVD, mode", i,"E.V.",varexp, "%"),x =paste("Argument of component", i), y = paste("Modulus of component", i))+
            theme(text=element_text(size=text_size)))
    
    print(qplot(Arg(SVD$v[,i])%%(2*pi)*12/pi,Mod(SVD$v[,i]),colour=colz, shape=shepz)+geom_point(size=dot_size)+ 
            theme_minimal()+scale_color_manual(values=colorandum[colz])+coord_polar()+aes(ymin=0, xmin=0,xmax=2*pi)+
            theme(legend.position = "none")+scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+
            scale_y_continuous(breaks=seq(0, 1, by=0.25),expand=c(0,0), lim=c(0, 1.1))+ 
            labs(title=paste("Tissues SVD, mode", i, "E.V.",varexp, "%"),x =paste("Argument of component", i), y = paste("Modulus of component", i))+
            theme(text=element_text(size=text_size)))
  }
  dev.off()
}

Plot_cumulatives<-function(OUT,full_col, pcut=1, qcut=0.2, Rcut=0.5, loc=NULL, vals=c("R", "pval","qval"), sz=20, th=1){
  full_col=full_col[full_col$Class!="Cells",]
  dec_names=full_col$`Short name`
  names(dec_names)=full_col$`Full name`
  colorandum=full_col$`# color`
  names(colorandum)=full_col$`Full name`
  nmz=unique(full_col$Class)
  if(is.null(loc)) loc=paste(getwd(), "/", sep="")
  for(val in vals){
    all=NULL
    tbt=NULL
    #val="pval"
    
    pdf(paste(loc, val,"_R", gsub("\\.", "",toString(Rcut)),"_p", gsub("\\.", "",toString(pcut)),"_q", gsub("\\.", "",toString(qcut)),"_all.pdf", sep=""))
    for(j in names(OUT)){
      out=OUT[[j]]
      fit=out$data.fit
      #interesting.genes=inter.genes
      inf.phi=out$phi
      exprx=out$E
      gene.list=gsub("\\|.*$","",gsub("^.*_", "",fit$genes))
      clock.coord=sapply(out$geni,function(x){match(x, gene.list)})
      clock.coord=clock.coord[!is.na(clock.coord)]
      full=fit[,c("amp","a","b","R2", "genes","qval", "pval")]
      full=subset(full, pval<pcut)
      full=subset(full, qval<qcut)
      full$R=2*sqrt(full$a^2+full$b^2)
      full=subset(full, R>Rcut)
      if(nrow(full)>1){
        full$kind=names(OUT)[j]
        all=rbind(all,full)
        breaks = seq(0, max(full[,val]), by=0.01) 
        if(val!="R") breaks=breaks^7
        amp.cut = cut(full[,val], breaks, right=FALSE) 
        freq.cut = table(amp.cut) 
        if(val=="R") freq.cut = rev(freq.cut) #also tbz=tibble(R=rev(breaks[-1]), n.genes=cumsum.frq, kind=names(OUT)[j])
        cumsum.frq=c(cumsum(freq.cut),nrow(full[,val]))+1
        tbz=tibble(R=breaks[-1], n.genes=cumsum.frq, kind=j)
        if(val=="R") tbz=tibble(R=rev(breaks[-1]), n.genes=cumsum.frq, kind=j)
        tbt=rbind(tbt,tbz)
      }
    }
    if(!is.null(tbt)){
      if (val=="R") print(ggplot(tbt)+geom_line(aes(x=R,y=n.genes, color=kind), size=th)+scale_y_log10(limits=c(1,2500))+
                            labs(title="All tissues",y="# of genes", x="log2(peak to trough)", color="Category")+theme_minimal()+
                            scale_color_manual(values=colorandum[unique(tbt$kind)])+theme(legend.position = "none",text = element_text(size=sz)))
      else print(ggplot(tbt)+geom_line(aes(x=R,y=n.genes, color=kind), size=th)+scale_y_log10()+scale_x_continuous(trans=reverselog_trans(10))+
                   labs(title="All tissues",y="# of genes", x=val)+theme_minimal()+scale_color_manual(values=colorandum[unique(tbt$kind)])+
                   theme(legend.position = "none",text = element_text(size=sz)))
    }
    for(i in nmz){ 
      ct=which(full_col$Class==i)
      tbt=NULL
      for(j in full_col$`Full name`[ct]){
        out=OUT[[j]]
        fit=out$data.fit
        #interesting.genes=inter.genes
        inf.phi=out$phi
        exprx=out$E
        gene.list=gsub("\\|.*$","",gsub("^.*_", "",fit$genes))
        clock.coord=sapply(out$geni,function(x){match(x, gene.list)})
        clock.coord=clock.coord[!is.na(clock.coord)]
        full=fit[,c("amp","a","b","R2", "genes","qval", "pval")]
        full=subset(full, pval<pcut)
        full=subset(full, qval<qcut)
        full$R=2*sqrt(full$a^2+full$b^2)
        full=subset(full, R>Rcut)
        if(nrow(full)>1){
          full$kind=names(OUT)[j]
          all=rbind(all,full)
          breaks = seq(0, max(full[,val]), by=0.01) 
          if(val!="R") breaks=breaks^7
          amp.cut = cut(full[,val], breaks, right=FALSE) 
          freq.cut = table(amp.cut) 
          if(val=="R") freq.cut = rev(freq.cut) #also tbz=tibble(R=rev(breaks[-1]), n.genes=cumsum.frq, kind=names(OUT)[j])
          cumsum.frq=c(cumsum(freq.cut),nrow(full[,val]))+1
          tbz=tibble(R=breaks[-1], n.genes=cumsum.frq, kind=j)
          if(val=="R") tbz=tibble(R=rev(breaks[-1]), n.genes=cumsum.frq, kind=j)
          tbt=rbind(tbt,tbz)
        }
      }
      if(!is.null(tbt)){
        if(val=="R") print(ggplot(tbt)+geom_line(aes(x=R,y=n.genes, color=kind, linetype=kind), size=th)+
                             scale_y_log10(limits=c(1,2500))+labs(title=i,y="# of genes", x="log2(peak to trough)",  color="Category", linetype="Category")+
                             theme_minimal()+scale_color_manual(values=colorandum[unique(tbt$kind)])+
                             theme(text = element_text(size=sz), legend.position = "none"))#+theme(legend.position = "none"))
        else print(ggplot(tbt)+geom_line(aes(x=R,y=n.genes, color=kind, linetype=kind), size=th)+
                     scale_y_log10()+scale_x_continuous(trans=reverselog_trans(10))+labs(title=i,y="# of genes", x=val)+
                     theme_minimal()+scale_color_manual(values=colorandum[unique(tbt$kind)])+theme(text = element_text(size=sz), legend.position = "top"))#+theme(legend.position = "none"))
      }
    }
    # tbt=NULL
    # all$div=gsub("^.*-","",  all$kind)
    # for(dv in unique(all$div) ){
    #   full=filter(all, div==dv)
    #   breaks = seq(0, max(full[,val]), by=0.01) 
    #   amp.cut = cut(full[,val], breaks, right=FALSE) 
    #   freq.cut = table(amp.cut) 
    #   if(val=="R"){freq.cut = rev(freq.cut)} #also tbz=tibble(R=rev(breaks[-1]), n.genes=cumsum.frq, kind=names(OUT)[j])
    #   cumsum.frq=c(cumsum(freq.cut),nrow(full[,val]))+1
    #   tbz=tibble(R=breaks[-1], n.genes=cumsum.frq, kind=names(OUT)[j],div=dv)
    #   if(val=="R"){tbz=tibble(R=rev(breaks[-1]), n.genes=cumsum.frq, kind=names(OUT)[j], div=dv)}
    #   tbt=rbind(tbt,tbz)
    # }
    #print(ggplot(tbt)+geom_line(aes(x=R,y=n.genes, color=div))+labs(title="all",y="# of genes", x=val)+theme_minimal()+scale_color_manual(values=colorandum[unique(tbt$div)])+theme(legend.position = "none"))#+scale_color_manual(values=c("green3","blue","red","darkred","lightblue","purple")))
    dev.off()
  }
  
  
}


Plot_circular_density<-function(phi, kapp=20, x.t="time of day", y.t="density", title.t="Donor circadian phase", norm.area=FALSE, h24=FALSE, colorandum=NULL, colz="Category", legend.pos="right", th=1, sz=20){
  if(length(phi[[1]])==1){
    phi.df=tibble(phi=rep((1:10000-1)/5000*pi,1))
    phi.df$hour=phi.df$phi/pi*12
    phi.df$count=phi.df$phi
    if(h24) phi=phi/12*pi
    for(s in 1:length(phi.df$phi)){
      #if(s<length(phi.df$phi)/2)
      phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(phi))))
    }
    phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2] #normalize in the linear sense
    if(norm.area) phi.df$dens=phi.df$count/sqrt(sum(phi.df$count^2)*phi.df$phi[2]/2)#to romalize the radar plot integral to one
    print(ggplot(phi.df)+geom_line(aes(x=hour,y=dens), colour="darkblue",size=th)+coord_polar()+scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ylim(0,NA)+theme_minimal()+labs(x=x.t, y=y.t, title=title.t))
  }
  else{
    sphi=phi
    phi.full=NULL
    for(name in names(sphi)){
      phi=sphi[[name]]
      phi.df=tibble(phi=rep((1:1000-1)/500*pi,1))
      phi.df$hour=phi.df$phi/pi*12
      phi.df$count=phi.df$phi
      kapp=20
      for(s in 1:length(phi.df$phi)){
        #if(s<length(phi.df$phi)/2)
        phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(phi))))
      }
      phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2]
      phi.df$div=name
      phi.full=rbind(phi.full, phi.df)
    }
    if(!is.null(colorandum)){
      ppt=ggplot(phi.full)+geom_line(aes(x=hour,y=dens,colour=div),size=th)+coord_polar()+scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ylim(0,NA)+theme_minimal()+
        labs(x=x.t, y=y.t, title=title.t, colour=colz)+scale_color_manual(values=colorandum[unique(phi.full$div)])+theme(text= element_text(size=sz), legend.position = legend.pos)
    }
    else{
      ppt=ggplot(phi.full)+geom_line(aes(x=hour,y=dens,colour=div),size=th)+coord_polar()+scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ylim(0,NA)+theme_minimal()+labs(x=x.t, y=y.t, title=title.t, colour=colz)+scale_color_brewer(palette = "Set1")
    }
    print(ppt)
  }
}


manh.dist<-function(mat, by.col=FALSE, norm=FALSE){
  if(by.col==TRUE){mat=t(mat)}
  rz=nrow(mat)
  mat=sign(mat)
  dist.mat=matrix(0,nrow =rz, ncol = rz)
  for(i in 1:rz){
    for(j in (i):rz){
      delt=mat[i,]-mat[j,]
      dist.mat[j,i]=sum(abs(delt))
      if(norm){
        dist.mat[j,i]=dist.mat[j,i]/(sum(abs(mat[i,]))+sum(abs(mat[i,])))
      }
    }
  }
  dist.mot=as.dist(dist.mat)
  #dimnames(dist.mat)=list(rownames(mat), rownames(mat))
  return(dist.mot)
}



CM_from_matrix<-function(mat, col, x.t=NULL, y.t=NULL, title.t=NULL, sz=5, ang=90,hjst=1, plt.size=T, ordering=NULL, leg.pos="none"){
  meltmat=melt(mat)
  colnames(meltmat)=c("gene1", "gene2", "correlation")
  a=hclust(manh.dist(mat))
  colord=colnames(mat)[a$order]
  if(!is.null(ordering)) colord=ordering
  meltmat$gene1=factor(meltmat$gene1, colord)
  meltmat$gene2=factor(meltmat$gene2, colord)
  if(!plt.size){
    print(ggplot(meltmat, aes(y = gene1,
                              x = gene2)) +        ## global aes
            #geom_tile(aes(fill = phase)) +         ## to get the rect filled
            geom_tile(alpha=0)+
            geom_point(aes(colour = correlation), size=5)  +    ## geom_point for circle illusion
            #geom_point(aes(colour = correlation,
            #              size =abs(correlation)))  +    ## geom_point for circle illusion
            #scale_color_gradient(low = "red",high = "violet")+       ## color of the corresponding aes
            scale_colour_gradientn(colours= col, limits=c(-1,1), guide="none")+
            #scale_size(range = c(0, 5), guide="none")+             ## to tune the size of circles
            theme_minimal()+
            theme(axis.text.x = element_text(angle = ang, hjust=hjst),
                  axis.text.y = element_text(),text = element_text(size=sz*1.5),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())+
            labs(x=x.t, y=y.t, title= title.t))
  }
  if(plt.size){
    print(ggplot(meltmat, aes(y = gene1,
                              x = gene2)) +        ## global aes
            #geom_tile(aes(fill = phase)) +         ## to get the rect filled
            geom_tile(alpha=0)+
            geom_point(aes(colour = correlation,
                           size =abs(correlation)))  +    ## geom_point for circle illusion
            #scale_color_gradient(low = "red",high = "violet")+       ## color of the corresponding aes
            scale_colour_gradientn(colours= col, limits=c(-1,1))+
            scale_size(range = c(0, 7.5), guide="none")+             ## to tune the size of circles
            theme_minimal()+
            theme(axis.text.x = element_text(angle = ang, hjust=hjst),
                  axis.text.y = element_text(),text = element_text(size=sz*1.5),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())+
            labs(x=x.t, y=y.t, title= title.t))
  }
}


corr.genes<- function(E, genes, ENSG=FALSE, sampz=NULL, dot_size=1.5){
  gl=rownames(E)
  gene.list=gsub("\\|.*$","",gsub("^.*_", "",gl))
  if(ENSG){
    gls=gene.list
    gene.list=gsub("\\..*$","",gsub("_.*$", "",gl))
  }
  if(is.null(sampz)) sampz=colnames(E)
  idx=match(genes, gene.list)
  idx=idx[!is.na(idx)]
  E=as.matrix(E)
  E=E[,intersect(sampz,colnames(E))]
  E=sweep(E,1,rowMeans(E),FUN="-")
  Ex=E[idx,]
  corrm= cor(t(Ex))
  dimnames(corrm)=list(gene.list[idx],gene.list[idx])
  return(corrm)
}

Plot_radar<-function(OUT, genes, tissue="unspecified tissue",ymax=1.5, dot_size=1.5, label_size=5, text_size=12, WP=""){
  if(tissue %in% names(OUT)) tissues=tissue
  else if(any(startsWith(names(OUT), tissue))) tissues=names(OUT)[which(startsWith(names(OUT), tissue))]
  for (tissue in tissues){
    out=OUT[[tissue]]
    fit=out$data.fit
    inf.phi=out$phi
    exprx=out$E
    gene.list=gsub("\\|.*$","",gsub("^.*_", "",fit$genes))
    clock.coord=sapply(out$geni,function(x){match(x, gene.list)})
    clock.coord=clock.coord[!is.na(clock.coord)]
    geni=gene.list[clock.coord]
    full=fit[,c("amp","a","b","R2", "genes","qval", "pval", "phase")]
    full=full[clock.coord,]
    full$R=2*sqrt(full$a^2+full$b^2)
    full$geni=geni
    
    print(ggplot(full, aes(x=phase, y=R))+geom_point(size=dot_size)+
            geom_label_repel(aes(label=geni),size = label_size)+coord_polar()+ylim(0,ymax)+ theme_minimal()+
            scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ 
            labs(title=paste(WP, "in", tissue),x =paste("Peak phase"), y = paste("Amplitude"))+
            theme(text=element_text(size=text_size)))
  }
  
}


Multitix_radar<-function(OUT, genes, colorandum, tissues=NULL, qcut=0.2){
  if(!is.null(tissues))OUT=OUT[tissues]
  ptt=NULL
  gene.r=matrix(NA,ncol=length(OUT), nrow = length(genes))
  dimnames(gene.r)=list(genes,names(OUT))
  
  for (i in names(OUT)){
    out=OUT[[i]]
    df=out$data.fit
    df=subset(df, qval<qcut)
    gt=complex(modulus=df$amp, argument=df$phase/12*pi)
    names(gt)=gsub("\\|.*$","",gsub("^.*_", "",df$genes))
    com=intersect(names(gt),genes)
    gene.r[com,i]=gt[com]
  }
  geni=genes
  rg=Arg(gene.r[,1])%%(2*pi)
  geni=geni[base::order(rg)]
  ptt=NULL
  for(g in geni){
    G=gene.r[g,]
    G=G[!is.na(G)]
    G=G[G!=0]
    if(any(endsWith(names(OUT),"Male")) || any(endsWith(names(OUT),"young"))){
      colz=gsub("-old", "",gsub("-young", "",gsub("-Female", "",gsub("-Male", "", names(G)))))
      shepz=gsub("^.*-", "", names(G))
      
    }
    else{
      colz=names(G)
      shepz="all"
    }
        pt=tibble(tix=names(G), dist=Mod(G), ang=(Arg(G)%%(2*pi)*12/pi), gene=g, sex=gsub("^.*-","",names(G)),col=colorandum[colz],Category=shepz)
    ptt=rbind(ptt, pt)
  }
  ptt$gene_g=factor(ptt$gene, levels = geni)
  {
    if(length(unique(ptt$Category))==1){
    g1=ggplot(subset(ptt,gene%in%geni[1:4])) +geom_point(aes(x=ang, y=dist, color=col, shape=Category))+scale_color_identity()+coord_polar(start=0, direction=1)+labs(x="Peak phase", y="Amplitude")+
      scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+facet_grid(~gene_g)+ylim(0,NA)+ theme_minimal()+ theme(legend.position = "none")
    g2=ggplot(subset(ptt,gene%in%geni[5:8]))+geom_point(aes(x=ang, y=dist, color=col, shape=Category))+scale_color_identity() +coord_polar(start=0, direction=1)+labs(x="Peak phase", y="Amplitude")+
      scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+facet_grid(~gene_g)+ ylim(0,NA)+ theme_minimal()+theme(legend.position = "none")
    g3=ggplot(subset(ptt,gene%in%geni[9:12]))+geom_point(aes(x=ang, y=dist, color=col, shape=Category)) +scale_color_identity()+coord_polar(start=0, direction=1)+labs(x="Peak phase", y="Amplitude")+
      scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+facet_grid(~gene_g)+ ylim(0,NA)+ theme_minimal()+theme(legend.position = "none")
    }
    else{
    g1=ggplot(subset(ptt,gene%in%geni[1:4])) +geom_point(aes(x=ang, y=dist, color=col, shape=Category))+scale_color_identity()+coord_polar(start=0, direction=1)+labs(x="Peak phase", y="Amplitude")+
      scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+facet_grid(~gene_g)+ylim(0,NA)+ theme_minimal()
    g2=ggplot(subset(ptt,gene%in%geni[5:8]))+geom_point(aes(x=ang, y=dist, color=col, shape=Category))+scale_color_identity() +coord_polar(start=0, direction=1)+labs(x="Peak phase", y="Amplitude")+
      scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+facet_grid(~gene_g)+ ylim(0,NA)+ theme_minimal()
    g3=ggplot(subset(ptt,gene%in%geni[9:12]))+geom_point(aes(x=ang, y=dist, color=col, shape=Category)) +scale_color_identity()+coord_polar(start=0, direction=1)+labs(x="Peak phase", y="Amplitude")+
      scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+facet_grid(~gene_g)+ ylim(0,NA)+ theme_minimal()
    }
    
    grid.arrange(g1,g2,g3)
  }
  
}
cplx.dist<-function(mat, by.col=FALSE, norm=FALSE){
  if(by.col==TRUE){mat=t(mat)}
  rz=nrow(mat)
  dist.mat=matrix(0,nrow =rz, ncol = rz)
  for(i in 1:rz){
    for(j in (i):rz){
      delt=mat[i,]-mat[j,]
      dist.mat[j,i]=sqrt(sum(Mod(delt)^2))
      if(norm){
        dist.mat[j,i]=dist.mat[j,i]/(sqrt(sum(Mod(mat[i,])^2))+sqrt(sum(Mod(mat[i,])^2)))
      }
    }
  }
  dist.mat=as.dist(dist.mat)
  #dimnames(dist.mat)=list(rownames(mat), rownames(mat))
  return(dist.mat)
}

manh.dist<-function(mat, by.col=FALSE, norm=FALSE){
  if(by.col==TRUE){mat=t(mat)}
  rz=nrow(mat)
  mat=sign(mat)
  dist.mat=matrix(0,nrow =rz, ncol = rz)
  for(i in 1:rz){
    for(j in (i):rz){
      delt=mat[i,]-mat[j,]
      dist.mat[j,i]=sum(abs(delt))
      if(norm){
        dist.mat[j,i]=dist.mat[j,i]/(sum(abs(mat[i,]))+sum(abs(mat[i,])))
      }
    }
  }
  dist.mot=as.dist(dist.mat)
  #dimnames(dist.mat)=list(rownames(mat), rownames(mat))
  return(dist.mot)
}


Heatmap_WP<-function(int.WP=NULL, div="all", MS=F, pcut=0.01, min_genes=3, save_plot=F, path=NULL, plt_name=NULL, text_size=10, text_angle=90, shape_size=c(0,5), mask=NULL){
  colroma=vroom("/scratch/For_cedric_with_pizza/data/romaO.txt", col_names = FALSE, delim=" ")
  colroma$hex=paste("#", as.hexmode(round(colroma$X1*255,0)),as.hexmode(round(colroma$X2*255,0)),as.hexmode(round(colroma$X3*255,0)), sep="")
  WP_full=vroom("/data/shared/FelixProjects/XYmodels/Data/WikiPathway_2021_Human", col_names = FALSE, delim=" ")
  SS.age=get(load("/scratch/For_cedric_with_pizza/27082020/Aging/rhythmic_table_withthresh_NA5.RData"))
  SS.MF=get(load("/scratch/For_cedric_with_pizza/27082020/MF/rhythmic_table_withthresh_NA5.RData"))
  if(save_plot){
    if(is.null(path)){
      if(is.null(plt_name)){
        if(MS) pdf(paste("/data/shared/FelixProjects/XYmodels/new_plots_meanphi/heatmaps/WP_MS_", div, ".pdf", sep=""))
        else pdf(paste("/data/shared/FelixProjects/XYmodels/new_plots_meanphi/heatmaps/WP_", div, ".pdf", sep=""))
      }
      else{
        pdf(paste("/data/shared/FelixProjects/XYmodels/new_plots_meanphi/heatmaps/", plt_name, ".pdf", sep=""))
      }
    }
    else{
      pdf(path)
    }
  }
  
  
  load(paste("/scratch/For_cedric_with_pizza/OUT_files/OUT_big_mean_final_", div, "_musc_no_cells.RData", sep=""))
  
  if(div=="MF") SS=SS.MF
  if(div=="age_n5") SS=SS.age
  
  
  WP=WP_full$X1
  for(i in 1:nrow(WP_full)){
    tmp=WP_full[i,]
    tmp=tmp[which(!is.na(tmp))]
    tmp.str=tmp[1][1]
    for (j in 2:ncol(tmp)) {
      tmp.str=paste(tmp.str, tmp[j][1])
    }
    WP[i]=tmp.str
  }
  
  WP.ls=lapply(WP, function(x){
    a=strsplit(x, "\t")[[1]]
    return(a[3:length(a)])
  })
  
  names(WP.ls)=unlist(lapply(WP, function(x){a=strsplit(x, "\t")[[1]][1]}))
  
  if(is.null(int.WP)) int.WP=sample(names(WP.ls), 20)
  MM=F
  if(!is.null(mask)){
    if(is.null(nrow(mask))){
      if (!MS) {
        unspaced=gsub(" ", "", names(OUT))
        if(length(intersect(mask, unspaced))>length(intersect(mask, names(OUT)))){
          mask=names(OUT)[match(mask,unspaced)]
        }
        OUT=OUT[intersect(mask, names(OUT))]
      }
      else {
        unspaced=gsub(" ", "", names(SS))
        if(length(intersect(mask, unspaced))>length(intersect(mask, names(SS)))){
          mask=names(SS)[match(mask,unspaced)]
        }
        SS=SS[intersect(mask, names(SS))]
      }
    }
    else{
      MM=T
      if (!MS) {
        unspaced=gsub(" ", "", names(OUT))
        if(length(intersect(colnames(mask), unspaced))>length(intersect(colnames(mask), names(OUT)))){
          interm=intersect(colnames(mask), unspaced)
          mask=mask[,interm]
          colnames(mask)=names(OUT)[match(colnames(mask),unspaced)]
        }
        OUT=OUT[intersect(colnames(mask), names(OUT))]
      }
      else {
        unspaced=gsub(" ", "", names(SS))
        if(length(intersect(colnames(mask), unspaced))>length(intersect(colnames(mask), names(SS)))){
          interm=intersect(colnames(mask), unspaced)
          mask=mask[,interm]
          colnames(mask)=names(SS)[match(colnames(mask),unspaced)]
        }
        SS=SS[intersect(colnames(mask), names(SS))]
      }
      int.WP=intersect(int.WP, rownames(mask))
    }
  }
  
  
  
  if(MS){ 
    mat_order=matrix(0, length(int.WP), length(names(SS)))
    dimnames(mat_order)=list(int.WP,names(SS))
  }
  else{
    mat_order=matrix(0, length(int.WP), length(names(OUT)))
    dimnames(mat_order)=list(int.WP,names(OUT))
  }
  
  uterms=NULL
  for(n in int.WP){
    
    if(! n %in% names(WP.ls)){
      idx=which(grepl(n, names(WP.ls), fised=T))[1]
      n=names(WP.ls)[idx]
      if(is.null(n)){next}
    }
    
    genes=WP.ls[[n]]
    
    if(div=="all" | !MS){
      for (ou in names(OUT)) {
        if(MM){
          if(mask[n,ou]==0 | !mask[n,ou]){next}
        }
        out=OUT[[ou]]
        fit=out$data.fit
        
        gene.list=gsub("\\|.*$","",gsub("^.*_", "",fit$genes))
        clock.coord=sapply(genes,function(x){match(x, gene.list)})
        clock.coord=clock.coord[!is.na(clock.coord)]
        full=fit[,c("a","b","R2", "genes", "pval")]
        full=full[clock.coord,]
        full=subset(full, pval<pcut)
        gene.lista=gsub("\\|.*$","",gsub("^.*_", "",full$genes))
        if(nrow(full)> min_genes){
          full$R=2*sqrt(full$a^2+full$b^2)
          full$cpx=complex(argument = atan2(full$b,full$a))
          scpx=sum(full$cpx)
          tmp.tb=tibble(term=n,tissue=ou, phase=Arg(scpx), mean_ampl=mean(full$R), mean_R2=mean(full$R2), genes=paste(gene.lista, collapse=","))
          mat_order[n,ou]=scpx
        }
        else tmp.tb=tibble(term=n,tissue=ou, phase=0, mean_ampl=0, mean_R2=0, genes=character())
        uterms=rbind(uterms, tmp.tb)
      }
    }
    else{
      for (su in names(SS)) {
        if(MM){
          if(mask[n,su]==0 | !mask[n,su]){next}
        }
        ss=SS[[su]]
        idxs=which(grepl(su, names(OUT), fixed=T))
        if(length(idxs)!=2){next}
        for (j in idxs) {
          
          nm=names(OUT)[j]
          md=l=gsub("^.*-", "", names(OUT)[j])
          if(md %in% c("old", "Female")){mds=c(3,4,5)}
          else{mds=c(2,4,5)}
          
          out=OUT[[nm]]
          
          st=subset(ss, model %in% mds)
          
          fit=out$data.fit
          fit=fit[intersect(rownames(fit), rownames(st)),]
          
          gene.list=gsub("\\|.*$","",gsub("^.*_", "",fit$genes))
          clock.coord=sapply(genes,function(x){match(x, gene.list)})
          clock.coord=clock.coord[!is.na(clock.coord)]
          full=fit[,c("a","b","R2", "genes", "pval")]
          full=full[clock.coord,]
          full=subset(full, pval<pcut)
          gene.lista=gsub("\\|.*$","",gsub("^.*_", "",full$genes))
          if(nrow(full)> min_genes){
            full$R=2*sqrt(full$a^2+full$b^2)
            full$cpx=complex(argument = atan2(full$b,full$a))
            scpx=sum(full$cpx)
            tmp.tb=tibble(term=n,tissue=su, phase=Arg(scpx), mean_ampl=mean(full$R), mean_R2=mean(full$R2),Category= md, genes=paste(gene.lista, collapse=","))
            mat_order[n,su]=mat_order[n,su]+scpx
          }
          else tmp.tb=tibble(term=n,tissue=su, phase=0, mean_ampl=0, mean_R2=0, Category= character(), genes=character())
          uterms=rbind(uterms, tmp.tb)
        }
      }
    }
  }
  
  if(is.null(uterms)){
    cat("No pathways passed cutoffs \n")
    return(NA)
  }
  if(nrow(uterms)==0){
    cat("No pathways passed cutoffs \n")
    return(NA)
  }
  
  sz=text_size
  ang=text_angle
  hjst=1
  uterms$phase=uterms$phase%% (2*pi) /pi*12
  a=hclust(d=cplx.dist(t(mat_order), norm=FALSE))
  colord=colnames(mat_order)[a$order]
  a=hclust(d=cplx.dist((mat_order), norm=FALSE))
  roword=rownames(mat_order)[a$order]
  uterms$term=factor(uterms$term, roword)
  uterms$tissue=factor(uterms$tissue, colord)
  
  
  if(div=="all" | !MS){
    print(ggplot(uterms, aes(y = term,
                             x = tissue)) +        ## global aes
            #facet_grid(Generic~., scales = "free", space = "free")+
            #geom_tile(aes(fill = phase)) +         ## to get the rect filled
            geom_tile(alpha=0)+
            geom_point(aes(colour = phase,
                           size =(mean_R2)))  +    ## geom_point for circle illusion
            #scale_color_gradient(low = "red",high = "violet")+       ## color of the corresponding aes
            scale_colour_gradientn(colours= colroma$hex)+
            scale_size(range = shape_size)+             ## to tune the size of circles
            theme_minimal()+
            theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=sz*1.5), strip.text.y = element_text(angle = 0)) +
            labs(colour="Mean peak phase", size="Mean R2")) 
    
    print(ggplot(uterms, aes(y = term,
                             x = tissue)) +        ## global aes
            #facet_grid(Generic~., scales = "free", space = "free")+
            #geom_tile(aes(fill = phase)) +         ## to get the rect filled
            geom_tile(alpha=0)+
            geom_point(aes(colour = phase,
                           size =(mean_ampl)))  +    ## geom_point for circle illusion
            #scale_color_gradient(low = "red",high = "violet")+       ## color of the corresponding aes
            scale_colour_gradientn(colours= colroma$hex)+
            scale_size(range = shape_size)+             ## to tune the size of circles
            theme_minimal()+
            theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=sz*1.5), strip.text.y = element_text(angle = 0)) +
            labs(colour="Mean peak phase", size="Mean amplitude")) 
  }
  else{
    print(ggplot(uterms, aes(y = term,
                             x = tissue)) +        ## global aes
            #facet_grid(Generic~., scales = "free", space = "free")+
            #geom_tile(aes(fill = phase)) +         ## to get the rect filled
            geom_tile(alpha=0)+
            geom_point(aes(colour = phase,
                           size =(mean_R2),
                           shape=Category))  +    ## geom_point for circle illusion
            #scale_color_gradient(low = "red",high = "violet")+       ## color of the corresponding aes
            scale_colour_gradientn(colours= colroma$hex)+
            scale_size(range = shape_size)+             ## to tune the size of circles
            scale_shape_manual(values=c("\u25E3","\u25E5"))+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=sz), strip.text.y = element_text(angle = 0)) +
            labs(colour="Mean peak phase", size= "Mean R2")) 
    
    print(ggplot(uterms, aes(y = term,
                             x = tissue)) +        ## global aes
            #facet_grid(Generic~., scales = "free", space = "free")+
            #geom_tile(aes(fill = phase)) +         ## to get the rect filled
            geom_tile(alpha=0)+
            geom_point(aes(colour = phase,
                           size =(mean_ampl),
                           shape=Category))  +    ## geom_point for circle illusion
            #scale_color_gradient(low = "red",high = "violet")+       ## color of the corresponding aes
            scale_colour_gradientn(colours= colroma$hex)+
            scale_size(range = shape_size)+ ## to tune the size of circles
            scale_shape_manual(values=c("\u25E3","\u25E5"))+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=sz), strip.text.y = element_text(angle = 0)) +
            labs(colour="Mean peak phase", size="Mean amplitude")) 
    
  }
  if(save_plot) dev.off()
}


outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

which.max.matrix <- function(z) {
  if (!is.matrix(z)) {
    stop("Not a matrix")
  }
  m <- nrow(z)
  n <- ncol(z)
  # take care of NAs
  ind <- which.max(z)
  iy <- trunc((ind - 1)/m) + 1
  ix <- ind - (iy - 1) * m
  return(cbind(ix, iy))
}


delta.phi<-function(phi.0, phi,mode="forgotten"){
  mad<-12
  sdel=NULL
  for(j in 1:length(phi)){
    offset<-phi[j]-phi.0[j]
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
    offset<-phi[j]-phi.0[j]
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
  if(mode=="say"){
    cat("median:", mad*12/pi, "\n")
    return(bestphi)}
  else if(mode=="return"){
    return(list(phi=bestphi, median=mad*12/pi))
  }
  else if(mode=="no_median"){
    return(bestphi)
  }
  else if(mode=="deltas"){
    return(list(phi=bestphi, median=(mad*12/pi), deltas=sdel))
  }
  else{
    cat("median:", mad*12/pi, "\n")
    return(list(phi=bestphi, median=(mad*12/pi)))
  }
}


adjust.phases<-function(realphi, infphi, h24=FALSE){
  realphi=realphi*12/pi
  infphi=infphi*12/pi
  for(i in 1:length(realphi)){
    if(realphi[i]-infphi[i]>12){infphi[i]=24-infphi[i]}
    if(realphi[i]-infphi[i]<(-12)){infphi[i]=infphi[i]-24}
  }
  if(h24==FALSE){infphi=infphi*pi/12}
  return(infphi)
}