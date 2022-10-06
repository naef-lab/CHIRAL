rm(list=ls())
gc()
### Libraries ###
source("./nconds_functions.R")
source("./nconds.R")
list.of.packages <- c("lmtest", "data.table", "parallel", "combinat")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages, force=TRUE)
library(lmtest)
library(data.table)
library(parallel)
library(combinat)
#### Functions ####
#Divide E object by sex of donors
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
      A[[paste(i, "MALE", sep="-")]]=e
    }
    fidx=match(females, cn)
    fidx=fidx[!is.na(fidx)]
    if (length(fidx)>24){
      e$E=ee[,fidx]
      A[[paste(i, "FEMALE", sep="-")]]=e
    }
  }
  return(A)
}
#Divide E object by age of donors
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
      A[[paste(i, "YOUNG", sep="-")]]=e
    }
    fidx=match(females, cn)
    fidx=fidx[!is.na(fidx)]
    if (length(fidx)>24){
      e$E=ee[,fidx]
      A[[paste(i, "OLD", sep="-")]]=e
    }  }
  return(A)
}

spliti= function(x, sp, nb){
  v=sapply(strsplit(x,sp),"[[",nb)
  return(v)
}

# Reorganizes the CPM matrix into a list (E), removing unwanted tissues
CPM_to_E<-function(CPM.all.norm, min.samp=24, sep=NULL, samp=NULL){
  E.matrix=list()
  if(length(which(startsWith(names(CPM.all.norm), "Cells")))>0) CPM.all.norm=CPM.all.norm[-which(startsWith(names(CPM.all.norm), "Cells"))]
  if(length(which(startsWith(names(CPM.all.norm), "Whole")))>0) CPM.all.norm=CPM.all.norm[-which(startsWith(names(CPM.all.norm), "Whole"))]
  if(!is.null(samp)){
    if(tolower(sep)=="sex"){
      male=samp$sub.id[samp$SEX==1]
      female=samp$sub.id[samp$SEX==2]
      for(name in names(CPM.all.norm)){
        E=list()
        E$E=CPM.all.norm[[name]]
        E$E=E$E[,!is.na(E$E[1,])]
        if(!is.null(E$E)){
          E$tissue=name
          ml=fm=list()
          colnames(E$E)=spliti(colnames(E$E), "\\.", 2)
          ml$E=E$E[,colnames(E$E) %in% male]
          ml$type="GTEX-male"
          ml$tissue=name
          if(ncol(ml$E)>=min.samp){E.matrix[[paste(name, "MALE", sep="-")]]=ml}
          fm$E=E$E[,colnames(E$E) %in% female]
          fm$type="GTEX-female"
          fm$tissue=name
          if(ncol(fm$E)>=min.samp){E.matrix[[paste(name, "FEMALE", sep="-")]]=fm}
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
          ml$E=E$E[,colnames(E$E) %in% young]
          ml$type="GTEX-old"
          ml$tissue=name
          if(ncol(ml$E)>=min.samp){E.matrix[[paste(name, "YOUNG", sep="-")]]=ml}
          fm$E=E$E[,colnames(E$E) %in% old]
          fm$type="GTEX-young"
          fm$tissue=name
          if(ncol(fm$E)>=min.samp){E.matrix[[paste(name, "OLD", sep="-")]]=fm}
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
          E.matrix[[name]]=E}
      }
    }
  }
  
  return(E.matrix)
}
#Create the OUT structured file from any CPM and set of phases
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

#Fit each gene with Harmonic regression in all tissues
Fit_OUT<-function(OUT,period=24, NA5=T, N.cores){
  for(i in names(OUT)){ 
    E=OUT[[i]]$E
    de=sweep(E,1,rowMeans(E),FUN="-")
    if(NA5) E[abs(de)>5]=NA
    phase=OUT[[i]]$phi
    if(length(phase)<24){
      OUT[[i]]=NULL
      next
    }
    
    genes=as.list(rownames(E))
    names(genes)=genes
    dat.fit=mclapply(genes,function(x) harm_reg(as.numeric(E[x,]), 12*as.numeric(phase)/pi, period=period),mc.cores = N.cores)
    dat.fit=do.call(rbind,dat.fit)
    dat.fit=as.data.frame(dat.fit)
    
    dat.fit$qval=p.adjust(dat.fit$pval, "BH")
    OUT[[i]]$data.fit=data.frame(dat.fit,E,genes=rownames(E))
    colnames(OUT[[i]]$data.fit)=c("pval","phase", "amp","mu", "a","b", "period", "R2","qval", colnames(E), "genes")
    OUT[[i]]$E=OUT[[i]]$data.fit[,colnames(E)]
  }
  return(OUT)
}
#Harmonic regression funcion
harm_reg<-function(x, t, period){
  n=length(x)
  fit0=lm(x~1)
  c=cos(2*pi*t/period)
  s=sin(2*pi*t/period)
  fit1=lm(x~c+s)
  mu=coef(fit1)[1]
  a=coef(fit1)[2]
  b=coef(fit1)[3]
  p.val=lrtest(fit1, fit0)$Pr[2]
  amp=2*sqrt(a^2+b^2)
  phase=atan2(b,a)%%(2*pi)
  phase=period*phase/(2*pi)
  return(c(pval=p.val,phase=phase,amp=amp,mu=mu,a=a,b=b, period=period, R2=summary(fit1)$r.squared))
}

########

#### Main ####
if(!exists("N.cores")) N.cores = 18 
if(!exists("as.paper")) as.paper=FALSE
if(!exists("use.paper.DIP")) use.paper.DIP=FALSE

if(!file.exists(file="./paper_data/raw/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")){ 
  samp = fread('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
  fwrite(samp, file = "./paper_data/raw/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
}else{
  samp = fread("./paper_data/raw/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
}
samp$sub.id=spliti(samp$SUBJID,"-",2)
samp$AGE=as.numeric(spliti(samp$AGE,"-",1))+5

#Crate appropriate age and sex dependent files for subsequent analysis
#It's a bit slow but its normal

if(use.paper.DIP)   phi=get(load("./paper_data/DIPs.RData")) else phi=get(load("./data/DIPs.RData"))

if(!as.paper){
  
  dat.raw=get(load("./data/CPM/CPM_full.RData"))
  CPM.all.norm.large = get(load("./data/CPM/CPM.all.norm_large.RData"))
  E=CPM_to_E(CPM.all.norm.large)
  
  E.sex=split_E_sex(E, samp)
  
  OUT.MF=Make_big_OUT(E.sex, phi)
  
  OUT.MF=Fit_OUT(OUT.MF, N.cores = N.cores)
  
  save(OUT.MF, file="./data/OUT/OUT_MF.RData")
  
  E.age=split_E_age(E, samp)
  
  OUT.age=Make_big_OUT(E.age, phi)
  
  OUT.age=Fit_OUT(OUT.age, N.cores = N.cores)
  
  save(OUT.age, file="./data/OUT/OUT_AGE.RData")
  OUT.all=get(load("./data/OUT/OUT_ALL.RData"))
  
}

if(as.paper){
  OUT.MF=get(load("./paper_data/OUT_paper/OUT_MF.RData"))
  OUT.age=get(load("./paper_data/OUT_paper/OUT_AGE.RData"))
  OUT.all=get(load("./paper_data/OUT_paper/OUT_ALL.RData"))
  dat.raw=get(load("./paper_data/CPM_full.RData"))
}

Meta=read.table("./paper_data/sample_metadata.txt", header=TRUE)

qcut=0.2 # BH corrected p-value > 0.2
Rcut=0.5 # log2 peak trough > 0.5


#####
#Create SS (SubSample) files and perform the model selection on all genes; then verify if the pass the qvalue and R treshold
for(dv in c("MF", "AGE")){
  if(dv=="MF"){
    OUT= OUT.MF 
    MU=subset(Meta, mfSS==1)
  }else{
    OUT= OUT.age
    MU=subset(Meta, ageSS==1)
  }
  tix.c=gsub("-OLD", "",gsub("-YOUNG", "",gsub("-FEMALE", "",gsub("-MALE", "", names(OUT)))))
  
  
  SS=list()
  for (tx in tix.c[duplicated(tix.c)]){
    MT=subset(MU, tissue==tx)
    MT$fullID=sapply(strsplit(as.character(unlist(MT$fullID)),".",fixed=TRUE),"[[",2)
    idx=which(tix.c==tx)
    nms=gsub("^.*-", "", names(OUT)[idx])
    nm1=ifelse(dv=="MF", paste(tx,"-MALE", sep=""),  paste(tx,"-YOUNG", sep=""))
    nm2=ifelse(dv=="MF", paste(tx,"-FEMALE", sep=""),  paste(tx,"-OLD", sep=""))
    
    T1=OUT[[nm1]]$E
    T2=OUT[[nm2]]$E
    P1=OUT[[nm1]]$phi
    P2=OUT[[nm2]]$phi
    S1=match(intersect(MT$fullID, colnames(T1)), colnames(T1))
    S2=match(intersect(MT$fullID, colnames(T2)), colnames(T2))
    NS=length(S1)
    IRN=intersect(rownames(T1),rownames(T2))
    FM=cbind(T1[IRN,S1], T2[IRN,S2])
    FP=c(P1[S1],P2[S2])
    raw.E=dat.raw[[tx]]
    colnames(raw.E)=sapply(strsplit(as.character(unlist(colnames(raw.E))),".",fixed=TRUE),"[[",2)
    raw.E=raw.E[IRN,colnames(FM)]
    ii=which(rowMeans(raw.E) > 2 & apply(raw.E,1,function(x) length(x[x<0])) <20 )
    FM=FM[ii,]
    if(dv=="MF") conds=c(rep("MALE",NS),rep("FEMALE",NS)) else conds=c(rep("YOUNG",NS),rep("OLD",NS))
    ss=nconds(FM,conds=conds,t=FP*12/pi, out.prefix = NULL, N.cores = N.cores)

    gn=NULL
    for(nm in c(nm1, nm2)){
      dvt=gsub("^.*-", "", nm)
      out=OUT[[nm]]
      pvals=out$data.fit[,c("qval","amp")]
      rownames(pvals)=out$data.fit[,"genes"]
      pvalus=subset(pvals, qval<qcut & amp>Rcut)
      gn=c(gn, rownames(pvalus))
      coms=intersect(rownames(pvalus), rownames(ss))
      ss$qval=1
      ss$R=0
      ss[coms, c("qval", "R")]=pvals[coms,]
      colnames(ss)[(length(colnames(ss))-1):length(colnames(ss))]=paste(c("qvals", "R"), dvt, sep="_")
    }
    dvt="all"
    out=OUT.all[[tx]]
    pvals=out$data.fit[,c("qval","amp")]
    rownames(pvals)=out$data.fit[,"genes"]
    pvalus=subset(pvals, qval<qcut & amp>Rcut)
    gn=c(gn, rownames(pvalus))
    coms=intersect(rownames(pvals), rownames(ss))
    ss$qval=1
    ss$R=0
    ss[coms, c("qval", "R")]=pvals[coms,]
    colnames(ss)[(length(colnames(ss))-1):length(colnames(ss))]=paste(c("qvals", "R"),dvt, sep="_")
    gnu=unique(gn)
    ss$accepted=0
    ss$accepted[match(gnu,rownames(ss))]=1
    ss$model.c=ss$model*ss$accepted
    ss=ss[,-(ncol(ss)-1)]
    SS[[tx]]=ss
  }
  save(SS, file=paste("./data/OUT/SS_", dv, ".RData", sep=""))
}

