source("/scratch/For_cedric_with_pizza/model_selection/nconds_functions.R")
source("/scratch/For_cedric_with_pizza/model_selection/nconds.R")

source("/data/shared/FelixProjects/XYmodels/em-useful.R")
source("/scratch/For_cedric_with_pizza/clean_script/f24_R2_cycling.R")

dir.create(file.path("./Paper/out_MS"), showWarnings = FALSE)

as.paper=TRUE

if(!as.paper){
  
  phio=phi=get(load("./Paper/data/DIPs.RData"))
  
  dat.raw=get(load("./Paper/data/CPM_full.RData"))
  
  E=CPM_to_E(CPM.all.norm.large)
  
  E=split_E_sex(E, samp)
  
  OUT.MF=Make_big_OUT(E, phi)
  
  OUT.MF=Fit_OUT(OUT.MF)
  
  save(OUT.MF, file="./Paper/data/OUT/OUT_MF.RData")
  
  E=split_E_age(E, samp)
  
  OUT.age=Make_big_OUT(E, phi)
  
  OUT.age=Fit_OUT(OUT.age)
  
  save(OUT.age, file="./Paper/data/OUT/OUT_age_n5.RData")
}

if(as.paper){
  OUT.MF=get(load("./Paper/paper_data/OUT_paper/OUT_MF.RData"))
  OUT.age=get(load("./Paper/paper_data/OUT_paper/OUT_age_n5.RData"))
  OUT.all=get(load("./Paper/paper_data/OUT_paper/OUT_all.RData"))
  dat.raw=get(load("./Paper/paper_data/CPM_full.RData"))
  phi=get(load("./Paper/paper_data/DIPs.RData"))
}

Meta=read.table("./Paper/paper_data/sample_metadata.txt", header=TRUE)


qcut=0.2
Rcut=0.5


#####

for(dv in c("MF", "age_n5")){
  if(dv=="MF"){
    OUT= OUT.MF 
    MU=subset(Meta, mfSS==1)
  }
  else{
    OUT= OUT.age
    MU=subset(Meta, ageSS==1)
  }
  tix.c=gsub("-old", "",gsub("-young", "",gsub("-Female", "",gsub("-Male", "", names(OUT)))))
  
  
  SS=list()
  for (tx in tix.c[duplicated(tix.c)]){
    MT=subset(MU, tissue==tx)
    idx=which(tix.c==tx)
    nms=gsub("^.*-", "", names(OUT)[idx])
    nm1=ifelse(dv=="MF", paste(tx,"-Male", sep=""),  paste(tx,"-young", sep=""))
    nm2=ifelse(dv=="MF", paste(tx,"-Female", sep=""),  paste(tx,"-old", sep=""))
    T1=OUT[[nm1]]$E
    T2=OUT[[nm2]]$E
    P1=OUT[[nm1]]$phi
    P2=OUT[[nm2]]$phi
    S1=match(intersect(MT$fullID, colnames(T1)), colnames(T1))
    S2=match(intersect(MT$fullID, colnames(T2)), colnames(T2))
    #CG=intersect(rownames(T1), rownames(T2))
    NS=length(S1)
    IRN=intersect(rownames(T1),rownames(T2))
    FM=cbind(T1[IRN,S1], T2[IRN,S2])
    FP=c(P1[S1],P2[S2])
    raw.E=dat.raw[[tx]]
    raw.E=raw.E[IRN,colnames(FM)]
    ii=which(rowMeans(raw.E) > 2 & apply(raw.E,1,function(x) length(x[x<0])) <20 )
    FM=FM[ii,]
    if(dv=="MF") conds=c(rep("MA",NS),rep("FE",NS)) else cond=c(rep("YOUNG",NS),rep("OLD",NS))
    dat=nconds(FM,conds=conds,t=FP*12/pi,out.prefix = "./Paper/out_MS/test_MF.pdf")
    colnames(dat)[colnames(dat)=="BICW"]="AICW"
    ss=dat
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
    SS[[i]]=ss
  }
  save(SS, file=paste("./Paper/data/OUT/SS_", dv, ".RData", sep=""))
}

sso=SS.MF[[tx]]
for (nnn in names(SS)) {
  ss=SS[[nnn]]
  sso=SS.MF[[nnn]]
  for (i in ncol(ss)) {
    if(max(abs(ss[i]-sso[i]))>0.0001){
      cat(i)
      stop()
    }
    
  }
  
}




dplx=names(OUT)[duplicated(gsub())]

set.seed(42)
phi=get(load("./scratch/For_cedric_with_pizza/phi_all_musc.RData"))
phi=phi%%(2*pi)

dat.2=get(load("/scratch/For_cedric_with_pizza/OUT_files/OUT_big_mean_final_MF_musc.RData"))
dat.raw=get(load("/scratch/For_cedric_with_pizza/data/CPM_full.RData"))

all.tiss=lapply(dat.2,function(x) x$tissue)

tiss.2= unique(names(which(table(unlist(all.tiss))==2)))
tiss.2=tiss.2[-grep('Cells',(tiss.2))]

for(k in tiss.2){
  dat.raw.tiss=dat.raw[[k]]
  posi=which(is.na(match(all.tiss,k))==F)
  
  MA=dat.2[[posi[1]]]
  FE=dat.2[[posi[2]]]
  
  A=MA$data.fit[,grep('GTEX',colnames(MA$data.fit))]
  B=FE$data.fit[,grep('GTEX',colnames(FE$data.fit))]
  
  
  nami.MA=sapply(strsplit(colnames(A),split=".",fixed=T),"[[",2)
  nami.FE=sapply(strsplit(colnames(B),split=".",fixed=T),"[[",2)
  colnames(A)=nami.MA
  colnames(B)=nami.FE
  
  ii.MA=intersect(nami.MA,names(phi))
  ii.FE=intersect(nami.FE,names(phi))
  
  
  A=A[,ii.MA]
  B=B[,ii.FE]
  t.m=24*phi[ii.MA]/(2*pi)
  t.f=24*phi[ii.FE]/(2*pi)
  
  mm=which.min(c(ncol(A),ncol(B)))
  
  if(mm==2){
    ii=sort(sample(ncol(A),ncol(B)))
    A=A[,ii]
    t.m=t.m[ii]
  }
  if(mm==1){
    ii=sort(sample(ncol(B),ncol(A)))
    B=B[,ii]
    t.f=t.f[ii]
  }
  AB=cbind(A,B)
  t.mf=c(t.m,t.f)
  
  dat.raw.tiss=dat.raw.tiss[,grep('GTEX',colnames(dat.raw.tiss))]
  colnames(dat.raw.tiss)=sapply(strsplit(colnames(dat.raw.tiss),split=".",fixed=T),"[[",2)
  raw.E=dat.raw.tiss[rownames(AB),colnames(AB)]
  ii=which(rowMeans(raw.E) > 2 & apply(raw.E,1,function(x) length(x[x<0])) <20 )
  AB=AB[ii,]
  dat=nconds(AB,conds=c(rep("MA",ncol(A)),rep("FE",ncol(B))),t=t.mf,out.prefix = "/home/cgobet/lorenzo/plot/test_MF.pdf")
  save(dat,file=paste0("/home/cgobet/lorenzo/data/final_27082020_MF/",k,"_common_phase_corr_fit_MF__AIC_subsample_NA5.RData"))
  
} 


















SS=SS.MF
OUT=OUT.MF
SS=SS.age
OUT=OUT.age
qcut=0.2
Rcut=0.5
Gtix=list()
for (i in names(SS)) {
  gn=NULL
  ss=SS[[i]]
  nms=colnames(ss)
  idxs=which(grepl(i, names(OUT)))
  for (j in idxs) {
    nm=names(OUT)[j]
    div=gsub("^.*-", "", nm)
    out=OUT[[nm]]
    pvals=out$data.fit[,c("qval","amp")]
    rownames(pvals)=out$data.fit[,"genes"]
    pvalus=subset(pvals, qval<qcut & amp>Rcut)
    gn=c(gn, rownames(pvalus))
    
  }
  out=OUT.all[[i]]
  pvals=out$data.fit[,c("qval","amp")]
  rownames(pvals)=out$data.fit[,"genes"]
  pvalus=subset(pvals, qval<qcut & amp>Rcut)
  gn=c(gn, rownames(pvalus))
  Gtix[[i]]=unique(gn)
  ss=ss[intersect(Gtix[[i]], rownames(ss)),]
  SS[[i]]=ss
}

ss=SS[["Liver"]]
table(ss$model)

SS=SS.MF
OUT=OUT.MF
SS=SS.age
OUT=OUT.age
Gtix=list()
for (i in names(SS)){
  gn=NULL
  ss=SS[[i]]
  nms=colnames(ss)
  it=gsub("\\(.*", "", i)
  idxs=which(grepl(it, names(OUT)))
  cat(length(idxs))
  for (j in idxs){
    nm=names(OUT)[j]
    div=gsub("^.*-", "", nm)
    out=OUT[[nm]]
    pvals=out$data.fit[,c("qval","amp")]
    rownames(pvals)=out$data.fit[,"genes"]
    pvalus=subset(pvals, qval<qcut & amp>Rcut)
    gn=c(gn, rownames(pvalus))
    coms=intersect(rownames(pvals), rownames(ss))
    ss$qval=1
    ss$R=0
    ss[coms, c("qval", "R")]=pvals[coms,]
    colnames(ss)[(length(colnames(ss))-1):length(colnames(ss))]=paste(c("qvals", "R"),div, sep="_")
  }
  div="all"
  out=OUT.all[[i]]
  pvals=out$data.fit[,c("qval","amp")]
  rownames(pvals)=out$data.fit[,"genes"]
  pvalus=subset(pvals, qval<qcut & amp>Rcut)
  gn=c(gn, rownames(pvalus))
  coms=intersect(rownames(pvals), rownames(ss))
  ss$qval=1
  ss$R=0
  ss[coms, c("qval", "R")]=pvals[coms,]
  colnames(ss)[(length(colnames(ss))-1):length(colnames(ss))]=paste(c("qvals", "R"),div, sep="_")
  Gtix[[i]]=unique(gn)
  ss$accepted=0
  ss$accepted[match(Gtix[[i]],rownames(ss))]=1
  ss$model.c=ss$model*ss$accepted
  ss=ss[,-(ncol(ss)-1)]
  SS[[i]]=ss
}


SS.MF=SS
#save(SS.MF, file="/scratch/For_cedric_with_pizza/28042022/SS_MF.RData")
SS.age=SS
#save(SS.age, file="/scratch/For_cedric_with_pizza/28042022/SS_age_n5.RData")
