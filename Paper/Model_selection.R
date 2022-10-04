source("./nconds_functions.R")
source("./nconds.R")

#### Functions ####

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


########

#### Main ####
N.cores = 18

dir.create(file.path("./Paper/out_MS"), showWarnings = FALSE)

as.paper=FALSE

if(!as.paper){
  
  phio=phi=get(load("./data/DIPs.RData"))
  
  dat.raw=get(load("./data/CPM/CPM_full.RData"))
  
  E=CPM_to_E(CPM.all.norm.large)
  
  E=split_E_sex(E, samp)
  
  OUT.MF=Make_big_OUT(E, phi)
  
  OUT.MF=Fit_OUT(OUT.MF, N.cores = N.cores)
  
  save(OUT.MF, file="./data/OUT/OUT_MF.RData")
  
  E=split_E_age(E, samp)
  
  OUT.age=Make_big_OUT(E, phi)
  
  OUT.age=Fit_OUT(OUT.age, N.cores = N.cores)
  
  save(OUT.age, file="./data/OUT/OUT_AGE.RData")
}

if(as.paper){
  OUT.MF=get(load("./paper_data/OUT_paper/OUT_MF.RData"))
  OUT.age=get(load("./paper_data/OUT_paper/OUT_AGE.RData"))
  OUT.all=get(load("./paper_data/OUT_paper/OUT_ALL.RData"))
  dat.raw=get(load("./paper_data/CPM_full.RData"))
  phi=get(load("./paper_data/DIPs.RData"))
}

Meta=read.table("./paper_data/sample_metadata.txt", header=TRUE)


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
    NS=length(S1)
    IRN=intersect(rownames(T1),rownames(T2))
    FM=cbind(T1[IRN,S1], T2[IRN,S2])
    FP=c(P1[S1],P2[S2])
    raw.E=dat.raw[[tx]]
    raw.E=raw.E[IRN,colnames(FM)]
    ii=which(rowMeans(raw.E) > 2 & apply(raw.E,1,function(x) length(x[x<0])) <20 )
    FM=FM[ii,]
    if(dv=="MF") conds=c(rep("MALE",NS),rep("FEMALE",NS)) else cond=c(rep("YOUNG",NS),rep("OLD",NS))
    dat=nconds(FM,conds=conds,t=FP*12/pi, out.prefix = NULL)
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
    SS[[tx]]=ss
  }
  save(SS, file=paste("./data/OUT/SS_", dv, ".RData", sep=""))
}

