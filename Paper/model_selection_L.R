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
  
  save(OUT.age, file="./Paper/data/OUT/OUT_AGE.RData")
}

if(as.paper){
  OUT.MF=get(load("./Paper/paper_data/OUT_paper/OUT_MF.RData"))
  OUT.age=get(load("./Paper/paper_data/OUT_paper/OUT_AGE.RData"))
  OUT.all=get(load("./Paper/paper_data/OUT_paper/OUT_ALL.RData"))
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
    if(dv=="MF") conds=c(rep("MALE",NS),rep("FEMALE",NS)) else cond=c(rep("YOUNG",NS),rep("OLD",NS))
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
    SS[[tx]]=ss
  }
  save(SS, file=paste("./Paper/data/OUT/SS_", dv, ".RData", sep=""))
}

