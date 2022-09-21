#file of 01/07/2020
gene_inf<-c("DBP"   ,  "PER3"  ,  "TEF"   ,  "NR1D2" ,  "PER1" ,   "PER2"  ,  "NPAS2" ,  "ARNTL" ,  "NR1D1",  "CRY1"  ,  "CRY2","CIART" )
setwd("/scratch/For_cedric_with_pizza/git_scripts")
load("/scratch/For_cedric_with_pizza/git_scripts/data/CPM.all.norm.RData")
source("/scratch/For_cedric_with_pizza/git_scripts/CHIRAL/functions_1.R")
E.matrix=list()
#phase inferenece tissue by tissue
for(name in names(CPM.all.norm)){
  E=list()
  E$E=CPM.all.norm[[name]]
  E$E=E$E[,!is.na(E$E[1,])]
  if(!is.null(E$E)){
    E$tissue=name
    E$type="GTEx"
    #gene="PER3"
    #E$E=as.matrix(subset(E$E,rowMeans(E$E)>0))
    #cat(name, " ", dim(E$E), gene, "is present", any(gsub("^.*_", "",rownames(E$E))==gene),  "\n")
    E.matrix[[name]]=E}
}

OUT= mclapply(E.matrix,infer_l, gene_inf, mc.cores=16)

#phase averaging -> 1 phase per individual
{
  all_people=NULL
  for(i in names(OUT)){
    out=OUT[[i]]
    out=order.out.setgene(out)
    OUT[[i]]=out
    sampz=gsub("\\..*$","", gsub("GTEX.","", colnames(out$E)))
    all_people=c(all_people, sampz)
    cat(i,"ordered: ", out$has.been.flipped, "\n")
  }
  people=unique(all_people)
  phi_mat=matrix(nrow = length(people), ncol= length(names(OUT)))
  dimnames(phi_mat)=list(people, names(OUT))
  for(i in names(OUT)){
    out=OUT[[i]]
    sampz=gsub("\\..*$","", gsub("GTEX.","", colnames(out$E)))
    phi_mat[sampz,i]=out$phi
  }
  dim(phi_mat)
  tix.study=names(OUT)
  phi.study=phi_mat[,tix.study]
  phi.comp=phi.study
  for(i in 1:ncol(phi.study)){
    phi.comp[,i]=complex(argument=phi.study[,i])
  }
  #hist(Arg(phi.comp[4,]))
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
        pc[k]=sum(Mod(pp+pp[k])>1.95)
      }
      kstar=which.max(pc)
      gp=sum(pp[Mod(pp+pp[kstar])>1.95])
      phis[i]=gp/Mod(gp)}
  }
  phit=phis
}
phi=Arg(phit)%%(2*pi)




