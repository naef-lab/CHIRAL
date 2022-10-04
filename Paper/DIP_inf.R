rm(list=ls())

 
##### Functions ####
spliti= function(x, sp, nb){
  v=sapply(strsplit(x,sp),"[[",nb)
  return(v)
}

# CPM to E
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
          E.matrix[[name]]=E}
      }
    }
  }
  
  return(E.matrix)
}

#infer with CHIRAL 
infer_l=function(k,clockgenes=NULL){
  tix=k$tissue
  tp=k$type
  v=k$E
  out=CHIRAL(v, 500, clockgenes = clockgenes, standardize = TRUE, GTEx_names=TRUE) 
  out$tissue=tix
  out$type=tp
  return(out)
}

#Harmonic regression
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

#Fit each gene with Harmonic regression
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

Unique_donors<-function(OUT){
  all_people=NULL
  for(i in names(OUT)){
    sampz=colnames(OUT[[i]]$E)
    all_people=c(all_people, sampz)
  }
  people=unique(all_people)
  return(people)
}
order.from.hc.ref<-function(x){
  sx=x
  #hard coded reference
  full.ref=c(-0.43043911+0.02646768i, -0.21520829-0.08481562i, -0.18976262-0.06323798i, -0.26402067+0.03233651i, 
             -0.15460540-0.17813898i,  0.34503079+0.12988785i,  0.39656252+0.03842312i, -0.30343360+0.29074942i,  0.06226694-0.12542420i, -0.06159846-0.04409744i, -0.18018038-0.01939843i)
  names(full.ref)=c("DBP"   ,  "PER3"  ,  "TEF"   ,  "NR1D2"  ,   "PER2"  ,  "NPAS2" ,  "ARNTL" ,  "NR1D1"  , "CRY1"  ,  "CRY2",  "PER1")
  
  shifts=c(1:1000)*pi/500
  ts=complex(real=cos(shifts), imaginary=sin(shifts))
  common=intersect(names(full.ref), names(x))
  ref.mat=full.ref[common]%o%ts
  x=x[common]
  x=x/(sqrt(sum(Re(x*Conj(x)))))
  gen.p.scal=max(Re(t(ref.mat)%*%Conj(x))) #remember the inversion 
  inv.p.scal=max(Re(t(ref.mat)%*%x))
  
  if(gen.p.scal>inv.p.scal){return(sx)}
  return(Conj(sx))
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
  good.s=names(phit[score>0.25])
  phi.fin=phi_mat[good.s,]
  samp.tix=nrow(phi.fin)-colSums(is.na(phi.fin))
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
########

dir.create(file.path("./data/OUT"), showWarnings = FALSE)
as.paper=FALSE

samp <- fread('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
samp$sub.id=spliti(samp$SUBJID,"-",2)
samp$AGE=as.numeric(spliti(samp$AGE,"-",1))+5

if(as.paper){
  CPM.all.norm.large=get(load("./paper_data/CPM/CPM.all.norm_large.RData"))
  CPM.all.norm=get(load("./paper_data/CPM/CPM.all.norm.RData"))
}

if(!as.paper){
  CPM.all.norm.large=get(load("./data/CPM/CPM.all.norm_large.RData"))
  CPM.all.norm=get(load("./data/CPM/CPM.all.norm.RData"))
}

E=CPM_to_E(CPM.all.norm)

gene_inf=get(load("./paper_data/CGRs.RData"))

OUT=mclapply(E, infer_l, gene_inf, mc.cores=16)
OUT=Fit_OUT(OUT, N.cores=16)
OUT=Set_OUT(OUT)

donors=Unique_donors(OUT)

phi_matrix=Create_phi_matrix(OUT, donors)

phi_paper=get(load("./paper_data/DIPs.RData"))

phi=Phi.from.phi_mat(phi_matrix, ct=1.95) #these phases will not be identical to the DIP for the various existing stochastic steps

save(phi, file="./data/DIPs.RData")

ii = intersect(names(phi_paper),names(phi))
qplot(phi_paper[ii],phi[ii])

E=CPM_to_E(CPM.all.norm.large)

if(as.paper){phi=phi_paper}

OUT=Make_big_OUT(E, phi)

OUT=Fit_OUT(OUT, N.cores=16)

save(OUT, file="./data/OUT/OUT_ALL.RData")


