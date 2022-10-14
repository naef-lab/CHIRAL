#rm(list=ls())
#gc()
#### Libraries ####
list.of.packages <- c("vroom","ggplot2", "tibble", "gridExtra", "ggrepel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages, force=TRUE)
library(vroom)
library(ggplot2)
library(tibble)
library(gridExtra)
library(ggrepel)
########

#### Functions ####
#Plot the polar density of genes/donors for the out file given 
Plot_density<-function(OUT, phi, R.plot=FALSE, R.df=FALSE,cut=0.1, varz="pval", comp="small", title_param="", compet="", sz=20, th=1){
  phi.df=tibble(phi=rep((1:1000-1)/500*pi,1))
  phi.df$hour=phi.df$phi/pi*12
  phi.df$count=phi.df$phi
  kapp=20
  for(s in 1:length(phi.df$phi)){
    phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(phi))))
  }
  phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2] #normalize in the linear sense
  phi.df$dens=phi.df$count/sqrt(sum(phi.df$count^2)*phi.df$phi[2]/2) #to normalize the radar plot integral to one
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
    else { cat("comp variable can be set to either big or small")
      stop()
    }
    pho=c(pho,phi)
  }
  for(s in 1:length(phi.df$phi)){
    phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(pho))))
  }
  phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2] #normalize in the linear sense
  phi.df$dens=phi.df$count/sqrt(sum(phi.df$count^2)*phi.df$phi[2]/2)#to normalize the radar plot integral to one
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

#Plot the cSVD representation of the given OUT file for the selected genes up to specified component
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


#Plot the profiles of the selected genes in selected tissues and organizes them in a grid
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
    GE=complex(real = GE[,"a"], imaginary =  GE[,"b"])
    names(GE)=gsub("\\|.*$","",gsub("^.*_", "",GL))
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
  if(length(pltlist)==1){f.plot=pltlist[[1]]}
  else f.plot=do.call(grid.arrange,c(grobs=pltlist, as.table=FALSE))
  print(f.plot)
  if(R.plot) return(f.plot)
}

#Calculate the svd for an OUT file
svd.from.out<-function(OUT, inter.genes, ENSG=F){
  gene.c=matrix(0, nrow = length(inter.genes),ncol=length(names(OUT)))
  dimnames(gene.c)=list(inter.genes,names(OUT))
  gtot=NULL
  for(i in names(OUT)){ 
    out=OUT[[i]]
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
    
    gene.c[,i]=A[inter.genes]
  }
  gene.c[is.na(gene.c)]=0
  if(ncol(gene.c)*nrow(gene.c)==0){return(NULL)}
  SVD=svd(gene.c)
  
  for(k in 1:length(SVD$d)){
    mn=sum(SVD$v[,k])
    rot=Conj(mn)/Mod(mn)
    SVD$u[,k]=SVD$u[,k]*rot*max(Mod(SVD$v[,k]))*SVD$d[k]
    SVD$v[,k]=Conj(SVD$v[,k]*rot/max(Mod(SVD$v[,k])))
  }
  
  rownames(SVD$u)=rownames(gene.c)
  if(ENSG){
    gtot=unique(gtot)
    pox=match(rownames(SVD$u), gsub("\\..*$","",gsub("_.*$", "",gtot)))
    rownames(SVD$u)=gsub("\\|.*$","",gsub("^.*_", "",gtot))[pox]
  }
  rownames(SVD$v)=names(OUT)
  return(SVD)
}

########


colroma=vroom("./paper_data/roma.txt",  col_names = FALSE)
colroma$hex=paste("#", as.hexmode(round(colroma$X1*255,0)),as.hexmode(round(colroma$X2*255,0)),as.hexmode(round(colroma$X3*255,0)), sep="")
full_col=vroom("./paper_data/GO_full-colorandum.csv", show_col_types = FALSE)

if(!exists("N.cores")) N.cores = 18 
if(!exists("as.paper")) as.paper=FALSE
if(!exists("use.paper.DIP")) use.paper.DIP=FALSE

if(as.paper) load("./paper_data/OUT_paper/OUT_ALL.RData") else load("./data/OUT/OUT_ALL.RData")
if(use.paper.DIP) load("./paper_data/DIPs.RData") else load("./data/DIPs.RData")
  
dir.create("./plot", showWarnings = FALSE)
dir.create("./plot/Figure1", showWarnings = FALSE)

pdf("./plot/Figure1/Fig1_F.pdf")
Plot_density(OUT, phi,cut=0.2, varz="qval")
dev.off()

full_col=full_col[full_col$Class!="Cells",]
dec_names=full_col$`Short name`
names(dec_names)=full_col$`Full name`
colorandum=full_col$`# color`
names(colorandum)=full_col$`Full name`
full_col$Class[full_col$Class=="Respiratory"]="Metabolic"
nmz=unique(full_col$Class)
sz=20
th=1
pcut=2
qcut=0.2
Rcut=0.5


pdf("./plot/Figure1/Fig1_G.pdf")
  phi.full=NULL
  pho=NULL
  #One density plot per tissue
  for(name in names(OUT)){
    out=OUT[[name]]
    full=out$data.fit
    full=subset(full, pval<pcut)
    full=subset(full, qval<qcut)
    full$R=2*sqrt(full$a^2+full$b^2)
    full=subset(full, R>Rcut)
    phi=full$phase*pi/12
    pho=c(pho,phi)
    phi.df=tibble(phi=rep((1:1000-1)/500*pi,1))
    phi.df$hour=phi.df$phi/pi*12
    phi.df$count=phi.df$phi
    kapp=20
    for(s in 1:length(phi.df$phi)){
      phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(phi))))
    }
    phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2]
    ppt=ggplot(phi.df)+geom_line(aes(x=hour,y=dens), colour=colorandum[name], size=th)+coord_polar()+
      scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ylim(0,NA)+theme_minimal()+
      labs(x="time of day", y="density", title=name)+theme(text = element_text(size=sz))
    print(ppt)
    phi.df$tissue=name
    phi.full=rbind(phi.full, phi.df)
  }
  #One density plot per group of tissues
  for(i in nmz){ 
    ct=which(full_col$Class==i)
    tbt=NULL
    curr_nms=full_col$`Full name`[ct]
    curr_tb=subset(phi.full, tissue %in% curr_nms)
    ppt=ggplot(curr_tb)+geom_line(aes(x=hour,y=dens, colour=tissue), size=th)+coord_polar()+scale_color_manual(values = colorandum)+
      scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ylim(0,NA)+theme_minimal()+
      labs(x="time of day", y="density", title=i)+theme(legend.position = "none", text = element_text(size=sz))
    print(ppt)
  }
  #All tissues together
  ppt=ggplot(phi.full)+geom_line(aes(x=hour,y=dens,colour=tissue), size=th)+coord_polar()+scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ylim(0,NA)+
    theme_minimal()+labs(x="time of day", y="density", title="All Tissues")+scale_color_manual(values=colorandum[names(OUT)])+
    theme(legend.position = "none", text = element_text(size=sz))
  print(ppt)
  dev.off()


val="R"
all=NULL
tbt=NULL


pdf("./plot/Figure1/Fig1_H.pdf")
  for(j in names(OUT)){
    out=OUT[[j]]
    fit=out$data.fit
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
      breaks = seq(0, max(max(full[,val]),1), by=0.01) 
      if(val!="R") breaks=breaks^11
      amp.cut = cut(full[,val], breaks, right=FALSE) 
      freq.cut = table(amp.cut) 
      if(val=="R") freq.cut = rev(freq.cut) 
      cumsum.frq=c(cumsum(freq.cut),nrow(full[,val]))+1
      tbz=tibble(R=breaks[-1], n.genes=cumsum.frq, kind=j)
      if(val=="R") tbz=tibble(R=rev(breaks[-1]), n.genes=cumsum.frq, kind=j)
      tbt=rbind(tbt,tbz)
    }
  }
  #Cumulative plots for all tissues
  if(!is.null(tbt)){
    if (val=="R") print(ggplot(tbt)+geom_line(aes(x=R,y=n.genes, color=kind), size=th)+scale_y_log10(limits=c(1,2500))+
                          labs(title="All tissues",y="# of genes", x="log2(peak to trough)", color="Category")+theme_minimal()+
                          scale_color_manual(values=colorandum[unique(tbt$kind)])+theme(legend.position = "none",text = element_text(size=sz)))
    else print(ggplot(tbt)+geom_line(aes(x=R,y=n.genes, color=kind), size=th)+scale_y_log10()+scale_x_continuous(trans=reverselog_trans(10))+
                 labs(title="All tissues",y="# of genes", x=val)+theme_minimal()+scale_color_manual(values=colorandum[unique(tbt$kind)])+
                 theme(legend.position = "none",text = element_text(size=sz)))
  }
  #Cumulative plots for groups of tissues
  for(i in nmz){ 
    ct=which(full_col$Class==i)
    tbt=NULL
    for(j in full_col$`Full name`[ct]){
      out=OUT[[j]]
      fit=out$data.fit
      inf.phi=out$phi
      exprx=out$E
      gene.list=gsub("\\|.*$","",gsub("^.*_", "",fit$genes))
      clock.coord=sapply(out$geni,function(x){match(x, gene.list)})
      clock.coord=clock.coord[!is.na(clock.coord)]
      full=fit[,c("amp","a","b","R2", "genes","qval", "pval")]
      if (val=="R") {full=subset(full, pval<pcut)}
      full=subset(full, pval<pcut)
      full=subset(full, qval<qcut)
      full$R=2*sqrt(full$a^2+full$b^2)
      full=subset(full, R>Rcut)
      if(nrow(full)>1){
        full$kind=names(OUT)[j]
        all=rbind(all,full)
        breaks = seq(0, max(max(full[,val]),1), by=0.01) 
        if(val!="R") breaks=breaks^11 
        amp.cut = cut(full[,val], breaks, right=FALSE) 
        freq.cut = table(amp.cut) 
        if(val=="R") freq.cut = rev(freq.cut) 
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
                           theme(text = element_text(size=sz), legend.position = "none"))
      else print(ggplot(tbt)+geom_line(aes(x=R,y=n.genes, color=kind, linetype=kind), size=th)+
                   scale_y_log10()+scale_x_continuous(trans=reverselog_trans(10))+labs(title=i,y="# of genes", x=val)+
                   theme_minimal()+scale_color_manual(values=colorandum[unique(tbt$kind)])+theme(text = element_text(size=sz), legend.position = "top"))
    }
  }
  dev.off()


#Plot the SVD of interesting gene groups
lb=6
pt=3
sz=18
gene_inf=get(load("./paper_data/CGRs.RData"))
HSF1_genes=get(load("./paper_data/HSF1_g.RData"))
{
  Plot_cSVD(OUT, gene_inf, full_col,loc ="./plot/Figure1/Fig1_B.pdf" , CT=15, dot_size =pt, label_size = lb, text_size = sz)
  Plot_cSVD(OUT, HSF1_genes, full_col,loc ="./plot/Figure1/Fig1_D.pdf" , CT=5, dot_size =pt, label_size = lb, text_size = sz, ymax=0.6)
}
#Plot the various profiles
{
  pdf("./plot/Figure1/Fig1_C.pdf")
  tissuex=c("Brain - Cortex", "Adipose - Visceral (Omentum)")
  geni=gene_inf[c(2,8)]
  Plot_profiles(OUT, tissuex, geni, val="qval")
  dev.off()
}

{
  pdf("./plot/Figure1/Fig1_E.pdf")
  tissuex=c("Brain - Amygdala", "Spleen")
  geni=c("HSPH1", "HSPA1B")
  Plot_profiles(OUT, tissuex, geni, val="qval")
  dev.off()
}
