#rm(list=ls())
#gc()
### Libraries ###
list.of.packages <- c("vroom", "ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages, force=TRUE)
library(vroom)
library(ggplot2)
### Functions ###
spliti= function(x, sp, nb){
  v=sapply(strsplit(x,sp),"[[",nb)
  return(v)
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
    if(any(endsWith(rownames(SVD$v),"MALE")) || any(endsWith(rownames(SVD$v),"YOUNG"))){
      colz=gsub("-OLD", "",gsub("-YOUNG", "",gsub("-FEMALE", "",gsub("-MALE", "", rownames(SVD$v)))))
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

#Plots cumulative distribution of genes whose rhythms are validated by model selection
Plot_MS_cumulative<-function(OUT, SS, MS=T, div="MF", val="R", pth=NULL, strict=F, qcut=0.2, size=20, th=1){
  pox=unique(gsub("^.*-", "", names(OUT)))
  nn=names(OUT)[duplicated(gsub("-OLD", "",gsub("-YOUNG", "",gsub("-FEMALE", "",gsub("-MALE", "", names(OUT))))))]
  nn=gsub("-OLD", "",gsub("-YOUNG", "",gsub("-FEMALE", "",gsub("-MALE", "", nn))))
  OUT=OUT[c(paste(nn,pox[1], sep="-"),paste(nn,pox[2], sep="-"))]
  if(is.null(pth)){
    if(div=="MF") pdf("./Figure2/Fig2_C-F-I.pdf")
    if(div=="AGE")pdf("./Figure3/Fig3_B-F-I.pdf")
  }else{
    pdf(pth)
  }
  
  all=NULL
  tbt=NULL
  
  id=list(div1=which(endsWith(names(OUT),paste("-", pox[1], sep=""))),div2=which( endsWith(names(OUT),paste("-", pox[2], sep=""))))
  for (ik in names(id) ) {
    fulll=NULL
    for(j in names(OUT)[id[[ik]]]){
      out=OUT[[j]]
      fit=out$data.fit
      jj=gsub("-OLD", "",gsub("-YOUNG", "",gsub("-FEMALE", "",gsub("-MALE", "", j))))
      md=l=gsub("^.*-", "", j)
      if(strict==FALSE){
        if(md %in% c("OLD", "FEMALE")){mds=c(3,4,5)}
        else{mds=c(2,4,5)}}
      else{
        if(md %in% c("OLD", "FEMALE")){mds=c(3)}
        else{mds=c(2)}}
      idt=paste("amp_", md, sep="")
      
      #interesting.genes=inter.genes
      inf.phi=out$phi
      exprx=out$E
      full=fit[,c("amp","a","b","R2", "genes","qval", "pval")]
      ss=SS[[jj]]
      sg=rownames(ss)[which(ss$model.c %in% mds)]
      full=full[sg,]
      full=cbind(full, ss[sg, idt])
      colnames(full)[ncol(full)]="amp_C"
      #if (val=="R") {full=subset(full, qval<qcut)}
      full$R=2*sqrt(full$a^2+full$b^2)
      fulll=rbind(fulll, full)
    }
    full=fulll
    if(nrow(full)>1){
      full$kind=names(OUT)[j]
      all=rbind(all,full)
      breaks = seq(0, max(max(full[,val])+0.1,1), by=0.01)
      if(val %in% c("pval","qval")) breaks=breaks^7
      amp.cut = cut(full[,val], breaks, right=FALSE)
      freq.cut = table(amp.cut)
      if(val %in% c("R", "amp_C")) freq.cut = rev(freq.cut) #also tbz=tibble(R=rev(breaks[-1]), n.genes=cumsum.frq, kind=names(OUT)[j])
      cumsum.frq=c(cumsum(freq.cut),nrow(full[,val]))+1
      tbz=tibble(R=breaks[-1], n.genes=cumsum.frq, kind=l)
      if(val%in% c("R", "amp_C")) tbz=tibble(R=rev(breaks[-1]), n.genes=cumsum.frq, kind=l)
      tbt=rbind(tbt,tbz)
    }
  }
  if(!is.null(tbt)){
    if(val %in% c("R", "amp_C")) print(ggplot(tbt)+geom_line(aes(x=R,y=n.genes, color=kind, linetype=kind), size=th)+
                                         scale_y_log10()+labs(title="Full",y="# of genes", x="log2(peak to trough)", color="", linetype="")+
                                         theme_minimal()+scale_color_manual(values=c("#008ECC","#111E6C","#1034A6","darkred"))+
                                         theme(text = element_text(size=sz)))#+theme(legend.position = "none"))
    else print(ggplot(tbt)+geom_line(aes(x=R,y=n.genes, color=kind, linetype=kind), size=th)+scale_y_log10()+
                 scale_x_continuous(trans=reverselog_trans(10))+labs(title="Full",y="# of genes", x=val)+
                 theme_minimal()+scale_color_manual(values=c("#008ECC","#111E6C","#1034A6","darkred"))+
                 theme(text = element_text(size=sz)))#+theme(legend.position = "none"))
  }
  for(i in nmz){
    ct=which(full_col$Class==i)
    tbt=NULL
    if(length(intersect(full_col$`Full name`[ct],nn))>0){
      for(l in pox){
        fulll=NULL
        for(jdl in intersect(full_col$`Full name`[ct],nn)){
          j=paste(jdl, l, sep="-")
          out=OUT[[j]]
          fit=out$data.fit
          jj=gsub("-OLD", "",gsub("-YOUNG", "",gsub("-FEMALE", "",gsub("-MALE", "", j))))
          md=gsub("^.*-", "", j)
          if(strict==FALSE){
            if(md %in% c("OLD", "FEMALE")){mds=c(3,4,5)}
            else{mds=c(2,4,5)}}
          else{
            if(md %in% c("OLD", "FEMALE")){mds=c(3)}
            else{mds=c(2)}}
          #interesting.genes=inter.genes
          idt=paste("amp_", md, sep="")
          inf.phi=out$phi
          exprx=out$E
          full=fit[,c("amp","a","b","R2", "genes","qval", "pval")]
          ss=SS[[jj]]
          sg=rownames(ss)[which(ss$model.c %in% mds)]
          full=full[sg,]
          full=cbind(full, ss[sg, idt])
          colnames(full)[ncol(full)]="amp_C"
          #if (val=="R") {full=subset(full, qval<qcut)}
          full$R=2*sqrt(full$a^2+full$b^2)
          fulll=rbind(fulll, full)
        }
        full=fulll
        if(nrow(full)>1){
          full$kind=names(OUT)[j]
          all=rbind(all,full)
          breaks=seq(0, max(max(full[,val])+0.1,1), by=0.01)
          if(val %in% c("pval","qval")) breaks=breaks^7
          amp.cut = cut(full[,val], breaks, right=FALSE)
          freq.cut = table(amp.cut)
          if(val %in% c("R", "amp_C")) freq.cut = rev(freq.cut) #also tbz=tibble(R=rev(breaks[-1]), n.genes=cumsum.frq, kind=names(OUT)[j])
          cumsum.frq=c(cumsum(freq.cut),nrow(full[,val]))+1
          tbz=tibble(R=breaks[-1], n.genes=cumsum.frq, kind=l)
          if(val%in% c("R", "amp_C")) tbz=tibble(R=rev(breaks[-1]), n.genes=cumsum.frq, kind=l)
          tbt=rbind(tbt,tbz)
        }
      }
      if(!is.null(tbt)){
        if(val %in% c("R", "amp_C")) print(ggplot(tbt)+geom_line(aes(x=R,y=n.genes, color=kind, linetype=kind), size=th)+
                                             scale_y_log10()+labs(title=i,y="# of genes", x="log2(peak to trough)", color="", linetype="")+
                                             theme_minimal()+scale_color_manual(values=c("#008ECC","#111E6C","#1034A6","darkred"))+
                                             theme(text = element_text(size=sz)))#+theme(legend.position = "none"))
        else print(ggplot(tbt)+geom_line(aes(x=R,y=n.genes, color=kind, linetype=kind), size=th)+
                     scale_y_log10()+scale_x_continuous(trans=reverselog_trans(10))+
                     labs(title=i,y="# of genes", x=val)+theme_minimal()+
                     scale_color_manual(values=c("#008ECC","#111E6C","#1034A6","darkred"))+
                     theme(text = element_text(size=sz)))#+theme(legend.position = "none"))
      }
    }
  }
  for(tx in nn){
    tbt=NULL
    for(j in c(paste(tx,pox[1], sep="-"),paste(tx,pox[2], sep="-"))){
      out=OUT[[j]]
      fit=out$data.fit
      jj=gsub("-OLD", "",gsub("-YOUNG", "",gsub("-FEMALE", "",gsub("-MALE", "", j))))
      md=l=gsub("^.*-", "", j)
      if(strict==FALSE){
        if(md %in% c("OLD", "FEMALE")){mds=c(3,4,5)}
        else{mds=c(2,4,5)}}
      else{
        if(md %in% c("OLD", "FEMALE")){mds=c(3)}
        else{mds=c(2)}}       
      idt=paste("amp_", md, sep="")
      inf.phi=out$phi
      exprx=out$E
      full=fit[,c("amp","a","b","R2", "genes","qval", "pval")]
      ss=SS[[jj]]
      sg=rownames(ss)[which(ss$model.c %in% mds)]
      full=full[sg,]
      full=cbind(full, ss[sg, idt])
      colnames(full)[ncol(full)]="amp_C"
      #if (val=="R") {full=subset(full, qval<qcut)}
      full$R=2*sqrt(full$a^2+full$b^2)
      if(nrow(full)>1){
        full$kind=names(OUT)[j]
        all=rbind(all,full)
        breaks = seq(0, max(max(full[,val])+0.1,1), by=0.01)
        if(val %in% c("pval","qval")) breaks=breaks^7
        amp.cut = cut(full[,val], breaks, right=FALSE)
        freq.cut = table(amp.cut)
        if(val %in% c("R", "amp_C")) freq.cut = rev(freq.cut) 
        cumsum.frq=c(cumsum(freq.cut),nrow(full[,val]))+1
        tbz=tibble(R=breaks[-1], n.genes=cumsum.frq, kind=l)
        if(val%in% c("R", "amp_C")) tbz=tibble(R=rev(breaks[-1]), n.genes=cumsum.frq, kind=l)
        tbt=rbind(tbt,tbz)
      }
    }
    if(!is.null(tbt)){
      if(val %in% c("R", "amp_C")) print(ggplot(tbt)+geom_line(aes(x=R,y=n.genes, color=kind, linetype=kind), size=th)+
                                           scale_y_log10()+labs(title=tx,y="# of genes", x="log2(peak to trough)", color="", linetype="")+
                                           theme_minimal()+scale_color_manual(values=c("#008ECC","#111E6C","#1034A6","darkred"))+
                                           theme(text = element_text(size=sz)))#+theme(legend.position = "none"))
      else print(ggplot(tbt)+geom_line(aes(x=R,y=n.genes, color=kind, linetype=kind), size=th)+
                   scale_y_log10()+scale_x_continuous(trans=reverselog_trans(10))+
                   labs(title=tx,y="# of genes", x=val)+theme_minimal()+
                   scale_color_manual(values=c("#008ECC","#111E6C","#1034A6","darkred"))+
                   theme(text = element_text(size=sz)))#+theme(legend.position = "none"))
    }
  } 
  dev.off()
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

#Plots the densities according to the model selection parameters
Plot_MS_density<-function(OUT, SS, phi, phenot, MS=T, div="MF", pth=NULL, strict=F, qcut=0.2, size=20, th=1){
  
  full_gene_phi=NULL
  if(is.null(pth)){
    if(div=="MF") pdf("./Figure2/Fig2_B-E-H2.pdf")
    if(div=="AGE")pdf("./Figure3/Fig3_C-E-H2.pdf")
  }else{
    pdf(pth)
  }
  
  phenot$SUBJID=gsub("^.*-","",phenot$SUBJID)
  phenot$age_cat="middle"
  phenot$sex=phenot$SEX
  phenot$sex[phenot$SEX==1]="MALE"
  phenot$sex[phenot$SEX==2]="FEMALE"
  phenot$age_cat[phenot$AGE>60]="OLD"
  phenot$age_cat[phenot$AGE<50]="YOUNG"
  phio=phi
  if(div=="MF"){
    phi.dff=NULL
    for(i in c(1:2)){
      phi=phio[phenot$SUBJID[phenot$sex==unique(phenot$sex)[i]]]
      phi=phi[!is.na(phi)]
      phi.df=tibble(phi=rep((1:1000-1)/500*pi,1))
      phi.df$hour=phi.df$phi/pi*12
      phi.df$count=phi.df$phi
      kapp=15
      for(s in 1:length(phi.df$phi)){
        phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(phi))))
      }
      phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2]
      phi.df$divs=unique(phenot$sex)[i]
      phi.dff=rbind(phi.dff, phi.df)
    }
    ppt=ggplot(phi.dff)+geom_line(aes(x=hour,y=dens, colour=divs), size=th)+coord_polar()+
      scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ylim(0,NA)+
      theme_minimal()+labs(x="Time of day [h]", y="Density", title="Full", colour= "")+scale_color_manual(values=c("#008ECC","#111E6C","#1034A6","darkred"))+
      theme(text = element_text(size=sz))
    print(ppt)
  }
  if(div=="AGE"){
    phi.dff=NULL
    for(i in c(1:2)){
      phi=phio[phenot$SUBJID[phenot$age_cat==c("YOUNG", "OLD")[i]]]
      phi=phi[!is.na(phi)]
      phi.df=tibble(phi=rep((1:1000-1)/500*pi,1))
      phi.df$hour=phi.df$phi/pi*12
      phi.df$count=phi.df$phi
      kapp=20
      for(s in 1:length(phi.df$phi)){
        phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(phi))))
      }
      phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2]
      phi.df$divs=c("YOUNG", "OLD")[i]
      phi.dff=rbind(phi.dff, phi.df)
    }
    ppt=ggplot(phi.dff)+geom_line(aes(x=hour,y=dens, colour=divs), size=th)+coord_polar()+
      scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ylim(0,NA)+
      theme_minimal()+labs(x="Time of day [h]", y="Density", title="Full", colour= "")+scale_color_manual(values=c("#008ECC","#111E6C","#1034A6","darkred"))+
      theme(text = element_text(size=sz))
    print(ppt)
  }
  phi.full=NULL
  phi.dfff=NULL
  gen.names=gsub("-OLD", "",gsub("-YOUNG", "",gsub("-FEMALE", "",gsub("-MALE", "", names(SS)))))
  for(name in unique(gen.names)){
    idx=which(startsWith(names(OUT), name))
    phi.dff=NULL
    ss=SS[[name]]
    for (id in idx) {
      md=l=gsub("^.*-", "", names(OUT)[id])
      if(strict==FALSE){
        if(md %in% c("OLD", "FEMALE")){mds=c(3,4,5)}
        else{mds=c(2,4,5)}}
      if(strict==TRUE){
        if(md %in% c("OLD", "FEMALE")){mds=c(3)}
        else{mds=c(2)}}  
      out=OUT[[id]]
      phi=out$phi
      df=out$data.fit
      sg=rownames(ss)[which(ss$model.c %in% mds)]
      if (MS) df=df[sg,]
      if(!MS) {df=subset(df,qval<qcut & amp>Rcut)}
      gene_phi=df$phase/12*pi
      if (MS) {
        ssg=ss[sg,]
        lnt=ncol(ss)
        if(md %in% c("OLD", "FEMALE")){gene_phi=ssg[,lnt-7]/12*pi}
        else{gene_phi=ssg[,lnt-11]/12*pi}
      }
      phi.df=tibble(phi=rep((1:1000-1)/500*pi,1))
      phi.df$hour=phi.df$phi/pi*12
      phi.df$count=phi.df$phi
      kapp=20
      for(s in 1:length(phi.df$phi)){
        phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(phi))))
      }
      phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2]
      phi.df$tissue=name
      phi.df$divs=gsub("^.*-", "", names(OUT)[id])
      phi.df$kind="Donors"
      s.phi.df=phi.df
      phi.df=tibble(phi=rep((1:1000-1)/500*pi,1))
      phi.df$hour=phi.df$phi/pi*12
      phi.df$count=phi.df$phi
      kapp=20
      for(s in 1:length(phi.df$phi)){
        phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(gene_phi))))
      }
      phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2]
      phi.df$tissue=name
      phi.df$divs=gsub("^.*-", "", names(OUT)[id])
      phi.df$kind="Genes"
      phi.df=rbind(phi.df, s.phi.df)
      phi.dff=rbind(phi.dff, phi.df)
      phi.dfff=rbind(phi.dfff, tibble(genes=gene_phi, div=tolower(gsub("^.*-", "", names(OUT)[id]))))
    }
    ppt=ggplot(phi.dff)+geom_line(aes(x=hour,y=dens, linetype=kind, colour=divs), size=th)+
      scale_color_manual(values= c("#008ecc", "#111e6c"))+coord_polar()+scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+
      ylim(0,NA)+theme_minimal()+labs(x="Time of day [h]", y="Density", title=name, colour="", linetype= "")+theme(text = element_text(size=sz))
    print(ppt)
    phi.full=rbind(phi.full, phi.dff)
  }
  if(div=="MF"){
    phi.dff=NULL
    for(i in c(1:2)){
      phi=phio[phenot$SUBJID[phenot$sex==unique(phenot$sex)[i]]]
      phi=phi[!is.na(phi)]
      phi.df=tibble(phi=rep((1:1000-1)/500*pi,1))
      phi.df$hour=phi.df$phi/pi*12
      phi.df$count=phi.df$phi
      kapp=20
      for(s in 1:length(phi.df$phi)){
        phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(phi))))
      }
      phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2]
      phi.df$divs=unique(phenot$sex)[i]
      phi.df$kind="Donors"
      phi.dff=rbind(phi.dff, phi.df)
      
      gene_phi=subset(phi.dfff, div==unique(phenot$sex)[i])
      phi.df=tibble(phi=rep((1:1000-1)/500*pi,1))
      phi.df$hour=phi.df$phi/pi*12
      phi.df$count=phi.df$phi
      kapp=20
      for(s in 1:length(phi.df$phi)){
        phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(gene_phi$genes))))
      }
      phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2]
      phi.df$divs=unique(phenot$sex)[i]
      phi.df$kind="Genes"
      phi.dff=rbind(phi.dff, phi.df)
      
      
    }
    ppt=ggplot(phi.dff)+geom_line(aes(x=hour,y=dens, linetype=kind, colour=divs), size=th)+
      scale_color_manual(values= c("#008ecc", "#111e6c"))+coord_polar()+scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+
      ylim(0,NA)+theme_minimal()+labs(x="Time of day [h]", y="Density", title="Full", colour="", linetype= "")+theme(text = element_text(size=sz))
    print(ppt)
  }
  if(div=="AGE"){
    phi.dff=NULL
    for(i in c(1:2)){
      phi=phio[phenot$SUBJID[phenot$age_cat==c("YOUNG", "OLD")[i]]]
      phi=phi[!is.na(phi)]
      phi.df=tibble(phi=rep((1:1000-1)/500*pi,1))
      phi.df$hour=phi.df$phi/pi*12
      phi.df$count=phi.df$phi
      kapp=20
      for(s in 1:length(phi.df$phi)){
        phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(phi))))
      }
      phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2]
      phi.df$divs=c("YOUNG", "OLD")[i]
      phi.df$kind="Donors"
      phi.dff=rbind(phi.dff, phi.df)
      
      gene_phi=subset(phi.dfff, div==c("YOUNG", "OLD")[i])
      phi.df=tibble(phi=rep((1:1000-1)/500*pi,1))
      phi.df$hour=phi.df$phi/pi*12
      phi.df$count=phi.df$phi
      kapp=20
      for(s in 1:length(phi.df$phi)){
        phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(gene_phi$genes))))
      }
      phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2]
      phi.df$divs=phi.df$divs=c("YOUNG", "OLD")[i]
      phi.df$kind="Genes"
      phi.dff=rbind(phi.dff, phi.df)
    }
    ppt=ggplot(phi.dff)+geom_line(aes(x=hour,y=dens, linetype=kind, colour=divs), size=th)+
      scale_color_manual(values= c("#008ecc", "#111e6c"))+coord_polar()+scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+
      ylim(0,NA)+theme_minimal()+labs(x="Time of day [h]", y="Density", title="Full", colour="", linetype= "")+theme(text = element_text(size=sz))
    print(ppt)
  }
  ppt=ggplot(phi.full)+geom_line(aes(x=hour,y=dens,colour=tissue,linetype=divs))+
    coord_polar()+scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ylim(0,NA)+theme_minimal()+
    labs(x="time of day", y="density", title="full")+scale_color_manual(values=colorandum[unique(phi.full$tissue)])+
    theme(legend.position = "none",text = element_text(size=sz))
  print(ppt)
  dev.off()
}

#Plots the cumulative barplot and the relative split between models
Plot_MS_barplot<-function(SS, div, pth=NULL){
  if(is.null(pth)){
    if(div=="MF") pth= "./Figure2/Fig2_D.pdf"
    if(div=="AGE") pth= "./Figure3/Fig3_D.pdf"
  }
  model_gn=tibble()
  totz=tibble()
  for(j in names(SS)){
    tmg=tibble(tissue=j,genes=0, scaled=0, model=0)
    ss=SS[[j]]
    tot=0
    for (mds in c(2:5)) {
      sg=rownames(ss)[which(ss$model.c %in% mds)]
      tot=tot+length(sg)
    }
    for (mds in c(2:5)) {
      sg=rownames(ss)[which(ss$model.c %in% mds)]
      tmg[2]=length(sg)
      tmg[3]=tmg[2]/tot
      tmg[4]=mds
      model_gn=rbind(model_gn, tmg)
    }
    tut=tibble(tissue=j,genes=tot)
    totz=rbind(totz, tut)
  }
  hjst=1
  ang=90
  model_gn$model=factor(model_gn$model, c(5,4,3,2))
  model2=model_gn[model_gn$model==2,]
  totuz=totz$tissue[base::order(totz$genes)]
  totz$tissue=factor(totz$tissue,totuz)
  model_gn$tissue=factor(model_gn$tissue,totuz)
  filz=c(colroma$hex[c(50,100,150,200)])
  p1=ggplot(model_gn, aes(x=tissue, y=scaled,fill = as.factor(model)))+geom_bar(stat="identity")+
    scale_fill_manual(values=filz,name="Model")+theme_void()+
    theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=15))+
    labs(x="Tissue", y="Relative fraction of genes in each model")
  
  ggsave(filename= gsub(".pdf", "1.pdf", pth),p1,width = 12, height = 10)    
  
  p2=ggplot(totz, aes(x=tissue, y=genes))+geom_bar(stat="identity", fill = colroma$hex[250])+
    theme_minimal()+theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=15),panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
    labs(x="Tissue", y="total number of genes in models 2 to 5")+theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  
  ggsave(filename= gsub(".pdf", "2.pdf", pth),p2,width = 12, height = 10)     
}

### Main ###
if(!file.exists(file="./paper_data/raw/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")){ 
  samp = fread('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
  fwrite(samp, file = "./paper_data/raw/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
}else{
  samp = fread("./paper_data/raw/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
}
samp$sub.id=spliti(samp$SUBJID,"-",2)
samp$AGE=as.numeric(spliti(samp$AGE,"-",1))+5

if(!exists("N.cores")) N.cores = 18 
if(!exists("as.paper")) as.paper=FALSE
if(!exists("use.paper.DIP")) use.paper.DIP=FALSE

if(use.paper.DIP)   phi=get(load("./paper_data/DIPs.RData")) else   phi=get(load("./data/DIPs.RData"))

if(!as.paper){
  OUT.MF=get(load("./data/OUT/OUT_MF.RData"))
  OUT.age=get(load("./data/OUT/OUT_AGE.RData"))
  OUT.all=get(load("./data/OUT/OUT_ALL.RData"))
  SS.age=get(load("./data/OUT/SS_AGE.RData"))
  SS.MF=get(load("./data/OUT/SS_MF.RData"))
  
}

if(as.paper){
  OUT.MF=get(load("./paper_data/OUT_paper/OUT_MF.RData"))
  OUT.age=get(load("./paper_data/OUT_paper/OUT_AGE.RData"))
  OUT.all=get(load("./paper_data/OUT_paper/OUT_ALL.RData"))
  SS.age=get(load("./paper_data/OUT_paper/SS_AGE.RData"))
  SS.MF=get(load("./paper_data/OUT_paper/SS_MF.RData"))
  
}


colroma=vroom("./paper_data/roma.txt",  col_names = FALSE, show_col_types = FALSE)
colroma$hex=paste("#", as.hexmode(round(colroma$X1*255,0)),as.hexmode(round(colroma$X2*255,0)),as.hexmode(round(colroma$X3*255,0)), sep="")
full_col=vroom("./paper_data/GO_full-colorandum.csv", show_col_types = FALSE)
full_col=full_col[full_col$Class!="Cells",]
dec_names=full_col$`Short name`
names(dec_names)=full_col$`Full name`
colorandum=full_col$`# color`
names(colorandum)=full_col$`Full name`
nmz=unique(full_col$Class)

dir.create("./plot/Figure2", showWarnings = FALSE)
dir.create("./plot/Figure3", showWarnings = FALSE)

phenot=get(load("./paper_data/phenotypes.RData"))

###################### model selection plots #########################


#Parameters for the subsequent plots


#MS: use model selection
if(!exists("MS")) MS=T
#Plot parameters
if(!exists("sz")) sz=20
if(!exists("th")) th=1
#q-value cut. Note that any qcut>0.2 has no bearing if MS==T
if(!exists("qcut")) qcut=0.2
#Parameter to determine if using genes only rhythmic in condition X or genes also rhythmic in condition X
if(!exists("strict")) strict=F
#Value for the cumulative plots, possible values:
#"R": amplitude
#"pval": p-value
#"qval": q-value

if(!exists("val")) val="R"

for (div in c("MF", "AGE")){
  if(div=="MF") {
    OUT= OUT.MF
    SS=SS.MF
    pthb="./plot/Figure2/Fig2_D.pdf"
    pthc="./plot/Figure2/Fig2_C-F-I.pdf"
    pthd="./plot/Figure2/Fig2_B-E-H.pdf"
    
  }
  if(div=="AGE"){
    OUT= OUT.age
    SS=SS.age
    pthb="./plot/Figure3/Fig3_D.pdf"
    pthc="./plot/Figure3/Fig3_B-F-I.pdf"
    pthd="./plot/Figure3/Fig3_C-E-H.pdf"
  } 
  Plot_MS_barplot(SS=SS, div=div, pth=pthb)
  
  Plot_MS_cumulative(OUT=OUT, SS=SS, MS=MS, div=div, val=val, pth=pthc, strict=strict, qcut=qcut, size=size, th=th)
  
  Plot_MS_density(OUT=OUT, SS=SS, phi=phi, phenot=phenot, MS=MS, div=div, pth=pthd, strict=strict, qcut=qcut, size=size, th=th)
  
}
####### cSVD ######

lb=6
pt=3
sz=18
gene_inf=get(load("./paper_data/CGRs.RData"))
for (div in c("MF", "AGE")){
  if(div=="MF") {
    OUT= OUT.MF
    Plot_cSVD(OUT, gene_inf, full_col,loc ="./plot/Figure2/Fig2_A" , CT=15, dot_size =pt, label_size = lb, text_size = sz)
  }
  if(div=="AGE") {
    OUT= OUT.age
    Plot_cSVD(OUT, gene_inf, full_col,loc ="./plot/Figure3/Fig3_A" , CT=15, dot_size =pt, label_size = lb, text_size = sz)
  }
}
