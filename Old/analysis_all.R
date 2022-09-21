#General analysis for all together
E.matrix=list()
for(name in names(CPM.all.norm)){
  E=list()
  E$E=CPM.all.norm[[name]]
  E$E=E$E[,!is.na(E$E[1,])]
  if(!is.null(E$E) && any(gsub("\\|.*$","",gsub("^.*_", "",rownames(E$E)))=="PER3")){
    E$tissue=name
    E$type="GTEx"
    gene="DBP"
    sampz=gsub("\\..*$","", gsub("GTEX.","", colnames(E$E)))
    cmn=intersect(names(phi), sampz)
    idx=match(cmn,sampz)
    if(length(idx)>48){
      E$E=E$E[,idx]
      E$phi=phi[cmn]
      E$score=score[cmn]
      E$geni=gene_inf
      E.matrix[[name]]=E
    }
  }
}
OUT=E.matrix


for(i in names(OUT)){ 
  out=OUT[[i]]
  E=out$E
  #E=E[,-ncol(E)]
  #E=subset(E,rowMeans(E)>0)
  common.genes=intersect(common.genes,rownames(E))
  phase=out$phi
  #phase=E.matrix[[i]]$phi
  # cat(i, "correlation", abs(cor.c(E.matrix[[i]]$phi, OUT[[i]]$phi)), "")
  # inf.phi=delta.phi(E.matrix[[i]]$phi, OUT[[i]]$phi, mode = "say")
  #print(qplot(E.matrix[[i]]$phi,phase))
  dat.fit=as.data.frame(t(apply(E,1,f24_R2_cycling,t=24*as.numeric(phase)/(2*pi))))
  dat.fit$qval=p.adjust(dat.fit$pval)
  out$data.fit=data.frame(dat.fit,E,genes=rownames(E))
  if(is.null(out$type)){out$type="mice"}
  out=order.out.setgene(out)
  OUT[[i]]=out
}
sout=OUT
SVD=svd.from.out(OUT,gene_inf)
{pdf("big_full_SVD_all.pdf")
  for(i in 1:ncol(SVD$u)){
    #i=1
    print(qplot(Arg(SVD$u[,i])%%(2*pi)*12/pi,Mod(SVD$u[,i]))+geom_label_repel(aes(label=rownames(SVD$u)))+coord_polar()+ylim(0,NA)+
            scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ labs(title=paste("Genes SVD, mode", i,"variance explained",round(SVD$d[i]/sum(SVD$d),4)),x =paste("Argument of component", i), y = paste("Modulus of component", i)))
    print(qplot(Arg(SVD$v[,i])%%(2*pi)*12/pi,Mod(SVD$v[,i]))+coord_polar()+aes(ymin=0, xmin=0,xmax=2*pi)+geom_label_repel(aes(label=rownames(SVD$v)))+ylim(0,NA)+
            scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ labs(title=paste("Tissues SVD, mode", i, "variance explained",round(SVD$d[i]/sum(SVD$d),4)),x =paste("Argument of component", i), y = paste("Modulus of component", i)))
    print(qplot(Arg(SVD$v[,i])%%(2*pi)*12/pi,Mod(SVD$v[,i]))+coord_polar()+aes(ymin=0, xmin=0,xmax=2*pi)+ylim(0,NA)+
            scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ labs(title=paste("Tissues SVD, mode", i,"variance explained",round(SVD$d[i]/sum(SVD$d),4)),x =paste("Argument of component", i), y = paste("Modulus of component", i)))
  }
  dev.off()}

inter.genes=gene_inf
common.genes=gsub("\\|.*$","",gsub("^.*_", "",rownames(E.matrix[[1]]$E)))
nmz=names(OUT)
gene.c=matrix(0, nrow = length(inter.genes),ncol=length(nmz))
gene.big=matrix(0, nrow = length(unique(common.genes)),ncol=length(nmz))
dimnames(gene.c)=list(inter.genes,nmz)
dimnames(gene.big)=list(common.genes,nmz)
types=nmz
names(types)=nmz
pdf("big_mean_all.pdf")
for(i in nmz){ 
  out=OUT[[i]]
  fit=out$data.fit
  interesting.genes=inter.genes
  inf.phi=out$phi
  exprx=out$E
  gene.list=gsub("\\|.*$","",gsub("^.*_", "",fit$genes))
  clock.coord=sapply(out$geni,function(x){match(x, gene.list)})
  clock.coord=clock.coord[!is.na(clock.coord)]
  full=fit[,c("amp","a","b","R2", "genes")]
  full$genz=gsub("\\|.*$","",gsub("^.*_", "",full$genes))
  full$R=sqrt(full$a^2+full$b^2)
  GE=fit[,c("a","b", "genes")]
  GL=GE$genes
  GE=complex(real = GE[,"a"], imaginary =  GE[,"b"])
  names(GE)=gsub("\\|.*$","",gsub("^.*_", "",GL))
  fit=fit[clock.coord,]
  A=fit[,c("mean","amp", "a","b","R2", "genes")]
  n1=which(colnames(fit)=="qval")+1
  B=as.matrix(fit[,c(n1:(ncol(fit)-1))])
  C=out$alpha[interesting.genes,]
  rownames(B)=gsub("\\|.*$","",gsub("^.*_", "",A$genes))
  B=sweep(B,1,rowMeans(B),FUN="-")
  rownames(A)=gsub("\\|.*$","",gsub("^.*_", "",A$genes))
  A=complex(real = A[,"a"], imaginary =  A[,"b"])
  names(A)=rownames(B)
  #nr=rep(-Conj(A["PER3"])/sqrt(A["PER3"]*Conj(A["PER3"])),length(A))
  #nr=rep(-Conj(A["TEF"])/sqrt(A["TEF"]*Conj(A["TEF"])),length(A))
  #A=A*nr
  alp=alpha_tot[alpha_tot$name==i,c("ensg", "genes", "use","R2", "name")]
  gene.cyc=alp[base::order(-alp$R2),]
  gene.good=gene.cyc[1:40,]
  out$pesi=Re(A)*0
  geni=names(A)
  pt=tibble(genes=names(A), dist=Mod(A), ang=(Arg(A)%%(2*pi)*12/pi), weight=out$pesi)
  #gtf=clock.integrity(ref="mice truephi", modo=20, x=out$E,  genes =  gsub("\\|.*$","",gsub("^.*_", "",rownames(out$E))), phi=out$phi)
  #gft=clock.integrity(ref="gtex infphi", modo=20, x=out$E,  genes =  gsub("\\|.*$","",gsub("^.*_", "",rownames(out$E))), phi=out$phi)
  print(ggplot(pt, aes(x=ang, y=dist, label=genes, color=weight))+coord_polar(start=0, direction=1)+
          ggrepel::geom_label_repel()+labs(x="peak phase", y="amplitude",colour = "Tissue")+ylim(0,NA)+
          geom_point()+ggtitle(paste(i, " mice ref:", round(gtf$value,4) , " n samples:", ncol(exprx) , sep=" "))+
          scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+
          scale_colour_gradient(limits=c(0,1)))
  pltlist=list()
  for(k in 1:length(geni)){
    fg=geni[k]
    idx=match(fg, gsub("\\|.*$","",gsub("^.*_", "",rownames(out$E))))
    ts=(1:1000)*pi/500
    pred=GE[fg]*Conj(complex(argument = ts))
    fit=tibble(gene=Re(pred), phase=ts,type="fit")#, conf=1)
    geneplt=tibble(gene=as.matrix(out$E)[idx,],phase=out$phi,type="data")#, conf=out$score)
    geneplt=rbind(geneplt,fit)
    pltlist[[k]]=ggplot(geneplt)+geom_point(aes(x=phase, y=gene, color=type, size=type))+labs(x="Phase",y=fg)+ theme(legend.position = "none")+scale_size_manual(values=c(1,0.00007))+scale_color_manual(values=c("green3","blue","red","darkred"))
  }
  do.call(grid.arrange,c(grobs=pltlist, top=i))
  fullR=(full[base::order(-full$R),])[1:40,]
  fullR2=(full[base::order(-full$R2),])[1:40,]
  g1=ggplot(fullR)+geom_point(aes(x=reorder(genz, -R), y=R))+theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=sz))
  g2=ggplot(fullR2)+geom_point(aes(x=reorder(genz, -R2), y=R2))+theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=sz))
  grid.arrange(g1,g2)
  geni=fullR$genz[1:12]
  pltlist=list()
  for(k in 1:length(geni)){
    fg=geni[k]
    idx=match(fg, gsub("\\|.*$","",gsub("^.*_", "",rownames(out$E))))
    ts=(1:1000)*pi/500
    pred=GE[fg]*Conj(complex(argument = ts))
    fit=tibble(gene=Re(pred), phase=ts,type="fit")#, conf=1)
    geneplt=tibble(gene=as.matrix(out$E)[idx,],phase=out$phi,type="data")#, conf=out$score)
    geneplt=rbind(geneplt,fit)
    pltlist[[k]]=ggplot(geneplt)+geom_point(aes(x=phase, y=gene, color=type, size=type))+labs(x="Phase",y=fg)+ theme(legend.position = "none")+scale_size_manual(values=c(1,0.00007))+scale_color_manual(values=c("green3","blue","red","darkred"))
  }
  do.call(grid.arrange,c(grobs=pltlist, top=i))
  
  breaks = seq(0, max(full$R), by=0.01) 
  amp.cut = cut(full$R, breaks, right=FALSE) 
  freq.cut = table(amp.cut) 
  freq.cut = rev(freq.cut) 
  cumsum.frq=c(cumsum(freq.cut),nrow(full$R))+1
  plot(rev(breaks[-1]),cumsum.frq,type='l',lwd=4,xlab='R (Log2)',ylab="Inverse Cumulative Sum",log="y")
  breaks = seq(0, max(full$R2), by=0.01) 
  amp.cut = cut(full$R2, breaks, right=FALSE) 
  freq.cut = table(amp.cut) 
  freq.cut = rev(freq.cut) 
  cumsum.frq=c(cumsum(freq.cut),nrow(full$R2))+1
  #p2=plot(rev(breaks[-1]),cumsum.frq,type='l',lwd=4,xlab='R2 (Log2)',ylab="Inverse Cumulative Sum",log="y")
  heatmap.3(cor(t(B)))
  # print(ggplot(gene.good)+geom_point(aes(x=reorder(genes, -R2), y=R2,color=use))+
  #         theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=sz))+
  #         scale_color_manual(values=c("green3","blue","red","darkred"))+
  #         ggtitle("most rythmic genes"))
  gene.c[,i]=A[interesting.genes]
  gene.big[,i]=A[common.genes]
  print(qplot(out$phi,bins=12))
  #RT=tibble(R.sq=sum(E$R2), tissue=out$tissue)
  #R.squared=rbind(R.squared,RT)
  #dev.off()
  types[i]=out$type
}
dev.off()
ptt=NULL
pdf("big_genes_all.pdf")
for(g in rownames(gene.c)){
  G=gene.c[g,]
  G=G[!is.na(G)]
  pt=tibble(tix=names(G), dist=Mod(G), ang=(Arg(G)%%(2*pi)*12/pi), gene=g)
  print(ggplot(pt, aes(x=ang, y=dist, label=tix))+coord_polar(start=0, direction=1)+
          ggrepel::geom_label_repel()+labs(x="peak phase", y="amplitude")+ylim(0,NA)+
          geom_point()+ggtitle(g)+scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24)))
  ptt=rbind(ptt, pt)
}
{g1=ggplot(subset(ptt,gene%in%c('PER1','PER2','PER3','ARNTL'))) +geom_point(aes(x=ang, y=dist, color=gene))+coord_polar(start=0, direction=1)+labs(x="peak phase", y="amplitude")+
    scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+facet_grid(~gene)+ylim(0,NA)+ theme(legend.position = "none")
  g2=ggplot(subset(ptt,gene%in%c('NPAS2','CIART','TEF','HLF','DBP')))+geom_point(aes(x=ang, y=dist, color=gene)) +coord_polar(start=0, direction=1)+labs(x="peak phase", y="amplitude")+
    scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+facet_grid(~gene)+ ylim(0,NA)+theme(legend.position = "none")
  g3=ggplot(subset(ptt,gene%in%c('CRY1','CRY2','NR1D1','NR1D2')))+geom_point(aes(x=ang, y=dist, color=gene)) +coord_polar(start=0, direction=1)+labs(x="peak phase", y="amplitude")+
    scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+facet_grid(~gene)+ ylim(0,NA)+theme(legend.position = "none")
  grid.arrange(g1,g2,g3)}
dev.off()


{sz=7
  ang=90
  hjst=1
  alpha_tot=NULL
  alpha_mean=NULL
  alpha_mid=NULL
  for (name in names(OUT)){
    ot=OUT[[name]]$data.fit
    out=OUT[[name]]
    n1=ncol(ot)-which(colnames(ot)=="qval")-1
    ot$ensg=ot$genes
    ot$genes=gsub("^.*_", "",ot$genes)
    ot$genes=gsub("\\|.*$","", ot$genes)
    #ot=ot[ot$genes %in% gene_inf,]
    ot=ot[,c("a","b", "R2", "genes","pval", "ensg")]
    ot$tissue=gsub("-.[^-]*$","", name)
    ot$type="GTEx"
    if(!is.null(OUT[[name]]$type)){ot$type=OUT[[name]]$type}
    #ot$type=names.kind$kind[names.kind$names==name]
    ot$R=sqrt(ot$a^2+ot$b^2)
    ot$WR=ot$R*ot$R2
    ot$logP=log(ot$pval)/length(OUT[[name]]$phi)
    ot$bic=-(n1*log(1-ot$R2)+2*log(n1))
    ot$bin_bic=H_theta(ot$bic, 10)
    ot$mean_BB=mean(ot$bin_bic)
    ot$use="Not used"
    ot$use[ot$genes %in% gene_inf]="Used"
    ot$div=gsub("^.*-","", name)
    ot$name=name
    #ot$W[match(out$geni,ot$genes)]=out$pesi
    alpha_tot=rbind(alpha_tot,ot)
    ot$W=0
    o=ot[ot$genes %in% out$geni,c("R2","bic","R", "W", "genes","a","b")]
    o$W=1#out$pesi
    o$WR2=o$W*o$R2
    #o=o[o$genes %in% names(full.ref),]# rememebr to load correct ref
    co=colMeans(o[,1:4])
    alpha_mean=rbind(alpha_mean,tibble(co[1],co[2],co[3],co[4], name,ot$type[1], gsub("^.*-","", name)))
    o$tissue=name
    o$div=gsub("^.*-","", name)
    alpha_mid=rbind(alpha_mid,o)
    #cat(name,"mean  binary BIC",mean(ot$bin_bic), "mean BIC",mean(ot$bic), "mean R2",mean(ot$R2), " mean log P",mean(ot$logP),"samples", ncol(E.matrix[[name]]$E), "\n")
  }
  colnames(alpha_mean)=c(colnames(o[,1:4]),"tissue", "type", "div")
}

ggplot(alpha_mean)+geom_point(aes(x=reorder(interaction(tissue,sep="-",lex.order=TRUE), -R2),y=R2,color=type))+theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=sz))+scale_color_manual(values=c("green3","blue","red","darkred"))
ggplot(alpha_mean)+geom_point(aes(x=reorder(interaction(tissue,div,sep="-",lex.order=TRUE), -R),y=R,color=type))+theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=sz))+scale_color_manual(values=c("green3","blue","red","darkred"))
ggplot(alpha_mid)+geom_boxplot(aes(x=reorder(genes, -R2, FUN = mean),y=R2))+theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=sz))+scale_fill_manual(values=c("green3","blue","red","darkred"))
ggplot(alpha_mid)+geom_boxplot(aes(x=reorder(genes, -R, FUN = mean),y=R))+theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=sz))+scale_fill_manual(values=c("green3","blue","red","darkred"))
ggplot(alpha_tot)+geom_boxplot(aes(x=reorder(interaction(tissue,div,sep="-",lex.order=TRUE), -R2, FUN = mean),y=R2,fill=type))+theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=sz))+scale_fill_manual(values=c("green3","blue","red","darkred"))
ggplot(alpha_tot)+geom_boxplot(aes(x=reorder(interaction(tissue,div,sep="-",lex.order=TRUE), logP, FUN = mean),y=logP,fill=type))+theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=sz))+scale_fill_manual(values=c("green3","blue","red","darkred"))
ggplot(alpha_tot)+geom_boxplot(aes(x=reorder(interaction(tissue,div,sep="-",lex.order=TRUE), -R, FUN = mean),y=R,fill=type))+theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=sz))+scale_fill_manual(values=c("green3","blue","red","darkred"))
ggplot(alpha_tot)+stat_density(aes(x=R,color=type),geom="line",position="identity")+scale_color_manual(values=c("green3","blue","red","darkred"))
ggplot(alpha_tot)+geom_density(aes(x=R2,color=type))+scale_color_manual(values=c("green3","blue","red","darkred"))
ggplot(alpha_tot)+geom_density(aes(x=-logP,color=type))+scale_color_manual(values=c("green3","blue","red","darkred"))


#alpha.s=alpha_tot[alpha_tot$tissue %in% dupl,]
#alp=alpha.s[alpha.s$div=="Female",c("ensg", "genes", "use","R2", "name")]
alp=alpha_tot[,c("ensg", "genes", "use","R2", "name")]
gene.cyc=spread(alp, key="name",value="R2")
gene.cyc$rem=1-rowSums(is.na(gene.cyc))/(ncol(gene.cyc)-4)
gene.cyc$mean_obj=rowMeans(gene.cyc[,4:(ncol(gene.cyc)-1)],na.rm = TRUE)
gene.cyc$min_obj=apply(gene.cyc[,4:(ncol(gene.cyc)-2)],1,gmin, 5)
gene.cyc=gene.cyc[gene.cyc$rem>0.75,]
gene.cyc=gene.cyc[base::order(-gene.cyc$mean_obj),]
gene.good=gene.cyc[1:40,]
poss.genes=gene.cyc[1:40,"genes"]
ggplot(gene.good)+geom_point(aes(x=reorder(genes, -mean_obj), y=mean_obj,color=use))+theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=sz))+scale_color_manual(values=c("green3","blue","red","darkred"))+ labs(title = "average gene R2", y="R2", x="genes")
