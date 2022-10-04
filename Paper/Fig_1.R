rm(list=ls())
gc()
source("supplementary_functions.R")
source("CHIRAL.R")

####### RUN AFTER from_cpm #######

colroma=vroom("./paper_data/roma.txt",  col_names = FALSE)
colroma$hex=paste("#", as.hexmode(round(colroma$X1*255,0)),as.hexmode(round(colroma$X2*255,0)),as.hexmode(round(colroma$X3*255,0)), sep="")
full_col=vroom("./paper_data/GO_full-colorandum.csv", show_col_types = FALSE)

as.paper=TRUE

if(as.paper){
  load("./paper_data/OUT_paper/OUT_ALL.RData")
  load("./paper_data/DIPs.RData")
  colroma=vroom("./paper_data/roma.txt",  col_names = FALSE, show_col_types = FALSE)
}
if(!as.paper){
  load("./data/OUT/OUT_ALL.RData")
  load("./data/DIPs.RData")
}

dir.create("./Figure1", showWarnings = FALSE)

pdf("./Figure1/Fig1_F.pdf", sep="")
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

{
  pdf("./Figure1/Fig1_G.pdf")
  phi.full=NULL
  pho=NULL
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
  
  ppt=ggplot(phi.full)+geom_line(aes(x=hour,y=dens,colour=tissue), size=th)+coord_polar()+scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ylim(0,NA)+
    theme_minimal()+labs(x="time of day", y="density", title="All Tissues")+scale_color_manual(values=colorandum[names(OUT)])+
    theme(legend.position = "none", text = element_text(size=sz))
  print(ppt)
  dev.off()
}



val="R"
all=NULL
tbt=NULL


pdf("./Figure1/Fig1_H.pdf")
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
                         theme(text = element_text(size=sz), legend.position = "none"))#+theme(legend.position = "none"))
    else print(ggplot(tbt)+geom_line(aes(x=R,y=n.genes, color=kind, linetype=kind), size=th)+
                 scale_y_log10()+scale_x_continuous(trans=reverselog_trans(10))+labs(title=i,y="# of genes", x=val)+
                 theme_minimal()+scale_color_manual(values=colorandum[unique(tbt$kind)])+theme(text = element_text(size=sz), legend.position = "top"))#+theme(legend.position = "none"))
  }
}
dev.off()


lb=6
pt=3
sz=18
gene_inf=get(load("./paper_data/CGRs.RData"))
HSF1_genes=get(load("./paper_data/HSF1_g.RData"))
{
  Plot_cSVD(OUT, gene_inf, full_col,loc ="./Figure1/Fig1_B" , CT=15, dot_size =pt, label_size = lb, text_size = sz)
  Plot_cSVD(OUT, HSF1_genes, full_col,loc ="./Figure1/Fig1_D" , CT=5, dot_size =pt, label_size = lb, text_size = sz, ymax=0.6)
}

{
  pdf("./Figure1/Fig1_C.pdf")
  tissuex=c("Brain - Cortex", "Adipose - Visceral (Omentum)")
  geni=gene_inf[c(2,8)]
  Plot_profiles(OUT, tissuex, geni, val="qval")
  dev.off()
}

{
  pdf("./Figure1/Fig1_E.pdf")
  tissuex=c("Brain - Amygdala", "Spleen")
  geni=c("HSPH1", "HSPA1B")
  Plot_profiles(OUT, tissuex, geni, val="qval")
  dev.off()
}
