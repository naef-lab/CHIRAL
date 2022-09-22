samp <- fread('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
samp$sub.id=spliti(samp$SUBJID,"-",2)
samp$AGE=as.numeric(spliti(samp$AGE,"-",1))+5

as.paper=TRUE

if(!as.paper){
  
  phio=phi=get(load("./data/DIP.RData"))
  
  E=CPM_to_E(CPM.all.norm.large)
  
  E=split_E_sex(E, samp)
  
  OUT.MF=Make_big_OUT(E, phi)
  
  OUT.MF=Fit_OUT(OUT.MF)
  
  E=split_E_age(E, samp)
  
  OUT.age=Make_big_OUT(E, phi)
  
  OUT.age=Fit_OUT(OUT.age)
}

if(as.paper){
  OUT.MF=get(load("./paper_data/OUT_paper/OUT_MF.RData"))
  OUT.age=get(load("./paper_data/OUT_paper/OUT_age_n5.RData"))
  OUT.all=get(load("./paper_data/OUT_paper/OUT_all.RData"))
  phio=phi=get(load("./paper_data/DIP.RData"))
}


SS.age=get(load("./paper_data/OUT_paper/SS_age_n5.RData"))
SS.MF=get(load("./paper_data/OUT_paper/SS_MF.RData"))

colroma=vroom("./paper_data/roma.txt",  col_names = FALSE, show_col_types = FALSE)
full_col=vroom("./paper_data/GO_full-colorandum.csv", show_col_types = FALSE)
full_col=full_col[full_col$Class!="Cells",]
dec_names=full_col$`Short name`
names(dec_names)=full_col$`Full name`
colorandum=full_col$`# color`
names(colorandum)=full_col$`Full name`
nmz=unique(full_col$Class)

dir.create("./Figure2", showWarnings = FALSE)
dir.create("./Figure3", showWarnings = FALSE)



###################### Cumulative model selection #########################


MS=MSS=T
sz=20
th=1
qcut=0.2
strict=F
for (strict in c(FALSE)) {
  for (div in c("MF", "age_n5")){
    qcut=2
    if(div=="MF") OUT= OUT.MF
    if(div=="age_n5") OUT= OUT.age
    pox=unique(gsub("^.*-", "", names(OUT)))
    nn=names(OUT)[duplicated(gsub("-old", "",gsub("-young", "",gsub("-Female", "",gsub("-Male", "", names(OUT))))))]
    nn=gsub("-old", "",gsub("-young", "",gsub("-Female", "",gsub("-Male", "", nn))))
    OUT=OUT[c(paste(nn,pox[1], sep="-"),paste(nn,pox[2], sep="-"))]
    for(val in c("R")){
      if(div=="MF") pdf("./Figure2/Fig2_C-F-I.pdf")
      if(div=="age_n5")pdf("./Figure3/Fig3_B-F-I.pdf")
      
      all=NULL
      tbt=NULL
      
      id=list(div1=which(endsWith(names(OUT),pox[1])),div2=which( endsWith(names(OUT),pox[2])))
      for (ik in names(id) ) {
        fulll=NULL
        for(j in names(OUT)[id[[ik]]]){
          out=OUT[[j]]
          fit=out$data.fit
          jj=gsub("-old", "",gsub("-young", "",gsub("-Female", "",gsub("-Male", "", j))))
          if(div=="MF"){SS=SS.MF}
          else{SS=SS.age}
          md=l=gsub("^.*-", "", j)
          if(strict==FALSE){
            if(md %in% c("old", "Female")){mds=c(3,4,5)}
            else{mds=c(2,4,5)}}
          else{
            if(md %in% c("old", "Female")){mds=c(3)}
            else{mds=c(2)}}
          if(md=="Male") md="ma"
          if(md=="Female") md="fe"
          idt=paste("amp_", toupper(md), sep="")
          
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
              jj=gsub("-old", "",gsub("-young", "",gsub("-Female", "",gsub("-Male", "", j))))
              if(div=="MF"){SS=SS.MF}
              else{SS=SS.age}
              md=gsub("^.*-", "", j)
              if(strict==FALSE){
                if(md %in% c("old", "Female")){mds=c(3,4,5)}
                else{mds=c(2,4,5)}}
              else{
                if(md %in% c("old", "Female")){mds=c(3)}
                else{mds=c(2)}}
              #interesting.genes=inter.genes
              if(md=="Male") md="ma"
              if(md=="Female") md="fe"
              idt=paste("amp_", toupper(md), sep="")
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
          jj=gsub("-old", "",gsub("-young", "",gsub("-Female", "",gsub("-Male", "", j))))
          if(div=="MF"){SS=SS.MF}
          else{SS=SS.age}
          md=l=gsub("^.*-", "", j)
          if(strict==FALSE){
            if(md %in% c("old", "Female")){mds=c(3,4,5)}
            else{mds=c(2,4,5)}}
          else{
            if(md %in% c("old", "Female")){mds=c(3)}
            else{mds=c(2)}}        #interesting.genes=inter.genes
          if(md=="Male") md="ma"
          if(md=="Female") md="fe"
          idt=paste("amp_", toupper(md), sep="")
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
            if(val %in% c("R", "amp_C")) freq.cut = rev(freq.cut) #also tbz=tibble(R=rev(breaks[-1]), n.genes=cumsum.frq, kind=names(OUT)[j])
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
  }
}


############### Density MS ############

for (div in c("MF", "age_n5")) {
  full_gene_phi=NULL
  if(div=="MF") pdf("./Figure2/Fig2_B-E-H.pdf")
  if(div=="age_n5")pdf("./Figure3/Fig3_C-E-H.pdf")
  if(div=="MF") OUT= OUT.MF
  if(div=="age_n5") OUT= OUT.age
  phenot=get(load("./paper_data/phenotypes.RData"))
  phenot$SUBJID=gsub("^.*-","",phenot$SUBJID)
  phenot$age_cat="middle"
  phenot$sex=phenot$SEX
  phenot$sex[phenot$SEX==1]="male"
  phenot$sex[phenot$SEX==2]="female"
  phenot$age_cat[phenot$AGE>60]="old"
  phenot$age_cat[phenot$AGE<50]="young"
  phio=get(load("./data/DIPs.RData"))
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
        #if(s<length(phi.df$phi)/2)
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
  if(div=="age_n5"){
    phi.dff=NULL
    for(i in c(1:2)){
      phi=phio[phenot$SUBJID[phenot$age_cat==c("young", "old")[i]]]
      phi=phi[!is.na(phi)]
      phi.df=tibble(phi=rep((1:1000-1)/500*pi,1))
      phi.df$hour=phi.df$phi/pi*12
      phi.df$count=phi.df$phi
      kapp=20
      for(s in 1:length(phi.df$phi)){
        #if(s<length(phi.df$phi)/2)
        phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(phi))))
      }
      phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2]
      phi.df$divs=c("young", "old")[i]
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
  gen.names=gsub("-old", "",gsub("-young", "",gsub("-Female", "",gsub("-Male", "", names(OUT)))))
  if(div=="MF"){SS=SS.MF}
  else{SS=SS.age}
  for(name in unique(gen.names)){
    idx=which(startsWith(names(OUT), name))
    phi.dff=NULL
    ss=SS[[name]]
    for (id in idx) {
      md=l=gsub("^.*-", "", names(OUT)[id])
      if(strict==FALSE){
        if(md %in% c("old", "Female")){mds=c(3,4,5)}
        else{mds=c(2,4,5)}}
      if(strict==TRUE){
        if(md %in% c("old", "Female")){mds=c(3)}
        else{mds=c(2)}}  
      out=OUT[[id]]
      phi=out$phi
      df=out$data.fit
      sg=rownames(ss)[which(ss$model.c %in% mds)]
      if (MS) df=df[sg,]
      if(!MS) df=df[df$qval<qcut,]
      #df=df[df$qval<qcut,]
      gene_phi=df$phase/12*pi
      if (MSS) {
        ssg=ss[sg,]
        lnt=ncol(ss)
        if(md %in% c("old", "Female")){gene_phi=ssg[,lnt-7]/12*pi}
        else{gene_phi=ssg[,lnt-11]/12*pi}
      }
      phi.df=tibble(phi=rep((1:1000-1)/500*pi,1))
      phi.df$hour=phi.df$phi/pi*12
      phi.df$count=phi.df$phi
      kapp=20
      for(s in 1:length(phi.df$phi)){
        #if(s<length(phi.df$phi)/2)
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
        #if(s<length(phi.df$phi)/2)
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
        #if(s<length(phi.df$phi)/2)
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
        #if(s<length(phi.df$phi)/2)
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
  if(div=="age_n5"){
    phi.dff=NULL
    for(i in c(1:2)){
      phi=phio[phenot$SUBJID[phenot$age_cat==c("young", "old")[i]]]
      phi=phi[!is.na(phi)]
      phi.df=tibble(phi=rep((1:1000-1)/500*pi,1))
      phi.df$hour=phi.df$phi/pi*12
      phi.df$count=phi.df$phi
      kapp=20
      for(s in 1:length(phi.df$phi)){
        #if(s<length(phi.df$phi)/2)
        phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(phi))))
      }
      phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2]
      phi.df$divs=c("young", "old")[i]
      phi.df$kind="Donors"
      phi.dff=rbind(phi.dff, phi.df)
      
      gene_phi=subset(phi.dfff, div==c("young", "old")[i])
      phi.df=tibble(phi=rep((1:1000-1)/500*pi,1))
      phi.df$hour=phi.df$phi/pi*12
      phi.df$count=phi.df$phi
      kapp=20
      for(s in 1:length(phi.df$phi)){
        #if(s<length(phi.df$phi)/2)
        phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(gene_phi$genes))))
      }
      phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2]
      phi.df$divs=phi.df$divs=c("young", "old")[i]
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


######### Barplot MS ########

for (div in c("MF", "age_n5")) {
  if(div== "MF") SS=SS.MF
  if(div=="age_n5") SS=SS.age
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
      #tmg=tibble(tissue=j,gene=sg, model=mds)
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
  #p1=
  filz=c(colroma$hex[c(50,100,150,200)])
  p1=ggplot(model_gn, aes(x=tissue, y=scaled,fill = as.factor(model)))+geom_bar(stat="identity")+
    scale_fill_manual(values=filz,name="Model")+theme_void()+
    theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=15))+
    labs(x="Tissue", y="Relative fraction of genes in each model")
  
  if(div=="MF") ggsave( filename= "./Figures2-3/Fig2_D1.pdf",p1,width = 12, height = 10)
  if(div=="age_n5") ggsave( filename= "./Figures2-3/Fig2_D1.pdf",p1,width = 12, height = 10)    
  
  p2=ggplot(totz, aes(x=tissue, y=genes))+geom_bar(stat="identity", fill = colroma$hex[250])+
    theme_minimal()+theme(axis.text.x = element_text(angle = ang, hjust=hjst),text = element_text(size=15),panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
    labs(x="Tissue", y="total number of genes in models 2 to 5")+theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  
  if(div=="MF") ggsave( filename= "./Figures2-3/Fig2_D2.pdf",p2,width = 12, height = 10)
  if(div=="age_n5") ggsave( filename= "./Figures2-3/Fig2_D2.pdf",p2,width = 12, height = 10) 
}

lb=6
pt=3
sz=18
gene_WP=WP.ls[[nm]]
for (div in c("MF", "age_n5")){
  if(div=="MF") {
    OUT= OUT.MF
    Plot_cSVD(OUT, gene_inf, full_col,loc =file.path(wd, "Figure2/Fig2_A") , CT=15, dot_size =pt, label_size = lb, text_size = sz)
  }
  if(div=="age_n5") {
    OUT= OUT.age
    Plot_cSVD(OUT, gene_inf, full_col,loc =file.path(wd, "Figure3/Fig3_A") , CT=15, dot_size =pt, label_size = lb, text_size = sz)
  }
}
