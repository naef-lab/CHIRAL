############# Load packages  ###################
library(colorRamps)
library(circlize)
library(viridis)
library(enrichR)
library(ComplexHeatmap)
library(lmtest)
library(scico)
library(parallel)
library(combinat)
library(plotrix)
library(rWikiPathways)
library(RColorBrewer)
library(RCy3)

################################################


################ Source external functions #####
setwd("/home/cgobet/CHIRAL/Paper/")
source("ridgeInR.R")
source("nconds_functions.R")
source("nconds.R")
################################################

################ Define functions ##############

# Moving average with a 4-hour window and 1-hour step
return_ma=function(x,t){
  
  r=rep(NA,24)
  vv=c(23:24,1:24,1:2)
  for(i in 1:24){
    if(vv[i+4] > vv[i]){
      r[i]=mean(as.numeric(x[which(t > vv[i] & t <= vv[i+4])]),na.rm=T)
      
    }else{
      r[i]=mean(as.numeric(x[which(t > vv[i] | t <= vv[i+4])]),na.rm=T)
      
    }
    
  }
  return(r)
}

# Extract gene names from GO term table returned from enrichR
gene2table=function(y){
  if(nrow(y) > 1 ){
    s <- sort(unique(unlist(strsplit(y$Genes, ';'))))
    m <- sapply(s, grepl, x=y$Genes, fixed=TRUE)
    m[] <- (as.integer(m))
    
    m=as.data.frame(t(m))
    colnames(m)=y$Term
    m
  }
  
}

# Harmonic regression with likelihood ratio test and p-value. Compute amplitude and phase.
harm_reg=function(x, t, period)
{
  n=length(x)
  fit0=lm(x~1)
  c=cos(2*pi*t/period)
  s=sin(2*pi*t/period)
  fit1=lm(x~c+s)
  a=coef(fit1)[2]
  b=coef(fit1)[3]
  p.val=lrtest(fit1, fit0)$Pr[2]
  amp=2*sqrt(a^2+b^2)
  phase=atan2(b,a)%%(2*pi)
  phase=period*phase/(2*pi)
  
  c(p.val_harm=p.val,phase_harm=phase,amp_harm=amp,period_arm=period)
}

################################################


# Main function to generate the complex heatmap including GO term and TF activity

generate_comp_heatmap=function(tiss.2, dat, dat.c1, dat.c2, subsampling, mo, output_path, cytoscape = FALSE){
  #tiss.2 = tissue
  #dat = the full tissue dataset (OUT files from model_selection.R (./paper_data/OUT_paper/OUT_MF.RData or ./data/OUT/OUT_MF.RData))
  #dat.c1 = subsampled tissue data and model selection for MALE or YOUNG (from ./paper_data/OUT_paper/OUT_paper/SS_AGE.RData or ./data/OUT/SS_MF.RData)
  #dat.c2 = subsampled tissue data and model selection for FEMALE or OLD (from ./paper_data/OUT_paper/SS_AGE.RData or ./data/OUT/SS_MF.RData)
  #subsampling = MF or AGE
  #mo = rhythmicity model (1-5)
  #output_path = path to save the heatmaps pdfs
  #cytoscape = TRUE or FALSE, allowing to plot the enriched wikipathways with amplitudes depicted in both conditions. Cytoscape is required.
  
  v.1=NULL
  
  if(subsampling == 'MF'){
    var.amp.1="amp_MALE"
    var.amp.2="amp_FEMALE"
    var.phase.1="phase_MALE"
    var.phase.2="phase_FEMALE"
    c.1="MALE"
    c.2="FEMALE"
  }else{
    var.amp.1="amp_YOUNG"
    var.amp.2="amp_OLD"
    var.phase.1="phase_YOUNG"
    var.phase.2="phase_OLD"
    c.1="YOUNG"
    c.2="OLD"
  }
  
  dat.mo =subset(dat, model.c == mo)
  dat.c1.mo = dat.c1$E[rownames(dat.mo),]
  dat.c2.mo = dat.c2$E[rownames(dat.mo),]
  
  if(nrow(dat.mo) >10 ){
    
    vv=order(dat.mo[[var.phase.1]],dat.mo[[var.phase.2]])
    dat.mo=dat.mo[vv,]
    dat.c1.mo=dat.c1.mo[vv,]
    dat.c2.mo=dat.c2.mo[vv,]
    
    time.c1=24*dat.c1$phi/(2*pi)
    ii=order(time.c1)
    dat.c1.mo=dat.c1.mo[,ii]
    time.c1=time.c1[ii]
    dat.ma.c1=t(apply(dat.c1.mo,1,return_ma,time.c1))
    
    time.c2=24*dat.c2$phi/(2*pi)
    ii=order(time.c2)
    dat.c2.mo=dat.c2.mo[,ii]
    time.c2=time.c2[ii]
    dat.ma.c2=t(apply(dat.c2.mo,1,return_ma,time.c2))
    
    
    dat.ma.c1=sweep(dat.ma.c1,1,rowMeans(dat.ma.c1,na.rm=T),FUN="-")
    dat.ma.c2=sweep(dat.ma.c2,1,rowMeans(dat.ma.c2,na.rm=T),FUN="-")
    
    TT=cbind(dat.ma.c1, dat.ma.c2)
    ref.c1=rep(0,24)
    names(ref.c1)=0:23
    ref.c2=rep(0,24)
    names(ref.c2)=0:23
    
    c1.phase=hist(dat.mo[[var.phase.1]][dat.mo[[var.phase.1]]!=0],breaks=seq(0,24,1))$count
    c2.phase=hist(dat.mo[[var.phase.2]][dat.mo[[var.phase.2]]!=0],breaks=seq(0,24,1))$count
    
    TT.anno=rownames(dat.mo)[order(-dat.mo[[var.amp.1]],-dat.mo[[var.amp.2]])][1:round(0.025*nrow(dat.mo))]
    
    ##### GO TERM #######
    dbs <- listEnrichrDbs()
    
    dbb=dbs$libraryName[grep('GO_Biological_Process_2021|KEGG_2021_Human|WikiPathway_2021_Human',dbs$libraryName)]
    GO.sub=list()
    
    g.name=gsub(".+_","",rownames(dat.mo))
    
    GO=enrichr(g.name,databases=dbb)
    GO.sub=lapply(GO,function(x) subset(x,Adjusted.P.value < 0.05 & Combined.Score > 20))
    
    Go.sub.m=lapply(GO.sub,gene2table)
    
    go.term=sort(unique(unlist(lapply(Go.sub.m,colnames))))
    GG=matrix(0,ncol=length(go.term),nrow=length(g.name))
    colnames(GG)=go.term
    rownames(GG)=g.name
    
    for(kk in 1:length(Go.sub.m)){
      if(!is.null(Go.sub.m[[kk]])){
        Go.sub.m[[kk]]=Go.sub.m[[kk]][rownames(Go.sub.m[[kk]])%in%rownames(GG),]
        GG[rownames(Go.sub.m[[kk]]),colnames(Go.sub.m[[kk]])]= as.matrix(Go.sub.m[[kk]])
      }
    }
    
    duplicated.columns <- duplicated(t(GG))
    GG<- GG[, !duplicated.columns]
    
    p.val=lapply(GO.sub,function(x) x[,c('Term','Adjusted.P.value')])
    p.val=do.call(rbind,p.val)
    p.val=-log10(p.val[match(colnames(GG),p.val[,1]),2])
    
    
    if(length(GG)!=0 & !is.null(nrow(GG))){
      colnames(GG)=sapply(strsplit(colnames(GG),split='\\('),"[[",1)
      
      GG=GG[,order(-p.val)]
      p.val=p.val[order(-p.val)]
      if(ncol(GG) > 15){
        GG=GG[,1:15]
        p.val=p.val[1:15]
      }else{
        GG=GG[,1:ncol(GG)]
        p.val=p.val[1:length(p.val)]
        
      }
    }else{
      if(is.null(nrow(GG))){
        
        GG=matrix(sample(length(GG)))
        colnames(GG)="No Go Term"
        p.val=1
      }else{
        GG=matrix(sample(nrow(GG)))
        colnames(GG)="No Go Term"
        p.val=1
      }
    }
    
    if(cytoscape){
      
      WPs=GO.sub$WikiPathway_2021_Human
      WPs_2=paste0("WP",gsub(".+WP","",WPs$Term))
      k=1
      if(!is.null(WPs_2)){ 
        for(i in WPs_2){ 
          commandsRun(paste0('wikipathways import-as-pathway id=',i)) 
          dat.s=subset(dat,model.c!=0)
          delta=dat.s[[var.amp.1]] - dat.s[[var.amp.2]]
          
          DD=data.frame(ensembl=gsub("\\..+","",rownames(dat.s)),gene_symbol=sub(".+_","",rownames(dat.s)),delta, amp1=dat.s[[var.amp.1]], amp2=dat.s[[var.amp.2]] )
          toggleGraphicsDetails()
          
          loadTableData(DD, data.key.column = "ensembl", table.key.column = "Ensembl") 
          
          # allN=unlist(getAllNodes())
          # allN=allN[allN!=""]
          # 
          # allE=unlist(getAllNodes())
          # allE=allE[allE!=""]
          # 
          # setNodeBorderColorBypass(allN, "#FFFFFF")
          # setNodeLabelColorBypass(allN, "#000000")
          # setNodeFontSizeBypass(allN,16)
          
          setNodeCustomBarChart(c("amp1","amp2"),style.name='WikiPathways',colors=c("#67A9CF","#EF8A62"))
          skip_to_next <- FALSE
          k=k+1
          tryCatch(exportImage(paste0(output_path,subsampling,"_",tiss.2, WPs_2[k]), type='PDF'), error = function(e) { skip_to_next <<- TRUE})
          
          if(skip_to_next) { next }     
          
          deleteNetwork()
        }
      }
    }
    
    ha.go.1 = HeatmapAnnotation(
      dist1 = anno_barplot(
        p.val,
        bar_width = 1,
        gp = gpar(col = "white", fill = scico(30, palette = 'vik')[30]),
        border = FALSE,
        baseline=0), show_annotation_name = FALSE)
    
    h.go=Heatmap((GG),name="Go",col=c('white','black'),show_column_names = T,top_annotation  = ha.go.1, cluster_rows=FALSE, cluster_columns = F,
                 row_names_side = "left", column_names_gp = gpar(fontsize = 5),
                 heatmap_legend_param = list(title = "GO Term",
                                             title_gp = gpar(fontsize = 5),labels_gp = gpar(fontsize = 5),
                                             legend_height = unit(2, "cm")))
    ##### TF Enrichment ########
    
    site=get(load("./paper_data/TF_hg19_average_max1_gs.RData"))
    N=as.matrix(site)
    
    dat.c1.all=dat.c1$E[rownames(dat),]
    dat.c2.all=dat.c2$E[rownames(dat),]
    
    time.c1.all= 24*dat.c1$phi/(2*pi)
    time.c2.all= 24*dat.c2$phi/(2*pi)
    
    i.c1=order(time.c1.all)
    time.c1.all=time.c1.all[i.c1]
    dat.c1.all=dat.c1.all[,i.c1]
    
    i.c2=order(time.c2.all)
    time.c2.all=time.c2.all[i.c2]
    dat.c2.all=dat.c2.all[,i.c2]
    
    E= cbind(dat.c1.all,dat.c2.all)
    E[is.na(E)]=0
    
    rownames(E)=gsub("\\..+_.+","",rownames(E))
    common.genes = intersect(rownames(E), rownames(N))
    
    E = E[common.genes, ]
    N = N[common.genes, ]
    
    E = as.matrix(E)
    N = as.matrix(N)
    
    opt = optimize.lambda(N, E)
    TF.ana = ridge.regression(N, E, opt$lambda.opt)
    
    ##### Add chip_results
    
    TF.ana.hat=TF.ana$Ahat
    
    TF.model=nconds(TF.ana.hat,conds=c(rep(c.1,ncol(dat.c1.all)),rep(c.2,ncol(dat.c2.all))),t=c(time.c1.all, time.c2.all),
                    out.prefix = NULL, N.cores = N.cores)
    
    TF.model.sub=subset(TF.model,model== mo &  AICW > 0.5)
    TF.model.sub$zscore=TF.ana$combined.Zscore[rownames(TF.model.sub)]
    
    top.TF=unique(c(rownames(TF.model.sub[order(-TF.model.sub$zscore)[1:10],]),
                    rownames(TF.model.sub[order(-TF.model.sub[[var.amp.1]])[1:10],])))
    if(mo==3){
      
      top.TF=unique(c(rownames(TF.model.sub[order(-TF.model.sub$zscore)[1:10],]),
                      rownames(TF.model.sub[order(-TF.model.sub[[var.amp.2]])[1:10],])))
    }
    
    TF.ana.hat=TF.ana.hat[order(TF.model[[var.phase.1]],TF.model[[var.phase.2]]),]
    TF.ana.hat=TF.ana.hat[rownames(TF.ana.hat)%in%top.TF,,drop=F]
    
    TF.c1=TF.ana.hat[,1:ncol(dat.c1.all),drop=F]
    TF.ma.c1=t(apply(TF.c1, 1, return_ma, time.c1.all))
    
    TF.c2=TF.ana.hat[,(ncol(dat.c1.all)+1):(ncol(TF.ana.hat)),drop=F]
    TF.ma.c2=t(apply(TF.c2, 1, return_ma, time.c2.all))
    
    TF.ana.me=cbind(TF.ma.c1,TF.ma.c2)
    
    if(ncol(TF.ana.me) == 0){TF.ana.me=t(as.matrix(rep(1,48))+ runif(48))}
    
    ha.TF.1 = HeatmapAnnotation(Motif= t(TF.ana.me),
                                annotation_name_gp=gpar(fontsize=5),
                                col = list(Motif=colorRamp2(c(range(TF.ana.me,na.rm=T)[1], 0, range(TF.ana.me,na.rm=T)[2]),
                                                            c('darkblue','white','darkred'))),
                                height = unit(3, "cm"), simple_anno_size_adjust = TRUE, annotation_legend_param = list(title = "TF activity", title_gp = gpar(fontsize = 5),labels_gp = gpar(fontsize = 5),
                                                                                                                       legend_height = unit(2, "cm")))
    
    
    ha = rowAnnotation(foo = anno_mark(at = match(TT.anno,rownames(TT)), labels = gsub(".+_","",TT.anno), labels_gp=gpar(fontsize=5)))
    
    ha1 = HeatmapAnnotation(name="phase",
                            dist1 = anno_barplot(
                              c(c1.phase,c2.phase),
                              bar_width = 1,
                              gp = gpar(col = "white", fill = scico(30, palette = 'vik')[30]),
                              border = FALSE,
                              baseline=0), show_annotation_name = FALSE)
    
    ha3 = HeatmapAnnotation(summary = anno_summary(height = unit(3, "cm"),width = unit(0.1, "cm"), pch = 1,size=unit(0.05, "cm")))
    
    
    pdf(paste0(output_path,subsampling,"_",tiss.2,"_",mo,".pdf"))
    
    mm=quantile(abs(TT),0.99,na.rm=T)
    
    v=Heatmap(TT, name = "Expression",
              col =circlize::colorRamp2(seq(-mm, mm, length = 50), scico(50,palette='vik')),
              cluster_rows=FALSE,
              cluster_columns = FALSE,
              show_row_dend = FALSE, show_column_dend = FALSE,
              show_row_names = FALSE, show_column_names = TRUE,
              column_split = factor(c(rep(c.1,24),rep(c.2,24)),levels=c(c.1,c.2)),
              row_title_gp = gpar(col = "#FFFFFF00"),
              column_labels =paste0('ZT',rep(0:23,2)),
              column_names_gp = gpar(fontsize = 6),
              row_gap = unit(5, "mm"),
              column_gap = unit(5, "mm"),
              top_annotation = ha1, bottom_annotation = ha.TF.1,heatmap_legend_param = list(title = "Expression",
                                                                                            title_gp = gpar(fontsize = 5),labels_gp = gpar(fontsize = 5),
                                                                                            legend_height = unit(2, "cm")))
    #                right_annotation = ha,
    
    if(sum(dat.mo[[var.amp.1]])!=0){
      v=v +  Heatmap(dat.mo[[var.amp.1]], name = paste0("Amp. ",c.1), width = unit(5, "mm"),
                     heatmap_legend_param = list(title = paste0("Amp. ",c.1,"(log2)"),
                                                 title_gp = gpar(fontsize = 5),labels_gp = gpar(fontsize = 5),
                                                 legend_height = unit(2, "cm")),show_row_names = FALSE,
                     top_annotation =ha3, column_names_gp = gpar(fontsize = 6),col = colorRamp2(c(0.5,max(dat.mo[[var.amp.1]])),c('white',"#601200"))) }
    
    if(sum(dat.mo[[var.amp.2]])!=0){
      v=v + Heatmap(dat.mo[[var.amp.2]], name =  paste0("Amp. ",c.2), width = unit(5, "mm"),
                    heatmap_legend_param = list(title = paste0("Amp. ",c.2,"(log2)"),title_gp = gpar(fontsize = 5),labels_gp = gpar(fontsize = 5),
                                                legend_height = unit(2, "cm")),show_row_names = FALSE,
                    top_annotation =ha3, column_names_gp = gpar(fontsize = 6),col = colorRamp2(c(0.5,max(dat.mo[[var.amp.2]])),c('white',"#601200")))
      
    }
    
    
    v.1=Heatmap(TT, name = paste("Expression",mo),
                col =circlize::colorRamp2(seq(-mm, mm, length = 50), scico(50,palette='vik')),
                cluster_rows=FALSE,
                cluster_columns = FALSE,
                show_row_dend = FALSE, show_column_dend = FALSE,
                show_row_names = FALSE, show_column_names = TRUE,
                column_split = factor(c(rep(c.1,24),rep(c.2,24)),levels=c(c.1,c.2)),
                row_title_gp = gpar(col = "#FFFFFF00"),
                column_labels =paste0('ZT',rep(0:23,2)),
                column_names_gp = gpar(fontsize = 6),
                row_gap = unit(5, "mm"),
                column_gap = unit(5, "mm"))#,right_annotation = ha)
    
    v= v + h.go
    draw(v,column_title =paste0(tiss.2," / ",nrow(dat.mo)," genes /"," Model: ", mo))
    
    
    for(i in 1){ decorate_heatmap_body("Go", {
      grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))
    },slice=i)}
    
  }else{
    pdf(paste0(output_path,subsampling,"_",tiss.2,"_",mo,".pdf"))
    
    grid.rect(gp=gpar(col="white"))
  }
  
  dev.off()
  
  return(v.1)
}

####### Run the heatmaps for every model in AGE or MF condition

dv='MF' # or 'AGE'
output_path="./plot/complex_heatmaps/"
N.cores = 18 
dir.create(output_path, showWarnings = FALSE,recursive = TRUE)
as.paper=FALSE

if(as.paper){ 
  path.dat="./paper_data/OUT_paper/"
}else{
  path.dat="./data/OUT/"
  
}
  if(dv =='MF'){
    dat.ss=get(load(paste0(path.dat,"SS_MF.RData")))
    dat.all=get(load(paste0(path.dat, "OUT_MF.RData")))
    c.1='MALE'
    c.2='FEMALE'
  }else{
    dat.ss=get(load(paste0(path.dat, "SS_AGE.RData")))
    dat.all=get(load(paste0(path.dat, "OUT_AGE.RData")))
    c.1='YOUNG'
    c.2='OLD'
  }


for(k in names(dat.ss)){
  DD=NULL
  
  for(mo in 2:5){
    
    print(k);print(mo);
    outi=generate_comp_heatmap(k, dat.ss[[k]], dat.all[[paste(k,c.1,sep="-")]],
                               dat.all[[paste(k,c.2,sep="-")]], dv, mo, output_path, FALSE)
    if(!is.null(outi)){
      DD=DD %v% generate_comp_heatmap(k, dat.ss[[k]], dat.all[[paste(k,c.1,sep="-")]],
                                      dat.all[[paste(k,c.2,sep="-")]], dv, mo, output_path, FALSE)
    }
  }
  
  if(!is.null(DD)){
    pdf(paste0(output_path,dv,"_",k,"_summary.pdf"))
    draw(DD,row_title=paste(dv,k), show_heatmap_legend = FALSE)
    dev.off()
  }
}






