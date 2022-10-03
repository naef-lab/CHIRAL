#### Libraries ####
library(edgeR)
library(data.table)
library(parallel)
########

#### General variables ####

your_path="/Users/cgobet/Documents/2022/09_GTEX_revision/CHIRAL/Paper/"
N.cores = 4 # Number of core to paralellize the different pre-processing functions. Be sure to have enough RAM if you increase the number of cores used.
setwd(your_path)

########

#### Functions #####

#Filter genes with less than 10 counts on average, normalized by library size (TMM) and log2 transform with 0.25 pseudo-count 
run_edgeR = function(k){
  
  posx=which(gtex.tissue == k)
  gtex.sub=as.matrix(gtex[,posx])
  rownames(gtex.sub)=paste(gtex$Name,gtex$Description,sep="_")
  gtex.sub=subset(gtex.sub,rowMeans(gtex.sub)>10)
  Ex.y <- DGEList(counts=gtex.sub)
  Ex.y <- calcNormFactors(Ex.y)
  CPM = data.frame(cpm(Ex.y, log=TRUE, prior.count = 0.25))
  return(CPM)
  
}

#Filter the samples based on the metadata (RNA quality (RIN) >6, Mapped read > 40K, Unique mapping rate > 80%, Positive ischemic time)
# Remove covariates
Norm.CPM<- function(CPM.all, high_filter=T, ncores=18, samp, samp.2){
    samp=samp[,c('SAMPID','SMRIN',"SMTSISCH","SMMAPRT","SMUNMPRT","SMMPPD","SMNABTCH","SMATSSCR","SMTSPAX","SMPTHNTS")]
  
  if (high_filter) samp=subset(samp,SMRIN > 6 & SMMPPD > 40000000 & SMMAPRT > 0.8 &  SMTSISCH > 0 & SMATSSCR %in%c(0:4))
  else samp = subset(samp,SMRIN > 4)
  
  samp$subj.id=paste(spliti(samp$SAMPID,"-",1),spliti(samp$SAMPID,"-",2),sep="-")
  samp.all= data.frame(samp,samp.2[match(samp$subj.id,samp.2$SUBJID),])
  samp.all$SAMPID=gsub('-',"\\.",samp.all$SAMPID)
  rownames(samp.all) =samp.all$SAMPID
  samp.all$SMSISH.f= cut(samp.all$SMTSISCH, c(0,200,500,800,2000))
  
  CPM.all.norm=mclapply(CPM.all,remove_covariates,samp.all,mc.cores = ncores)
  return(CPM.all.norm)
}

# Remove covariates using a linear regression with age, sex, ischemic time and death type as covariates
remove_covariates=function(tt, samp.all){
  
  tt=subset(tt,rowMeans(tt)>0)
  tt=tt[,names(tt)%in%samp.all$SAMPID]
  tt.info=samp.all[names(tt),]
  tt.norm=tt
  
  if(ncol(tt)>10){
    for(j in 1:nrow(tt)){
      for.fit=data.frame(y=as.numeric(tt[j,]), 
                         isch=tt.info$SMSISH.f, 
                         age=tt.info$AGE, 
                         sex=as.factor(tt.info$SEX),
                         death.type=as.factor(tt.info$DTHHRDY))
      
      if(length(unique(for.fit$sex))==1){
        resid=summary(lm(formula = y ~ isch + age  + death.type, data = for.fit))$residuals
      }else{
        resid=summary(lm(formula = y ~ isch + age + sex + death.type, data = for.fit))$residuals
      }
      tt.norm[j,]=NA
      tt.norm[j,as.numeric(names(resid))]=resid
    }
    
  }else{
    tt.norm= NULL
  }
  return(tt.norm)
}

spliti= function(x, sp, nb){
  v=sapply(strsplit(x,sp),"[[",nb)
  return(v)
}
########

##### Main ####

#Read raw rna-seq count data from the GTEX database
gtex = fread("https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")
gtex=as.data.frame(gtex)
names(gtex)=gsub("\\.","-", names(gtex))

#Read metadata from the GTEX database
samp <- fread('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')

#get sample position for each tissue
gtex.tissue=samp[match(colnames(gtex),samp$SAMPID),'SMTSD']
tiss=unique(gtex.tissue)
tiss=as.matrix(tiss[-1])
gtex.tissue=as.matrix(gtex.tissue)
tiss=as.list(tiss)
names(tiss)=tiss

#Filter read count and normalize by library size, log transform, tissue-by-tissue
CPM.all= mclapply(tiss, run_edgeR, mc.cores = N.cores)

dir.create(file.path("./data"), showWarnings = FALSE)
dir.create(file.path("./data/CPM"), showWarnings = FALSE)

save(CPM.all, file = "./data/CPM/CPM_full.RData")

#Filter samples based on RNA quality, sequencing depth, mapping quality and then remove covariates (bias removal). 
# The fit residuals are used to infer TIP and DIP.
meta.1 <- fread('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')
meta.2 <- fread('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')

CPM.all.norm=Norm.CPM(CPM.all, high_filter=TRUE, ncores=N.cores, meta.1, meta.2)

save(CPM.all.norm,file="./data/CPM/CPM.all.norm.RData")

#Filter samples based on RNA quality then remove covariates (bias removal). The fit residuals are used to assess rhythmicity using DIP from above.

CPM.all.norm_large=Norm.CPM(CPM.all, high_filter=FALSE, ncores=N.cores, meta.1, meta.2)

save(CPM.all.norm_large,file="./data/CPM/CPM.all.norm_large.RData")

########