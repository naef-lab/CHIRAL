source("./supplementary_functions.R")


gtex = vroom("/scratch/For_cedric_with_pizza/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads_2.gct")
GTEx = fread("https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")
gtex=as.data.frame(GTEx)
names(gtex)=gsub("\\.","-", names(gtex))
samp <- fread('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')

gtex.tissue=samp[match(colnames(gtex),samp$SAMPID),'SMTSD']

tiss=unique(gtex.tissue)
tiss=as.matrix(tiss[-1])
gtex.tissue=as.matrix(gtex.tissue)

CPM.all=list()
for(k in tiss){ 
  
  posx=which(gtex.tissue==k)
  gtex.sub=as.matrix(gtex[,posx])
  rownames(gtex.sub)=paste(gtex$Name,gtex$Description,sep="_")
  gtex.sub=subset(gtex.sub,rowMeans(gtex.sub)>10)
  Ex.y <- DGEList(counts=gtex.sub)
  Ex.y <- calcNormFactors(Ex.y)
  CPM = data.frame(cpm(Ex.y,log=TRUE, prior.count = 0.25))
  CPM.all[[k]]=CPM
}


dir.create(file.path("./data"), showWarnings = FALSE)

dir.create(file.path("./data/CPM"), showWarnings = FALSE)

save(CPM.all, file = "./data/CPM/CPM_full.RData")

CPM.all.norm=Norm.CPM(CPM.all, high_filter=T, ncores=6)

save(CPM.all.norm,file="./data/CPM/CPM.all.norm.RData")

CPM.all.norm_large=Norm.CPM(CPM.all, high_filter=F, ncores=18)

save(CPM.all.norm_large,file="./data/CPM/CPM.all.norm_large.RData")
