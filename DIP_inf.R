rm(list=ls())
gc()
####### RUN AFTER data_preparation

source("supplementary_functions.R")
source("CHIRAL.R")

dir.create(file.path("./data"), showWarnings = FALSE)

as.paper=TRUE

samp <- fread('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
samp$sub.id=spliti(samp$SUBJID,"-",2)
samp$AGE=as.numeric(spliti(samp$AGE,"-",1))+5
CPM.all.norm.large=get(load("./paper_data/CPM.all.norm_large.RData"))
CPM.all.norm=get(load("./paper_data/CPM.all.norm.RData"))

E=CPM_to_E(CPM.all.norm)

gene_inf=get(load("./paper_data/CGRs.RData"))

OUT= mclapply(E,infer_l, gene_inf, mc.cores=16)

OUT=Fit_OUT(OUT)

OUT=Set_OUT(OUT)

donors=Unique_donors(OUT)

phi_matrix=Create_phi_matrix(OUT, donors)

phi_paper=get(load("./paper_data/DIPs.RData"))

phi=Phi.from.phi_mat(phi_matrix, ct=1.95) #these phases will not be identical to the DIP for the various existing stochastic steps

save(phi, file="./data/DIPs.RData")

qplot(phi_paper,phi)

E=CPM_to_E(CPM.all.norm.large)

if(as.paper){phi=phi_paper}

OUT=Make_big_OUT(E, phi)

OUT=Fit_OUT(OUT)

save(OUT, file="./data/OUT_all.RData")


