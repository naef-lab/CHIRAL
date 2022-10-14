if(!("BiocManager" %in% installed.packages()[,"Package"])) install.packages("BiocManager")
list.of.packages <- c("edgeR", "data.table", "parallel", "R.utils", "lmtest", "ggplot2", "vroom", "tibble", "gridExtra", "ggrepel", "combinat",
                      "colorRamps", "circlize", "viridis", "enrichR","ComplexHeatmap", "scico", "rWikiPathways", "RColorBrewer", "RCy3")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages, force=TRUE)
if(!("CHIRAL" %in% installed.packages()[,"Package"])) devtools::install_github("/naef-lab/CHIRAL/tree/master/Pkg/CHIRAL")
library(edgeR)
library(data.table)
library(parallel)
library(R.utils)
library(CHIRAL)
library(lmtest)
library(ggplot2)
library(vroom)
library(tibble)
library(gridExtra)
library(ggrepel)
library(combinat)
library(colorRamps)
library(circlize)
library(viridis)
library(enrichR)
library(ComplexHeatmap)
library(lmtest)
library(scico)
library(parallel)
library(combinat)
library(rWikiPathways)
library(RColorBrewer)
library(RCy3)




#### General variables ####
#Define them to overwrite the default values

your_path=file.path(getwd(), "Paper") # should be the path were you cloned the repository, inside the Paper folder
N.cores = 18  # Number of core to paralellize the different pre-processing functions. Be sure to have enough RAM if you increase the number of cores used.
setwd(your_path)
as.paper=FALSE #If TRUE will not recompute files and used the ones you downloaded 
use.paper.DIP=TRUE #If TRUE the script recreates exactly the paper figures as it removes the stochasticity of the TIP inference

### Normalize data if not planning on using provided data ###
 
if(! as.paper) source("To_CPM.R")

### Infer and save the DIPs if not planning on using provided ones ###

source("DIP_inf.R")

### Generate Figure 1 ###

source("Fig_1.R")

### Model selection to assign gene rhythmicity ###

if(!as.paper) source("Model_selection.R") 

### Generate Figures 2 and 3 
#Parameters for Figures 2 and 3
#MS: use model selection
MS=T
#Plot parameters
sz=20
th=1
#q-value cutoff. Note that any qcut>0.2 has no bearing if MS==T
qcut=0.2
#Parameter to determine if using genes only rhythmic in condition X or genes also rhythmic in condition X
strict=F
#Value for the cumulative plots, possible values:
#"R": amplitude
#"pval": p-value
#"qval": q-value
val="R"

#Generate panels that are not heatmaps

source("Fig_2-3.R")

#Generate heatmaps

source("Complex_heatmaps.R")
