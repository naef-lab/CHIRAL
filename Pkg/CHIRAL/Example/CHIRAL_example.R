require(CHIRAL)
data<-read.table("./paper_data/example/Muscle_exon.txt", header = TRUE, row.names = 1)
#This matrix is formatted as it should be
#on the columns we have samples
#on the rows we have genes
#rownames are the gene names (in this case ENSG)

#Let's subselect a random sample of gene that we can use to infer the phases
#note that we can also run CHIRAL on the full matrix

gene_inf=sample(rownames(data), 20)

#Let's run the algorythm using the genes we just selected

out=CHIRAL(data, clockgenes = gene_inf)

##################################################################
##########################  Outputs  #############################
##################################################################
#We can look at all the output

names(out)

#we find the inferred phases, we can plot an histogram

hist(out$phi)

#This probably look really bad due to a random gene selection

#We also have the weights of the two state models

hist(out$weights)

#That can look very different depending what happened

#We pick at random a gene that we considered rhythmic and save also its position

gn=sample(out$geni[out$weights>0.7],1)

id_gn=which(out$geni==gn)

#look at his profile

plot(out$phi,out$clock[id_gn,], xlim=c(0, 2*pi))

#We can calculate  predictions, using the alphas

alphas=out$alpha

predicted=alphas$mu+alphas$a%o%cos(out$phi)+ alphas$b%o%sin(out$phi)

#we plot the predicted vs real

plot(out$clock[id_gn,], predicted[id_gn,])

#or look at all the predictions together

plot(c(out$clock), c(predicted))

##################################################################
################  circadian phase inference ######################
##################################################################

#Now we will take the good gene to infer circadian pahse of our data

ens_inf=get(load("./Paper/paper_data/CRGs_ensg.RData"))

#Load real phases and metadata

true_phi=read.table("./paper_data/example/true_phi.txt")$x

#Do the inference

out=CHIRAL(E=data, 100,clockgenes =ens_inf,standardize = TRUE )

#Plot the results: keep in mind that there is a shift and inversion that could have happened, to solve these intrinsic problems check 

plot(true_phi, out$phi)

