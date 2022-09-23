`CHIRAL` (Circular HIerarchical Reconstruction ALgorithm) is an R package that provides a Bayesian method to infer the circular coordinates of a set of samples. It has been developed for the circadian clock but can be used on any other circualr dataset, e.g. the cell cycle.

## Getting Started

These instructions will allow you to get `CHIRAL` running on your machine. 

### Prerequisites
You need to install R (see https://www.r-project.org/).

`CHIRAL` accepts a matrix containing on the rows the variable quantity (in general genes) and on the column different measures taken at different times in the cycle of interest (samples along the day for the circadian clock). Data that were used to implement `CHIRAL` are typically produced by bulk RNA-Seq, but it has been tested also on scRNA-seq. The variable quantities should follow a sinusoidal curve along the cycle, thus for RNA-seq and scRNA-seq log-transformed (and possibly smoothed in the case of sc) data are required to obtain relevant results. To use this on RNA-seq data we suggest either to only use exon count or to  input exon and intron measurements as two different genes. It come with helper function to calculate relevant parameters for circualr distributions, break the intrinsic symmetries of the algorythm and improve visual plotting.
Given the fully probabilistic and unsupervised approach to the problem `CHIRAL` has two symmetries: there is a rotational symmetry and a sign symmettry; these two can be broke either with additional knowledge on the process of interest or by comparison with the real pahses. 

### Installing

To install `CHIRAL` run the following code in R.
```
install.packages("devtools")
devtools::install_github("naef-lab/CHIRAL")
```
## Quick start
### Example dataset 
`CHIRAL` comes with example data in form of a list called example; these data are taken from "Transcriptomic analyses reveal rhythmic and CLOCK-driven pathways in human skeletal muscle", Perrin et al. eLife 2018. The list contains exon count data: example[["Muscle_exon"]], a vector with the real time when these samples were collected: example[["true_phi"]], and a vector indicating the subject from which each sample was taken simData[["time"]].

### Running an example
```
require(Chiral)

data<-example[["Muscle_exon"]]            #Import the data. This matrix is formatted as it should be: on the columns we have samples, on the rows we have genes, rownames are the gene names (in this case ENSG)

gene_inf=example[["CRG_ens"]]             #We select genes related to the clock. Note that we can also run CHIRAL on the full matrix but might pick up other sources of variation and intepret those as periodic.

out=CHIRAL(data, clockgenes = gene_inf)   #We run the algorythm using the genes we just selected
```
#Outputs

We can look at all the output
```
names(out)
```
we should find inferred phases in [0,2*pi], inferred standard deviation of the data, gene parameters gene weigths (relevant with two state model),
iteration of stop, inferred standard deviation of data of rhythmic genes, the input matrix, history of the Q function of the EM, genes selected, matrix used for inference.
A shorter output should coincide with an error message.
The main result of CHIRAL is the phase ordering, we will focus on that.
```
hist(out$phi)                                                           #We can start by plotting an histogram and look at the general distribution

true.phi=example[["true_phi"]]                                          #Load the real phases

plot(true.phi, out$phi)                                                 #Do the plot to compare inferred and true phases
```
#Helper functions

Now the plotting might look bad and might change from run to run due to the intrinsic nature of the algorythm.

For this reason we provide helper functions to break invariances and use the cyclic nature of the phases vor visually better plots.

```
inf.phi=delta.phi(true.phi, out$phi, mode = "say", median_scale=12/pi)  #Adjust the phases and get the median absolute error (in hours) using the helper function

abs(cor.c(true.phi, inf.phi))                                           #Get the correlation

adj.phi=adjust.phases(true.phi, inf.phi)                                #Adjust the phases for plotting

plot(true.phi, inf.phi)                                                 #Do the plot to compare inferred and true phases

```


## Help
A documentation using `?CHIRAL` is available as well as for the other helper functions. 
