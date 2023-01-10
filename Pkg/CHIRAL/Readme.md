`CHIRAL` (Circular HIerarchical Reconstruction ALgorithm) is an R package that provides a Bayesian method to infer the circular coordinates of a set of samples. It has been developed for the circadian clock but can be used on any other circular dataset, e.g. the cell cycle.

## Getting Started

These instructions will allow you to get `CHIRAL` running on your machine. 

### Prerequisites
You need to install R (see https://www.r-project.org/).

`CHIRAL` accepts a matrix containing on the rows the variable quantities (in general genes) and on the column the different measures taken at different times in the cycle of interest (i.e. samples along the day for the circadian clock). Data that were used to implement `CHIRAL` are typically produced by bulk RNA-Seq but has been tested also on scRNA-seq. The variable quantities should follow a sinusoidal curve along the cycle, thus for RNA-seq and scRNA-seq log-transformed and normalized (and possibly smoothed in the case of scRNA-seq) data are required to obtain relevant results. The package comes with helper functions to calculate relevant parameters for circular distributions, break the intrinsic symmetry of the algorithm, and visual plotting.
Given the fully probabilistic and unsupervised approach to the problem `CHIRAL` has two symmetries: there is a rotational symmetry and a sign symmetry; these two can be broken either with additional knownledge on the process of interest or by comparison with the real phases. 

### Installing

To install `CHIRAL` run the following code in R.
```
install.packages("devtools")
devtools::install_github("naef-lab/CHIRAL/Pkg/CHIRAL")
```
## Quick start
### Example dataset 
`CHIRAL` comes with an example dataset taken from "Transcriptomic analyses reveal rhythmic and CLOCK-driven pathways in human skeletal muscle", Perrin et al. eLife 2018. The data include time-series count RNA-seq from human biopsies `muscle_exon`, a vector with the real sample collection time `true_phi`, and the Ensembl clock reference genes used for the sample ordering `CRG_ens`.

### Running an example

The matrix is formatted as required: on the columns we have the samples, on the rows the genes, and the rownames are the gene names (in this case Ensembl format).
```
require(CHIRAL)
head(muscle_expm)        
```
We selected genes related to the clock as references. Note that, `CHIRAL` can also be run on the full matrix but might pick up other sources of variation and interpret those as periodic.
```
head(CRG_ens)  
```
Run the algorithm with the selected genes. 
```
out<-CHIRAL(data, clockgenes = CRG_ens)   
```
##### Outputs

The output contains multiple variables:
```
names(out)
```
In particular: the inferred phases in [0,2*pi], the inferred standard deviation of the data, the gene parameters gene weigths (relevant with two state model),
the stop iteration number, the inferred standard deviation of rhythmic gene, the input matrix, the history of the EM Q function, the selected genes, the matrix used for inference.
A shorter output should coincide with an error message. The main result of CHIRAL is the phase ordering, then we will focus on that:

Plotting an histogram and look at the general phase distribution.
```
hist(out$phi)    
```
The true phases:
```
head(true_phi)  
```
Plot to compare inferred and true phases.

```
plot(true_phi, out$phi)                                                 
```

##### Helper functions

Due to the intrinsic nature of the algorithm, the correlation might not seem correct in the plot. We provide helper functions to break invariances and use the cyclic nature of the phases for better visualization.

Adjust the phases and get the median absolute error (in hours) using the helper function. Get the circular correlation between real and inferred phases.
```
inf_phi=delta.phi(true_phi, out$phi, mode = "say", median_scale=12/pi)  
abs(cor.c(true_phi, inf_phi))       
```
Adjust the phases and plot the inferred and true phases.

```
adj_phi=adjust.phases(true_phi, inf_phi)                                
plot(true_phi, adj_phi)                                                
```


## Help
A documentation using `?CHIRAL` is available as well as for the other helper functions. 
