# Sex-dimorphic and age-dependent organization of 24 hour gene expression rhythms in human
In order to reproduce the figures and analysis of the paper follow the instruction following.

### Scripts and packages
Clone this github repository

Install the `CHIRAL` package from https://github.com/naef-lab/CHIRAL/tree/master/Pkg/CHIRAL

### Data
Download the external data to use pre-computed files
Visit https://doi.org/10.5281/zenodo.6637875 and 
Download the folders "OUT_paper" and "CPM". 
Copy that folder into the "Paper/paper_data" folder 

## Running the code

You should run the scripts in the order:

To_CPM.R

DIP_inf.R

Model_selection.R

Fig_1.R

Fig_2-3.R

Complex_heatmaps.R

## Reproducibility

There is in the code the as.paper variable set to TRUE by default.
This variable  controls if you load the pre computed data or use the ones you generate.
Note that there is both stochasticity and possible numerical errors, so if this parameter is set to FALSE the figures and findings might differ, although very very slightly, from what is shown in the paper.

### Helper functions

Please also find all the functions used to obtain the results.

# Licence 

The eintierty of these files and scripts are under the CC BY-NC license.

If you are unsure how to use products under creative commons license please visit https://creativecommons.org/licenses/



