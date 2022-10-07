# Sex-dimorphic and age-dependent organization of 24 hour gene expression rhythms in human
In order to reproduce the figures and analysis of the paper, follow the instructions below.

### Repository and packages
Clone this github repository:

```
git clone https://github.com/naef-lab/CHIRAL
```
Install the `CHIRAL` package: 

```
install.packages("devtools");
devtools::install_github("naef-lab/CHIRAL/Pkg/Chiral")
```

### Data
1. Download the external data to use pre-computed files. Visit https://doi.org/10.5281/zenodo.6637875 and download the folders "OUT_paper" and "CPM". 
2. Copy those folders into "Paper/paper_data". 
3. Make sure to set your working directory such that `./paper_data/OUT_paper/` and `./paper_data/CPM/` exist and are readable.
4. Set the your R working directory in the "wherever/you/cloned/Paper" folder

## Running the code

Check that `your_path` is correctly defined. Set `N.cores`, `as.paper`, and `use.paper.DIP` to fit your needs. Run the `Main.R` file which will perform the full analysis from scratch and reproduce the main figures of the paper :

```
Rscript Main.R
```

In a more control way, you can run each script individually with the following order:

```
Rscript To_CPM.R
Rscript DIP_inf.R
Rscript Model_selection.R
Rscript Fig_1.R
Rscript Fig_2-3.R
Rscript Complex_heatmaps.R
```
## Debugging

The code has been tested on a clean installation of R. It might require a restart of the R session for the installation of some packages.

## Reproducibility

The `as.paper` variable is set to `TRUE` by default.
This variable  controls if you load the pre-computed data or use the ones that you will generate.
Note that there is stochasticity in the TIP inference. If this parameter is set to `FALSE`, the figures might slightly different from what is shown in the paper.

### Helper functions

Find all the functions used to obtain the results.

# Licence 

The entirety of these files and scripts are under the CC BY-NC license.

If you are unsure how to use products under creative commons license please visit https://creativecommons.org/licenses/



