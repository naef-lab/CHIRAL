#' Place each gene in a certain model based on rhythmicity with other conditions.
#'
#' nconds enumerates all possible models on rhythmicity in different conditions and
#' places gene expression into the most likely model based on Bayesian Inference Criterion.
#'
#' The function most used is \pkg{nconds}. 
#'
#' @docType package
#' @name nconds
NULL

###################################   MAIN    ###############################################
#' Place gene in a certain model based on rhythmicity with other conditions.
#'
#' @param dat Data object of gene expression.
#' Expect sample names to be \code{cond_00} where \code{cond} is 
#' condition and \code{00} is time of sample.
#' @param conds List of conditions matching to sample names in fpath.
#' @param t List of times at which each sample in fpath was taken
#' @param period Period of an oscillation. Default is 24.
#' @param out.prefix Prefix for output file names
#' @param write.intermediates \code{[TRUE]} If true, write design matrices and fit to disk
#' @param prepare.data \code{[FALSE]} Use TRUE if \code{dat} is in format with colnames contain only sample names and no gene names.
#' If first column of \code{dat} is gene names, then set to \code{FALSE}.
#' @return Return design matrices and rhythmic fits in each model and plots.
#' @examples
#' ncond("exprs.txt", 
#' conds = c(rep("Liver", 24), c(rep("Kidney",24)), c(rep("WFAT",24))), # 3 conditions
#' t <- rep(seq(0, 46, 2), 3)  # samples from times 0 to 46, sampled every 2 hours for 3 conditions.
#' 

nconds <- function(dat, conds, t, 
                   period = 24, 
                   out.prefix = 'nconds_output', 
                   write.intermediates = TRUE){
                   
  if (missing(conds) & missing(t)){
    nam.split <- strsplit(names(dat),"_")
    conds <- as.factor(sapply(nam.split, "[[", 1))
    t <- as.numeric(sapply(nam.split, "[[", 2))
    n.co <- nlevels(conds)
  }else{
    nam <- paste(conds,t,sep="_")
    names(dat) <- nam
    n.co = length(unique(conds))
  }
  
  my_mat = creat_matrix_list(t, conds, n.co, period)
  if (write.intermediates){
    save(my_mat, file = paste0(out.prefix, ".matrix_output.Robj"))
  }
  
  fit = parallel::mclapply(split(dat, rownames(dat)),
                 do_all,
                 t=t,
                 n.co=n.co,
                 period=period,
                 my_mat=my_mat,
                 conds=conds,
                 mc.cores=1)
  if (write.intermediates){
    save(fit, file = paste0(out.prefix, ".fit_output.Robj"))
  }
  
  data <- InsertFitToMat(fit, dat,n.co)
  plot_models(data,
              paste0(out.prefix, ".plots"), 
              t, 
              n.co, 
              conds,
              period = 24)
              return(data)
}

##########################################################################################

