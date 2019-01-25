#' @title 
#' Runs the LandSCENT algorithm
#' 
#' @aliases DoLandSCENT
#'  
#' @description 
#' Main user function implement LandSCENT. This function is the typical 
#' workflow of the whole package for you to easily use.
#' 
#' @param exp.m
#' Can be three kinds of input.
#' One is a scRNA-Seq data matrix with rows labeling genes and columns 
#' labeling single cells, which is not log-transformed (can be FPKM, TPM etc.).
#' And rownames must be provided and must be the same with the 
#' rownames of \code{ppiA.m}.
#' The other kinds of two input can be either a "SingleCellExperiment"
#' class object or a "CellDataSet" class object
#' 
#' @param ppiA.m
#' The adjacency matrix of a user-given PPI network with rownames and 
#' colnames labeling genes (same gene identifier as in \code{exp.m}).
#' 
#' @param mc.cores
#' The number of cores to use, i.e. at most how many child processes will 
#' be run simultaneously. The option is initialized from environment variable 
#' MC_CORES if set. Must be at least one, and parallelization requires at 
#' least two cores.
#' 
#' @param pheno.v
#' A phenotype vector for the single cells, of same length and order as the 
#' columns of \code{exp.m}.
#' Function can also automatically extract phenotype information
#' from your original sce/cds data, please store the phenotype information
#' as name of \code{phenoInfo}.
#' 
#' @param reducedMatrix
#' The previous reduced dimension matrix, with rows lalabeling cells and two
#' colums labeling reduced dimensions.(Optional, only used when 
#' \code{reduceDim} set to \code{FALSE})
#' 
#' @param reduceDim
#' A logical, do deminsion reduction or not before generate the plots.
#' Default is TRUE.
#' 
#' @param PLOT
#' A logical. Decides whether to generate (default) the landSR figure 
#' and CellSR figure or not.
#' 
#' @param PDF
#' The figure file output format, via pdf (TRUE) file or not, default is TRUE.
#' 
#' @return Integrataion.l
#' A list contains input information and SR values, potency states and more
#' other results. Typically, it contains the same values as the output of
#' \code{PotencyInfer}.
#' 
#' @return PDF file
#' If PDF is TRUE(default), then it will automatically generate a pdf file
#' ploting cell density against potency states.
#' 
#' @export
#' 
DoLandSCENT <- function(exp.m, 
                        ppiA.m,
                        mc.cores = 1,
                        pheno.v = NULL,
                        reducedMatrix = NULL,
                        reduceDim = TRUE,
                        PLOT = TRUE,
                        PDF = TRUE)
{
    ### integrate scRNA-seq and PPI network
    Integrataion.l <- DoIntegPPI(exp.m = exp.m,
                                 ppiA.m = ppiA.m)
    
    ### compute maximum SR value
    Integrataion.l <- CompMaxSR(Integrataion.l)
    
    ### compute SR values for every cells
    Integrataion.l <- CompSRana(Integrataion.l, mc.cores = mc.cores)
    
    ### infer potency state and call landmarks
    Integrataion.l <- PotencyInfer(Integrataion.l, 
                                   pheno.v = pheno.v, 
                                   mixmod=TRUE,
                                   maxPS=5,
                                   pctG=0.01,
                                   kmax=9,
                                   reduceMethod = "PAM",
                                   num_clusters = NULL,
                                   epsMax = 10,
                                   minPts = 5,
                                   pctLM=0.05,
                                   pcorTH=0.1)
    
    ### generate figures
    if (PLOT == TRUE) {
        Plot_LandSR(Integrataion.l, reducedMatrix = reducedMatrix, reduceDim = reduceDim, PDF = PDF)
    }
    
    return(Integrataion.l)
}