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
#' Can be three major kinds of input:
#' One is a scRNA-Seq data matrix with rows labeling genes and columns 
#' labeling single cells. And it can be either a log-transformed data
#' matrix with minimal value around 0.1 (recommended), or an 
#' nonlog-transformed data matrix with minimal value 0.
#' The other two kinds of input can be either a "SingleCellExperiment"
#' class object or a "CellDataSet" class object
#' 
#' @param ppiA.m
#' The adjacency matrix of a user-given PPI network with rownames and 
#' colnames labeling genes (same gene identifier as in \code{exp.m})
#' 
#' @param log_trans
#' A logical. Whether to do log-transformation on the input data
#' matrix or not. Default is FALSE
#' 
#' @param mc.cores
#' The number of cores to use, i.e. at most how many child processes will 
#' be run simultaneously. The option is initialized from environment variable 
#' MC_CORES if set. Must be at least one (default), and parallelization 
#' requires at least two cores.
#' 
#' @param pheno.v
#' A phenotype vector for the single cells, of same length and order as the 
#' columns of \code{exp.m}.
#' Function can also automatically extract phenotype information
#' from your original sce/cds data, please store the phenotype information
#' as name of \code{phenoInfo}.
#' 
#' @param coordinates
#' Optional. The previous reduced dimension coordinates, with rows lalabeling cells 
#' and two colums labeling reduced dimensions
#' 
#' @param PLOT
#' A logical. Decides whether to generate (default) the landSR figure 
#' or not.
#' 
#' @param PDF
#' A logical. Output figure via pdf file or not, default is TRUE
#' 
#' @return Integration.l
#' A list contains input information and SR values, potency states and more
#' other results.
#' 
#' @return A PDF file
#' If PDF is TRUE(default), then it will automatically generate a pdf file
#' ploting cell density against potency states.
#' 
#' @export
#' 
DoLandSCENT <- function(exp.m, 
                        ppiA.m,
                        log_trans = FALSE,
                        mc.cores = 1,
                        pheno.v = NULL,
                        coordinates = NULL,
                        PLOT = TRUE,
                        PDF = TRUE)
{
    ### integrate scRNA-seq and PPI network
    Integration.l <- DoIntegPPI(exp.m = exp.m,
                                ppiA.m = ppiA.m,
                                log_trans = log_trans)
    
    ### compute SR values for every cells
    Integration.l <- CompSRana(Integration.l, mc.cores = mc.cores)
    
    ### infer potency states
    Integration.l <- InferPotency(Integration.l,
                                  pheno.v = pheno.v,
                                  diffvar = TRUE,
                                  maxPS = 5)
    
    ### infer landmarks
    Integration.l <- InferLandmark (Integration.l,
                                    pheno.v = pheno.v,
                                    pctG = 0.01,
                                    reduceMethod = "PCA",
                                    clusterMethod = "PAM",
                                    k_pam = 15,
                                    eps_dbscan = 10,
                                    minPts_dbscan = 5,
                                    pctLM = 0.05,
                                    pcorTH = 0.1)
    
    ### generate figures
    if (PLOT == TRUE) {
        Integration.l <- Plot_LandSR(Integration.l, coordinates = coordinates, PDF = PDF)
        Integration.l <- Plot_CellSR(Integration.l, coordinates = coordinates, PDF = PDF)
    }
    
    return(Integration.l)
}