#' @title 
#' Integration of gene expression matrix and PPI network
#' 
#' @aliases DoIntegPPI
#'  
#' @description 
#' This function finds the common genes between the scRNA-Seq data matrix 
#' and the genes present in the PPI network, and constructs the maximally 
#' connected subnetwork and reduced expression matrix for the computation 
#' of signaling entropy.
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
#' @return A list of two or three objects:
#' 
#' @return expMC
#' Reduced expression matrix with genes in the maximally connected subnetwork
#' 
#' @return adjMC
#' Adjacency matrix of the maximally connected subnetwork
#' 
#' @return data.sce/data.cds
#' Orginal input data objects
#' 
#' @details 
#' This function finds the common genes between the scRNA-Seq data matrix 
#' and the genes present in the PPI network, and constructs the maximally 
#' connected subnetwork and reduced expression matrix for the computation 
#' of signaling entropy.
#'  
#' @references 
#' Teschendorff AE, Tariq Enver. 
#' \emph{Single-cell entropy for accurate estimation of differentiation 
#' potency from a cell’s transcriptome.}
#' Nature communications 8 (2017): 15599.
#' doi:\href{https://doi.org/10.1038/ncomms15599}{
#' 10.1038/ncomms15599}.
#' 
#' Teschendorff AE, Banerji CR, Severini S, Kuehn R, Sollich P. 
#' \emph{Increased signaling entropy in cancer requires the scale-free 
#' property of protein interaction networks.}
#' Scientific reports 5 (2015): 9646.
#' doi:\href{https://doi.org/10.1038/srep09646}{
#' 10.1038/srep09646}.
#' 
#' Banerji, Christopher RS, et al. 
#' \emph{Intra-tumour signalling entropy determines clinical outcome 
#' in breast and lung cancer.}
#' PLoS computational biology 11.3 (2015): e1004115.
#' doi:\href{https://doi.org/10.1371/journal.pcbi.1004115}{
#' 10.1371/journal.pcbi.1004115}.
#' 
#' Teschendorff, Andrew E., Peter Sollich, and Reimer Kuehn.
#' \emph{Signalling entropy: A novel network-theoretical framework 
#' for systems analysis and interpretation of functional omic data.}
#' Methods 67.3 (2014): 282-293.
#' doi:\href{https://doi.org/10.1016/j.ymeth.2014.03.013}{
#' 10.1016/j.ymeth.2014.03.013}.
#' 
#' Banerji, Christopher RS, et al. 
#' \emph{Cellular network entropy as the energy potential in 
#' Waddington's differentiation landscape.}
#' Scientific reports 3 (2013): 3039.
#' doi:\href{https://doi.org/10.1038/srep03039}{
#' 10.1038/srep03039}.
#' 
#' @examples 
#' ### define a small network
#' ppiA.m <- matrix(0,nrow=10,ncol=10);
#' ppiA.m[1,] <- c(0,1,1,1,1);
#' for(r in 2:nrow(ppiA.m)){
#'   ppiA.m[r,1] <- 1;
#' }
#' rownames(ppiA.m) <- paste("G",1:10,sep="");
#' colnames(ppiA.m) <- paste("G",1:10,sep="");
#' 
#' ### define a positively valued expression matrix (20 genes x 10 samples)
#' exp.m <- matrix(rpois(20*10,8),nrow=20,ncol=10);
#' colnames(exp.m) <- paste("S",1:10,sep="");
#' rownames(exp.m) <- paste("G",1:20,sep="");
#' 
#' ### run integration function
#' Integrataion.l <- DoIntegPPI(exp.m,ppiA.m);
#' print(dim(Integrataion.l$expMC));
#' print(dim(Integrataion.l$adjMC));
#' 
#' @import Biobase
#' @import Matrix
#' @import scater
#' @importFrom scater normalise
#' @importFrom scater librarySizeFactors
#' @importFrom igraph graph.adjacency
#' @importFrom igraph clusters
#' @importFrom SummarizedExperiment assay
#' @importFrom DelayedArray isEmpty
#' @import monocle
#' @export
#'     
DoIntegPPI <- function(exp.m, 
                       ppiA.m)
{
    # set input data matrix class
    data.class <- class(exp.m)
    
    # select log-normalization methods based on data class
    if (data.class == "SingleCellExperiment") {
        sizeFactors(exp.m) <- scater::librarySizeFactors(exp.m)
        exp.m <- scater::normalise(exp.m, log_exprs_offset = 1.1)
        data.m <- Matrix::as.matrix(SummarizedExperiment::assay(exp.m, i = "logcounts"))
    }else if (data.class == "CellDataSet") {
        exp.m <- estimateSizeFactors(exp.m)
        data.m <- Matrix::as.matrix(t(t(Biobase::exprs(exp.m)) / 
                                          Biobase::pData(exp.m)[, 'Size_Factor']))
        data.m <- log2(data.m + 1.1)
    }else{
        data.m <- exp.m
    }
    
    if (min(data.m) == 0) {
        stop("Input matrix must have non-zero minimal value or set log_norm = TRUE!")
    }
    
    commonEID.v <- intersect(rownames(ppiA.m),rownames(data.m))
    
    if (DelayedArray::isEmpty(commonEID.v) == TRUE) {
        stop("scRNA-seq data should have the same gene identifier with the network!")
    }
    
    match(commonEID.v,rownames(data.m)) -> map1.idx
    expPIN.m <- data.m[map1.idx,]
    
    match(commonEID.v,rownames(ppiA.m)) -> map2.idx
    adj.m <- ppiA.m[map2.idx,map2.idx]
    
    gr.o <- igraph::graph.adjacency(adj.m,mode="undirected")
    comp.l <- igraph::clusters(gr.o)
    cd.v <- summary(factor(comp.l$member))
    mcID <- as.numeric(names(cd.v)[which.max(cd.v)])
    maxc.idx <- which(comp.l$member==mcID)
    adjMC.m <- adj.m[maxc.idx,maxc.idx]
    expMC.m <- expPIN.m[maxc.idx,]
    
    if (data.class == "SingleCellExperiment") {
        return(list(data.sce = exp.m, expMC = expMC.m, adjMC = adjMC.m))
    }else if (data.class == "CellDataSet") {
        return(list(data.cds = exp.m, expMC = expMC.m, adjMC = adjMC.m))
    }else{
        return(list(expMC = expMC.m, adjMC = adjMC.m))
    }
}
