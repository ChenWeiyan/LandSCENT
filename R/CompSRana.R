#' @title 
#' Computes the signaling entropy of given gene expression profiles 
#' over a given connected network
#' 
#' @aliases CompSRana
#'  
#' @description 
#' This is the main user function for computing signaling entropy of 
#' single cells. It takes as input the gene expression profile of a 
#' single cell and the adjacency matrix of a connected network. These 
#' inputs will be typically the output of the \code{DoIntePPI} function.
#' 
#' @param Integrataion.l
#' A list object ether from \code{DoIntePPI} function or \code{CompMaxSR} function.
#' 
#' @param local
#' A logical. If true (default), function computes the local signaling 
#' entropies of each gene in the network.
#' 
#' @param mc.cores
#' The number of cores to use, i.e. at most how many child processes will 
#' be run simultaneously. The option is initialized from environment variable 
#' MC_CORES if set. Must be at least one, and parallelization requires at 
#' least two cores.
#' 
#' @return A list incorporates the input list and four new elements:
#' 
#' @return SR
#' The global signaling entropy rate. If \code{maxSR} is provided,
#' then it is normalized by the maximum rate, hence a value between 0 and 1
#' 
#' @return inv
#' The stationary distribution of the sample
#' 
#' @return s
#' The unnormlised local entropies of each gene in every cell
#' 
#' @return ns
#' The normalised local entropies of each gene, so that each value is 
#' between 0 and 1
#' 
#' @details 
#' This function computes the signaling entropy for each individual 
#' sample (say a single cell) in the expression matrix provided as 
#' input. The input is an expression vector, not matrix, to allow 
#' user-flexibility for those wishing to run this in 
#' parallel (e.g. using parallel library).
#' 
#' @author Andrew E Teschendorff
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
#' ### integrate data and network
#' Integrataion.l <- DoIntegPPI(exp.m, ppiA.m);
#' 
#' ### compute SR values
#' Integrataion.l <- CompSRana(Integrataion.l);
#' 
#' ### output global signaling entropy
#' print(Integrataion.l$SR);
#' 
#' @import parallel
#' @import Biobase
#' @import SingleCellExperiment
#' @export
#'
CompSRana <- function(Integrataion.l,
                      local=TRUE,
                      mc.cores=1)
{
    ### see if maxSR provied or not
    if (is.null(Integrataion.l$maxSR) == TRUE) {
        maxSR <- NULL
    }else {
        maxSR <- Integrataion.l$maxSR
    }
    
    idx.l <- as.list(seq_len(ncol(Integrataion.l$expMC)))
    out.l <- mclapply(idx.l, CompSRanaPRL, 
                      exp.m=Integrataion.l$expMC, 
                      adj.m=Integrataion.l$adjMC,
                      local=local,
                      maxSR=maxSR,
                      mc.cores=mc.cores)
    SR.v <- sapply(out.l, function(v) return(v[[1]]))
    invP.v <- sapply(out.l, function(v) return(v[[2]]))
    S.v <- sapply(out.l, function(v) return(v[[3]]))
    NS.v <- sapply(out.l, function(v) return(v[[4]]))
    
    Integrataion.l$SR <- SR.v
    Integrataion.l$inv <- invP.v
    Integrataion.l$s <- S.v
    Integrataion.l$ns <- NS.v
    
    if (!is.null(Integrataion.l$data.sce)) {
        colData(Integrataion.l$data.sce)$SR <- SR.v
    }else if (!is.null(Integrataion.l$data.cds)) {
        pData(Integrataion.l$data.cds)$SR <- SR.v
    }
    
    return(Integrataion.l)
}

CompSRanaPRL <- function(idx,
                         exp.m,
                         adj.m,
                         local=TRUE,
                         maxSR=NULL)
{
    
    # compute outgoing flux around each node
    exp.v <- exp.m[,idx];
    sumexp.v <- as.vector(adj.m %*% matrix(exp.v,ncol=1));
    invP.v <- exp.v*sumexp.v;
    nf <- sum(invP.v);
    invP.v <- invP.v/nf;
    p.m <- t(t(adj.m)*exp.v)/sumexp.v;
    S.v <- apply(p.m,1,CompS);
    SR <- sum(invP.v*S.v);
    # if provided then normalise relative to maxSR
    if(is.null(maxSR)==FALSE){
        SR <- SR/maxSR;
    }
    if(local){
        NS.v <- apply(p.m,1,CompNS);
    }
    else {
        NS.v <- NULL;
    }
    return(list(sr=SR,inv=invP.v,s=S.v,ns=NS.v));
}
