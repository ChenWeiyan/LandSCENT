#' @title 
#' Computes the signaling entropy of given gene expression profiles 
#' over a given connected network
#' 
#' @aliases CompSRana
#'  
#' @description 
#' This is the main user function for computing signaling entropy of 
#' single cells. It takes as input the gene expression profile of 
#' single cells and the adjacency matrix of a connected network. These 
#' inputs will be typically the output of the \code{DoIntegPPI} function.
#' 
#' @param Integration.l
#' A list object from \code{DoIntegPPI} function.
#' 
#' @param local
#' A logical (default is FALSE). If TRUE, function computes the normalized 
#' local signaling entropies of each gene in the network.
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
#' The global signaling entropy rate. It is normalized by the 
#' maximum rate, hence a value between 0 and 1
#' 
#' @return inv
#' The stationary distribution of every sample
#' 
#' @return s
#' The unnormlised local entropies of each gene in every cell
#' 
#' @return ns
#' The normalised local entropies of each gene, so that each value is 
#' between 0 and 1
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
#' Integration.l <- DoIntegPPI(exp.m, ppiA.m);
#' 
#' ### compute SR values
#' Integration.l <- CompSRana(Integration.l);
#' 
#' ### output global signaling entropy
#' print(Integration.l$SR);
#' 
#' @import parallel
#' @import Biobase
#' @import SingleCellExperiment
#' @importFrom igraph arpack
#' @importFrom SummarizedExperiment colData<- 
#' @importFrom SummarizedExperiment colData
#' @export
#'
CompSRana <- function(Integration.l,
                      local = FALSE,
                      mc.cores=1)
{
    ### compute maxSR for SR normalization
    Integration.l <- CompMaxSR(Integration.l)
    maxSR <- Integration.l$maxSR
    
    idx.l <- as.list(seq_len(ncol(Integration.l$expMC)))
    out.l <- mclapply(idx.l, CompSRanaPRL, 
                      exp.m=Integration.l$expMC, 
                      adj.m=Integration.l$adjMC,
                      local=local,
                      maxSR=maxSR,
                      mc.cores=mc.cores)
    SR.v <- sapply(out.l, function(v) return(v[[1]]))
    invP.v <- sapply(out.l, function(v) return(v[[2]]))
    S.v <- sapply(out.l, function(v) return(v[[3]]))
    NS.v <- sapply(out.l, function(v) return(v[[4]]))
    
    Integration.l$SR <- SR.v
    Integration.l$inv <- invP.v
    Integration.l$s <- S.v
    Integration.l$ns <- NS.v
    
    if (!is.null(Integration.l$data.sce)) {
        colData(Integration.l$data.sce)$SR <- SR.v
    }else if (!is.null(Integration.l$data.cds)) {
        pData(Integration.l$data.cds)$SR <- SR.v
    }
    return(Integration.l)
}

CompMaxSR <- function(Integration.l){
    
    adj.m <- Integration.l$adjMC
    
    # find right eigenvector of adjacency matrix
    fa <- function(x,extra=NULL) {
        as.vector(adj.m %*% x)
    }
    ap.o <- igraph::arpack(fa,options=list(n=nrow(adj.m),nev=1,which="LM"), sym=TRUE)
    v <- ap.o$vectors
    lambda <- ap.o$values
    
    # maximum entropy
    MaxSR <- log(lambda)
    
    Integration.l$maxSR <- MaxSR
    
    return(Integration.l)
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

CompNS <- function(p.v){
    
    tmp.idx <- which(p.v>0);
    if(length(tmp.idx)>1){
        NLS <- -sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )/log(length(tmp.idx));
    }
    else {
        # one degree nodes have zero entropy, avoid singularity.
        NLS <- 0;
    }
    return(NLS);
}

CompS <- function(p.v){
    
    tmp.idx <- which(p.v>0);
    LS <-  - sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )
    return(LS);
}
