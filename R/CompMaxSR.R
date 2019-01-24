#' @title 
#' Computes the maximum possible signaling entropy rate of a network
#' 
#' @aliases CompMaxSR
#'  
#' @description 
#' For a given maximally connected unweighted network, specified by an 
#' adjacency matrix, there is a maximum possible signaling entropy rate, 
#' which this function computes.
#' 
#' @param Integrataion.l
#' A list object from \code{DoIntegPPI} function.
#' 
#' @return Integrataion.l
#' A list incorporates the input list and an new element \code{maxSR} contains 
#' the maximum possible signaling entropy rate.
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
#' ### construct adjacency matrix
#' adj.m <- matrix(0,nrow=5,ncol=5);
#' adj.m[1,] <- c(0,1,1,1,1);
#' 
#' ### set diagonal elements to be 1
#' for(r in 2:nrow(adj.m)){
#'   adj.m[r,1] <- 1;
#' }
#' 
#' Integrataion.l <- list(adjMC = adj.m)
#' 
#' ### compute maximum possible signaling entropy rate
#' Integrataion.l <- CompMaxSR(Integrataion.l);
#' print(Integrataion.l$maxSR);
#' 
#' @importFrom igraph arpack
#' @export
#' 
CompMaxSR <- function(Integrataion.l){
    
    adj.m <- Integrataion.l$adjMC
    
    # find right eigenvector of adjacency matrix
    fa <- function(x,extra=NULL) {
        as.vector(adj.m %*% x)
    }
    ap.o <- igraph::arpack(fa,options=list(n=nrow(adj.m),nev=1,which="LM"), sym=TRUE)
    v <- ap.o$vectors
    lambda <- ap.o$values
    
    # maximum entropy
    MaxSR <- log(lambda)
    
    Integrataion.l$maxSR <- MaxSR
    
    return(Integrataion.l)
}
