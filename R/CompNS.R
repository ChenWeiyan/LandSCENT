#' @title 
#' Computes the local signaling entropy normalised to lie between 0 and 1
#' 
#' @aliases CompNS
#'  
#' @description 
#' Computes the normalized local signaling entropy for a gene. This is an 
#' internal function which the user does not need to invoke.
#' 
#' @param p.v
#' A vector specifying the outgoing probabilities from a given node (gene).
#' Values should add to 1.
#' 
#' @return NLS
#' A value indicates the normalised local entropy.
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
#' ### specify a probability vector
#' p.v <- c(0.25,0.25,0.5);
#' NLS <- CompNS(p.v);
#' print(NLS);
#' 
#' @export
#' 
CompNS <- function(p.v){
    
    tmp.idx <- which(p.v>0);
    if(length(tmp.idx)>1){
        NLS <-  -sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )/log(length(tmp.idx));
    }
    else {
        # one degree nodes have zero entropy, avoid singularity.
        NLS <- 0;
    }
    return(NLS);
}
