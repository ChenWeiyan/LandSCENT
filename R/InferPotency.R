#' @title 
#' Infer Distinct potency states from cells' SR values
#' 
#' @aliases InferPotency
#'  
#' @description 
#' This function infers the discrete potency states of single cells and 
#' its distribution across the single cell population.
#' 
#' @param Integration.l
#' A list object from \code{CompSRana} function.
#' 
#' @param pheno.v
#' A phenotype vector for the single cells, of same length and order as the 
#' columns of \code{Integration.l$expMC}.
#' Function can also automatically extract phenotype information
#' from your original sce/cds data, please store the phenotype information
#' under the name of \code{phenoInfo}
#' 
#' @param diffvar
#' A logical. Default is TRUE.
#' Specifies whether the Gaussian mixture model to be fit assumes components 
#' to have different (default) or equal variance.
#' In the latter case, use *modelNames = c("E")*.
#' 
#' @param maxPS
#' Maximum number of potency states, when inferring discrete potency 
#' states of single cells. Default value is 5.
#' 
#' @return Integration.l
#' A list incorporates the input list with some new elements.
#' 
#' @return Integration.l$potencyState
#' Inferred discrete potency states for each single cell. It is indexed so 
#' that the index increases as the signaling entropy of the state decreases
#' 
#' @return Integration.l$distPSPH
#' If phenotype information provided, it will be a table giving the 
#' distribution of single-cells across potency states and 
#' phenotypes
#' 
#' @return Integration.l$prob
#' Table giving the probabilities of each potency state per phenotype value
#' 
#' @return Integration.l$hetPS
#' The normalised Shannon Index of potency per phenotype value
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
#' data(Example.m)
#' data(net13Jun12.m)
#' Integration.l <- DoIntegPPI(exp.m = Example.m[, c(1:58,61:84,86:98,100)], ppiA.m = net13Jun12.m)
#' data(SR.v)
#' Integration.l$SR <- SR.v
#' InferPotency.o <- InferPotency(Integration.l)
#' 
#' @import mclust
#' @import cluster
#' @import corpcor
#' @import MASS
#' @import parallel
#' @import Biobase
#' @import irlba
#' @import Rtsne
#' @importFrom dbscan dbscan
#' @importFrom isva EstDimRMT
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData<- 
#' @importFrom SummarizedExperiment colData
#' @importFrom stats as.dist
#' @importFrom stats cor
#' @importFrom stats median
#' @importFrom stats pnorm
#' @importFrom stats sd
#' @importFrom stats wilcox.test
#' @export
#'     
InferPotency <- function(Integration.l,
                         pheno.v = NULL,
                         diffvar = TRUE,
                         maxPS = 5)
{
    if (!is.null(Integration.l$data.sce)) {
        if (is.null(pheno.v)) {
            pheno.v <- 
                SingleCellExperiment::colData(Integration.l$data.sce)$phenoInfo
        }
        if (is.null(pheno.v)) {
            warning("No phenotype information, make sure it was stored as name of phenoInfo!")
        }
    }else if (!is.null(Integration.l$data.cds)) {
        if (is.null(pheno.v)) {
            pheno.v <- 
                Biobase::pData(Integration.l$data.cds)$phenoInfo
        }
        if (is.null(pheno.v)) {
            warning("No phenotype information, make sure it was stored as name of phenoInfo!")
        }
    }
    sr.v <- Integration.l$SR
    # Integration.l$phenoInfo <- pheno.v
    
    ### fit Gaussian Mixture Model for potency inference
    print("Fit Gaussian Mixture Model to Signaling Entropies")
    logitSR.v <- log2(sr.v / (1 - sr.v))
    if(diffvar == TRUE){ 
        ## default assumes different variance for clusters
        mcl.o <- Mclust(logitSR.v, G = seq_len(maxPS))
    }
    else {
        mcl.o <- Mclust(logitSR.v, G = seq_len(maxPS), modelNames = c("E"))
    }
    potS.v <- mcl.o$class
    nPS <- length(levels(as.factor(potS.v)))
    print(paste("Identified ",nPS," potency states",sep=""))
    for (i in seq_len(nPS)) {
        names(potS.v[which(potS.v == i)]) <- rep(paste("PS",i,sep=""), 
                                                 times = length(which(potS.v == i)))
    }
    
    mu.v <- mcl.o$param$mean
    sd.v <- sqrt(mcl.o$param$variance$sigmasq)
    avSRps.v <- (2^mu.v)/(1+2^mu.v)
    savSRps.s <- sort(avSRps.v, decreasing=TRUE, index.return=TRUE)
    spsSid.v <- savSRps.s$ix
    ordpotS.v <- match(potS.v,spsSid.v)
    
    if(!is.null(pheno.v)){
        nPH <- length(levels(as.factor(pheno.v)))
        distPSph.m <- table(pheno.v,ordpotS.v)
        print("Compute Shannon (Heterogeneity) Index for each Phenotype class")
        probPSph.m <- distPSph.m/apply(distPSph.m,1,sum)
        hetPS.v <- vector()
        for(ph in seq_len(nPH)){
            prob.v <- probPSph.m[ph,]
            sel.idx <- which(prob.v >0)
            hetPS.v[ph] <- 
                - sum(prob.v[sel.idx]*log(prob.v[sel.idx]))/log(nPS)
        }
        names(hetPS.v) <- rownames(probPSph.m)
        print("Done")
    }
    else {
        distPSph.m=NULL
        probPSph.m=NULL
        hetPS.v=NULL
    }
    
    if (!is.null(Integration.l$data.sce)) {
        colData(Integration.l$data.sce)$potencyState <- ordpotS.v
    }else if (!is.null(Integration.l$data.cds)) {
        pData(Integration.l$data.cds)$potencyState <- ordpotS.v
    }
    Integration.l$potencyState <- ordpotS.v
    
    if (!is.null(pheno.v)) {
        Integration.l$distPSPH <- distPSph.m
        Integration.l$prob <- probPSph.m
        Integration.l$hetPS <- hetPS.v
    }
    
    return(Integration.l)
}