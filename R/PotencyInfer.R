#' @title 
#' Infer Distinct potency states from cells' SR values
#' 
#' @aliases PotencyInfer
#'  
#' @description 
#' This function infers the discrete potency states of single cells, 
#' its distribution across the single cell population, identifies 
#' potency-coexpression clusters of single cells, called landmarks, 
#' and finally infers the dependencies of these landmarks which can 
#' aid in recontructing lineage trajectories in time course scRNA-Seq 
#' experiments.
#' 
#' @details 
#' This function infers the discrete potency states of single cells, 
#' its distribution across the single cell population, identifies 
#' potency-coexpression clusters of single cells, called landmarks, 
#' and finally infers the dependencies of these landmarks which can 
#' aid in recontructing lineage trajectories in time course scRNA-Seq 
#' experiments.
#' 
#' @param Integrataion.l
#' A list object from \code{CompSRana} function.
#' 
#' @param pheno.v
#' A phenotype vector for the single cells, of same length and order as the 
#' columns of \code{Integrataion.l$expMC}.
#' Function can also automatically extract phenotype information
#' from your original sce/cds data, please store the phenotype information
#' as name of \code{phenoInfo}
#' 
#' @param mixmod
#' A logical. Default is TRUE.
#' Specifies whether the Gaussian mixture model to be fit assumes components 
#' to have different (default) or equal variance.
#' In the latter case, use *mixmod=c("E")*.
#' 
#' @param maxPS
#' Maximum number of potency states to allow, when inferring discrete potency 
#' states of single cells. Default value is 5.
#' 
#' @param pctG
#' Percentage of all genes in \code{Integrataion.l$expMC} to select from
#' each principal component in an SVD/PCA of \code{Integrataion.l$expMC}.
#' The union set of all selected genes is then used for clustering. 
#' Default value is 0.01.
#' 
#' @param kmax
#' Maximum number of co-expression clusters allowed when performing clustering.
#' Default value is 9. Larger values are not allowed.
#' 
#' @param reduceMethod 
#' A character, either "tSNE" or "PAM"(default). Indicates using tSNE or
#' PAM method to do dimension reduction.
#' 
#' @param num_clusters
#' Only used when data reduced by tSNE, specify the number of clusters.
#' 
#' @param epsMax
#' Maximum value of scanning when implementing \code{dbscan} to determine
#' the cluster number.
#' 
#' @param minPts
#' Number of minimum points in the eps region (for core points).
#' Default is 5 points.
#' 
#' @param pctLM
#' Percentage of total number of single cells to allow as a minimum size for 
#' selecting interesting landmarks i.e. potency-coexpression clusters of single 
#' cells. Default value is 0.05.
#' 
#' @param pcorTH
#' Threshold for calling significant partial correlations. Default value is 0.1.
#' Usually, single-cell experiments profile large number of cells, so 0.1 is a 
#' sensible threshold.
#' 
#' @return potencyInfer.l
#' A list incorporates the input list and a new list \code{potencyInfer.l}
#' contains sixteen objects:
#' @return potS
#' Inferred discrete potency states for each single cell. It is indexed so 
#' that the index increases as the signaling entropy of the state decreases
#' @return distPSPH
#' Table giving the distribution of single-cells across potency states and 
#' phenotypes
#' @return prob
#' Table giving the probabilities of each potency state per phenotype value
#' @return hetPS
#' The normalised Shannon Index of potency per phenotype value
#' @return cl
#' The co-expression clustering index for each single cell
#' @return pscl
#' The potency coexpression clustering label for each single cell
#' @return distPSCL
#' The distribution of single cell numbers per potency state and coexpression 
#' cluster
#' @return medLM
#' A matrix of medoids of gene expression for the selected landmarks
#' @return srPSCL
#' The average signaling entropy of single cells in each potency coexpression 
#' cluster
#' @return srLM
#' The average signaling entropy of single cells in each landmark
#' @return distPHLM
#' Table giving the distribution of single cell numbers per phenotype and 
#' landmark
#' @return cellLM
#' Nearest landmark for each single cell
#' @return cellLM2
#' A vector specifying the nearest and next-nearest landmark for each single 
#' cell
#' @return adj
#' Weighted adjacency matrix between landmarks with entries giving the number 
#' of single 
#' cells mapping closest to the two landmarks
#' @return pcorLM
#' Partial correlation matrix of landmarks as estimated from the expression 
#' medoids
#' @return netLM
#' Adjacency matrix of landmarks specifying which partial correlations are 
#' significant
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
#' Integrataion.l <- DoIntegPPI(exp.m = Example.m[, c(1:58,61:84,86:98,100)], ppiA.m = net13Jun12.m)
#' Integrataion.l <- CompMaxSR(Integrataion.l)
#' data(SR.v)
#' Integrataion.l$SR <- SR.v
#' potencyInfer.o <- PotencyInfer(Integrataion.l)
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
PotencyInfer <- function(Integrataion.l,
                    pheno.v=NULL,
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
{
    if (!is.null(Integrataion.l$data.sce)) {
        if (is.null(pheno.v)) {
            pheno.v <- 
                SingleCellExperiment::colData(Integrataion.l$data.sce)$phenoInfo
        }
        if (is.null(pheno.v)) {
            warning("No phenotype information, make sure it was stored as name of phenoInfo!")
        }
    }else if (!is.null(Integrataion.l$data.cds)) {
        if (is.null(pheno.v)) {
            pheno.v <- 
                Biobase::pData(Integrataion.l$data.cds)$phenoInfo
        }
        if (is.null(pheno.v)) {
            warning("No phenotype information, make sure it was stored as name of phenoInfo!")
        }
    }
    exp.m <- Integrataion.l$expMC
    sr.v <- Integrataion.l$SR
    
    ### set an integer for gene selection
    ntop <- floor(pctG*nrow(exp.m))
    
    ### fit Gaussian Mixture Model for potency inference
    print("Fit Gaussian Mixture Model to Signaling Entropies")
    logitSR.v <- log2(sr.v / (1 - sr.v))
    if(mixmod == TRUE){ 
        ## default assumes different variance for clusters
        mcl.o <- Mclust(logitSR.v, G = seq_len(maxPS))
    }
    else {
        mcl.o <- Mclust(logitSR.v, G = seq_len(maxPS), modelNames = c("E"))
    }
    potS.v <- mcl.o$class
    nPS <- length(levels(as.factor(potS.v)))
    print(paste("Identified ",nPS," potency states",sep=""))
    # names(potS.v) <- paste("PS",seq_len(nPS),sep="")
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
        nPH <- length(levels(as.factor(pheno.v)));
        distPSph.m <- table(pheno.v,ordpotS.v)
        print("Compute Shannon (Heterogeneity) Index for each Phenotype class");
        probPSph.m <- distPSph.m/apply(distPSph.m,1,sum);
        hetPS.v <- vector();
        for(ph in seq_len(nPH)){
            prob.v <- probPSph.m[ph,];
            sel.idx <- which(prob.v >0);
            hetPS.v[ph] <- 
                - sum(prob.v[sel.idx]*log(prob.v[sel.idx]))/log(nPS);
        }
        names(hetPS.v) <- rownames(probPSph.m);
        print("Done");
    }
    else {
        distPSph.m=NULL; probPSph.m=NULL; hetPS.v=NULL;
    }
    
    ### now cluster cells independently of SR
    ### select genes over which to cluster
    print("Using RMT to estimate number of significant components of variation in scRNA-Seq data");
    tmp.m <- exp.m - rowMeans(exp.m);
    rmt.o <- isva::EstDimRMT(tmp.m, plot = FALSE);
    svd.o <- svd(tmp.m);
    tmpG2.v <- vector();
    print(paste("Number of significant components=",rmt.o$dim,sep=""));
    for(cp in seq_len(rmt.o$dim)){
        tmp.s <- sort(abs(svd.o$u[,cp]),decreasing=TRUE,index.return=TRUE);
        tmpG2.v <- union(tmpG2.v,rownames(exp.m)[tmp.s$ix[seq_len(ntop)]]);
    }
    selGcl.v <- tmpG2.v;
    
    ### now perform clustering of all cells over the selected genes
    if (reduceMethod == "tSNE") {
        print("Identifying co-expression clusters via tSNE")
        if (is.null(num_clusters)) {
            stop("Please specify cluster numbers if use tSNE for dimension reduction!")
        }
        map.idx <- match(selGcl.v,rownames(exp.m))
        data_tsne.m <- exp.m[map.idx ,]
        irlba_res <- irlba::prcomp_irlba(t(data_tsne.m), n = rmt.o$dim
                                         , center = TRUE)
        irlba_pca_res <- irlba_res$x
        topDim_pca <- irlba_pca_res
        tsne_res <- Rtsne::Rtsne(as.matrix(topDim_pca), dims = 2, 
                                 pca = FALSE)
        reducedMatrix <- tsne_res$Y[, 1:2]
        dbsc.v <- vector()
        eps.v <- vector()
        id.x <- 1
        for(k in seq(from = 1, to = epsMax, by = 0.2)){
            dbsc.o <- dbscan::dbscan(reducedMatrix, eps = k, minPts = minPts)
            clust.idx <- which(dbsc.o$cluster>0)
            label.v <- dbsc.o$cluster[clust.idx]
            if (length(unique(as.factor(label.v))) == num_clusters) {
                break()
            }
            dbsc.v[id.x] <- length(unique(as.factor(dbsc.o$cluster)))
            eps.v[id.x] <- k
            id.x <- id.x + 1
        }
        id.x <- which((dbsc.v - num_clusters) == min(dbsc.v - num_clusters))[1]
        eps.value <- min(k, eps.v[id.x])
        dbsc.o <- dbscan::dbscan(reducedMatrix, eps = eps.value, minPts = minPts)
        clust.idx <- dbsc.o$cluster
        k.opt <- length(unique(as.factor(dbsc.o$cluster)))
        print(paste("Inferred ",k.opt," clusters",sep=""))
        psclID.v <- paste("PS",ordpotS.v,"-CL",clust.idx,sep="")
    }else{
        print("Identifying co-expression clusters via PAM");
        map.idx <- match(selGcl.v,rownames(exp.m));
        distP.o <- as.dist( 0.5*(1-cor(exp.m[map.idx,])) );
        asw.v <- vector();
        for(k in 2:kmax){
            pam.o <- pam(distP.o,k,stand=FALSE);
            asw.v[k-1] <- pam.o$silinfo$avg.width
        }
        k.opt <- which.max(asw.v)+1;
        pam.o <- pam(distP.o,k=k.opt,stand=FALSE)
        clust.idx <- pam.o$cluster
        print(paste("Inferred ",k.opt," clusters",sep=""))
        psclID.v <- paste("PS",ordpotS.v,"-CL",clust.idx,sep="")
    }
    ### identify landmark clusters
    print("Now identifying landmarks (potency co-expression clusters)")
    distPSCL.m <- table(paste("CL",clust.idx,sep=""),
                        paste("PS",ordpotS.v,sep=""))
    sizePSCL.v <- as.vector(distPSCL.m)
    namePSCL.v <- vector()
    ci <- 1
    for(ps in seq_len(nPS)){
        for(cl in seq_len(k.opt)){
            namePSCL.v[ci] <- paste("PS",ps,"-CL",cl,sep="")
            ci <- ci+1
        }
    }
    names(sizePSCL.v) <- namePSCL.v
    ldmkCL.idx <- which(sizePSCL.v > pctLM*ncol(exp.m))
    print(paste("Identified ",length(ldmkCL.idx)," Landmarks",sep=""))
    
    ### distribution of phenotypes among LMs
    if(!is.null(pheno.v)){
        tab.m <- table(pheno.v,psclID.v)
        tmp.idx <- match(names(sizePSCL.v)[ldmkCL.idx],colnames(tab.m))
        distPHLM.m <- tab.m[,tmp.idx];
    }
    else {
        distPHLM.m <- NULL;
    }
    
    ### medoids
    print("Constructing expression medoids of landmarks");
    med.m <- matrix(0,nrow=length(selGcl.v),ncol=nPS*k.opt);
    srPSCL.v <- vector();
    ci <- 1;
    for(ps in seq_len(nPS)){
        for(cl in seq_len(k.opt)){
            tmpS.idx <- 
                intersect(which(ordpotS.v==ps),which(clust.idx==cl));
            med.m[,ci] <- apply(matrix(exp.m[map.idx,tmpS.idx],
                                       nrow=length(map.idx)),1,median);
            srPSCL.v[ci] <- mean(sr.v[tmpS.idx]);
            ci <- ci+1;
        }
    }
    names(srPSCL.v) <- namePSCL.v;
    srLM.v <- srPSCL.v[ldmkCL.idx];
    medLM.m <- med.m[,ldmkCL.idx];
    colnames(medLM.m) <- namePSCL.v[ldmkCL.idx];
    rownames(medLM.m) <- selGcl.v;
    ### now project each cell onto two nearest landmarks
    print("Inferring dependencies/trajectories/transitions between landmarks");
    cellLM2.v <- vector(); cellLM.v <- vector();
    for(c in seq_len(ncol(exp.m))){
        distCellLM.v <- 0.5*(1-as.vector(cor(exp.m[map.idx,c],medLM.m)));
        tmp.s <- sort(distCellLM.v,decreasing=FALSE,index.return=TRUE);
        cellLM2.v[c] <- paste("LM",tmp.s$ix[1],"-LM",tmp.s$ix[2],sep="");
        cellLM.v[c] <- colnames(medLM.m)[tmp.s$ix[1]];
    }
    
    adjLM.m <- matrix(0,nrow=ncol(medLM.m),ncol=ncol(medLM.m));
    rownames(adjLM.m) <- colnames(medLM.m);
    colnames(adjLM.m) <- colnames(medLM.m);
    for(lm1 in seq_len(ncol(medLM.m))){
        for(lm2 in seq_len(ncol(medLM.m))){
            adjLM.m[lm1,lm2] <- 
                length(which(cellLM2.v==paste("LM",lm1,"-LM",lm2,sep="")));
        }
    }
    sadjLM.m <- adjLM.m + t(adjLM.m);
    corLM.m <- cor(medLM.m);
    pcorLM.m <- cor2pcor(corLM.m);
    rownames(pcorLM.m) <- rownames(corLM.m);
    colnames(pcorLM.m) <- rownames(corLM.m);    
    netLM.m <- pcorLM.m;
    diag(netLM.m) <- 0;
    netLM.m[pcorLM.m < pcorTH] <- 0;
    netLM.m[pcorLM.m > pcorTH] <- 1;    
    
    if (!is.null(Integrataion.l$data.sce)) {
        colData(Integrataion.l$data.sce)$potencyState <- ordpotS.v
    }else if (!is.null(Integrataion.l$data.cds)) {
        pData(Integrataion.l$data.cds)$potencyState <- ordpotS.v
    }
    Integrataion.l$potencyState <- ordpotS.v
    
    if (!is.null(pheno.v)) {
        Integrataion.l$distPSPH <- distPSph.m
    }
    
    Integrataion.l$potencyInfer.l <- list(potS=ordpotS.v,distPSPH=distPSph.m,prob=probPSph.m,
                                          hetPS=hetPS.v,cl=clust.idx,pscl=psclID.v,
                                          distPSCL=distPSCL.m,medLM=medLM.m,srPSCL=srPSCL.v,
                                          srLM=srLM.v,distPHLM=distPHLM.m,cellLM=cellLM.v,
                                          cellLM2=cellLM2.v,adj=sadjLM.m,pcorLM=pcorLM.m,
                                          netLM=netLM.m)
    
    
    return(Integrataion.l)
}
