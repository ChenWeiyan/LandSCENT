#' @title 
#' Infer co-expression landmarks from cells' SR values
#' 
#' @aliases InferLandmark
#'  
#' @description 
#' This function identifies potency-coexpression clusters of single cells, 
#' called landmarks, and finally infers the dependencies of these landmarks 
#' which can aid in recontructing lineage trajectories in time course 
#' or development associated scRNA-Seq experiments.
#' 
#' @param Integration.l
#' A list object from \code{InferPotency} function.
#' 
#' @param pheno.v
#' A phenotype vector for the single cells, of same length and order as the 
#' columns of \code{Integration.l$expMC}.
#' Function can also automatically extract phenotype information
#' from your original sce/cds data, please store the phenotype information
#' under the name of \code{phenoInfo}.
#' 
#' @param pctG
#' A numeric. Percentage of all genes in \code{Integration.l$expMC} 
#' to select from each principal component in an SVD/PCA 
#' of \code{Integration.l$expMC}. The union set of all selected genes 
#' is then used for clustering. Default value is 0.01.
#' 
#' @param Component_use
#' A numeric. Specify the number of principal components in the PCA to use
#' in the downstream analysis. Default value is NULL.
#' 
#' @param reduceMethod
#' A character, either "PCA" or "tSNE". Indicates using PCA or
#' tSNE method to do dimension reduction.
#' 
#' @param clusterMethod 
#' A character, either "PAM" or "dbscan". Indicates using dbscan or
#' PAM method to do clustering.
#' 
#' @param k_pam
#' Only used when \code{clusterMethod} is set to be "PAM".
#' Maximum number of co-expression clusters, when performing clustering.
#' Default value is 9.
#' 
#' @param eps_dbscan
#' Only used when \code{clusterMethod} is set to be "dbscan". 
#' Size of the epsilon neighborhood. Default is 10.
#' 
#' @param minPts_dbscan
#' Only used when \code{clusterMethod} is set to be "dbscan".
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
#' @return Integration.l
#' A list incorporates the input list with a new list named 
#' \code{InferLandmark.l}.
#' 
#' @return InferLandmark.l
#' A list contains twelve objects:
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
#' of single cells mapping closest to the two landmarks
#' @return pcorLM
#' Partial correlation matrix of landmarks as estimated from the expression 
#' medoids
#' @return netLM
#' Adjacency matrix of landmarks specifying which partial correlations are 
#' significant
#' @return selectGene
#' Selected a group of genes for internal clustering
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
#' InferLandmark.o <- InferLandmark(InferPotency.o)
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
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData<- 
#' @importFrom SummarizedExperiment colData
#' @importFrom stats as.dist
#' @importFrom stats cor
#' @importFrom stats median
#' @importFrom stats pnorm
#' @importFrom stats sd
#' @importFrom stats var
#' @importFrom stats wilcox.test
#' @export
#'     
InferLandmark <- function(Integration.l,
                          pheno.v = NULL,
                          pctG = 0.01,
                          Component_use = NULL,
                          reduceMethod = c("PCA", "tSNE"),
                          clusterMethod = c("PAM", "dbscan"),
                          k_pam = 9,
                          eps_dbscan = 10,
                          minPts_dbscan = 5,
                          pctLM = 0.05,
                          pcorTH = 0.1)
{
    set.seed(2019)
    reduceMethod <- match.arg(reduceMethod)
    clusterMethod <- match.arg(clusterMethod)
    
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
    exp.m <- Integration.l$expMC
    sr.v <- Integration.l$SR
    ordpotS.v <- Integration.l$potencyState
    nPS <- length(levels(as.factor(ordpotS.v)))
    
    ### set an integer for gene selection
    ntop <- floor(pctG*nrow(exp.m))
    
    ### now cluster cells independently of SR
    ### select genes over which to cluster
    if (is.null(Integration.l$coordinates)) {
        
        tmp.m <- exp.m - rowMeans(exp.m)
        
        num_cell <- ncol(tmp.m)
        if (is.null(Component_use)) {
            if (num_cell >= 2e4) {
                print("Now estimating number of significant components of variation in scRNA-Seq data")
                print("Warning: It may take a long time since dimensionality of data matrix is big!")
                print("Consider specify number of compoents to use with 'Component_use'")
                rmtDim <- EstRMT(tmp.m)
            }else{
                print("Now estimating number of significant components of variation in scRNA-Seq data")
                rmtDim <- EstRMT(tmp.m)
            }
            print(paste("Number of significant components = ", rmtDim, sep=""))
        }else{
            Component_use <- min(Component_use, min(dim(tmp.m)))
            print(paste("Using top ", Component_use, " principal compoents for downstream analysis", sep = ""))
            rmtDim <- Component_use
        }
        
        if (rmtDim >= floor(0.5*min(dim(tmp.m)))) {
            svd.o <- svd(tmp.m)
        }else{
            svd.o <- irlba::irlba(tmp.m, nv = rmtDim)
        }
        
        ### select significant genes
        tmpG2.v <- vector()
        for(cp in seq_len(rmtDim)){
            tmp.s <- sort(abs(svd.o$u[, cp]), decreasing=TRUE, index.return=TRUE)
            tmpG2.v <- union(tmpG2.v, rownames(exp.m)[tmp.s$ix[seq_len(ntop)]])
        }
        selGcl.v <- tmpG2.v[which(tmpG2.v != "NA")]
        Integration.l$selectGene <- selGcl.v
        
        ### do dimension reduction
        if (reduceMethod == "tSNE") {
            print("Do dimension reduction via tSNE")
            pca.o <- irlba::prcomp_irlba(t(tmp.m), n = rmtDim, 
                                         center = TRUE, scale. = TRUE)
            irlba_pca_res <- pca.o$x
            tsne_res <- Rtsne::Rtsne(as.matrix(irlba_pca_res), dims = 2, 
                                     pca = FALSE)
            coordinates <- tsne_res$Y[, 1:2]
            Integration.l$coordinates <- coordinates
        }else {
            print("Do dimension reduction via PCA")
            coordinates <- svd.o$v[, 1:2]
            Integration.l$coordinates <- coordinates
        }
    }else {
        print("Dimension reduction done already")
        coordinates <- Integration.l$coordinates
        selGcl.v <- Integration.l$selectGene
    }
    
    map.idx <- match(selGcl.v,rownames(exp.m))
    
    ### now perform clustering of all cells over the selected genes
    if (clusterMethod == "dbscan") {
        print("Identifying co-expression clusters via dbscan")
        dbsc.o <- dbscan::dbscan(coordinates, eps = eps_dbscan, minPts = minPts_dbscan)
        
        Zero.idx <- length(which(dbsc.o$cluster == 0))
        if (Zero.idx != 0) {
            dbsc.o$cluster <- (dbsc.o$cluster + 1)
        }
        clust.idx <- dbsc.o$cluster
        
        names(clust.idx) <- colnames(Integration.l$expMC)
        k.opt <- length(unique(as.factor(dbsc.o$cluster)))
        print(paste("Inferred ",k.opt," clusters",sep=""))
        psclID.v <- paste("PS",ordpotS.v,"-CL",clust.idx,sep="")
    }else{
        print("Identifying co-expression clusters via PAM")
        distP.o <- as.dist( 0.5*(1-cor(exp.m[map.idx,])) )
        asw.v <- vector()
        for(k in 2:k_pam){
            pam.o <- pam(distP.o,k,stand=FALSE)
            asw.v[k-1] <- pam.o$silinfo$avg.width
        }
        k.opt <- which.max(asw.v)+1
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
    map.idx <- match(selGcl.v,rownames(exp.m))
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
    
    Integration.l$InferLandmark.l <- list(cl=clust.idx,pscl=psclID.v,
                                          distPSCL=distPSCL.m,medLM=medLM.m,
                                          srPSCL=srPSCL.v,srLM=srLM.v,
                                          distPHLM=distPHLM.m,cellLM=cellLM.v,
                                          cellLM2=cellLM2.v,adj=sadjLM.m,
                                          pcorLM=pcorLM.m,netLM=netLM.m)
    
    return(Integration.l)
}

EstRMT <- function(data.m)
{
    print("Centering and scaling matrix")
    M <- apply(data.m, 2, function(X) {
        (X - mean(X))/sqrt(var(X))
    })
    print("Done, now performing SVD")
    sigma2 <- stats::var(as.vector(M))
    Q <- nrow(data.m)/ncol(data.m)
    maxNEig <- min(dim(M))
    nEigEst <- floor(maxNEig/10)
    threshold.eigen <- sigma2 * (1 + 1/Q + 2 * sqrt(1/Q))
    
    c <- M/sqrt(nrow(M))
    
    if (min(dim(data.m)) >= 500) {
        print(paste("Using Fast IRLBA to approximate ",nEigEst," top singular values",sep=""))
        workDim <- 2 * nEigEst
        i.o <- irlba::irlba(c, nv = nEigEst, nu = nEigEst, maxit = 1000, work = workDim)
        evals <- (i.o$d^2)
        print("Done")
    }else{
        print("Performing full SVD since dimensionality of data matrix is not big")
        i.o <- svd(c)
        evals <- (i.o$d^2)
        print("Done")
    }
    
    intdim <- length(which(evals > threshold.eigen))
    
    return(intdim)
}
