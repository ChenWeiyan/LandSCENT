#' @title 
#' Construct diffusion map and select root cell
#' 
#' @aliases DoDiffusionMap
#'  
#' @description 
#' This function constructs diffusion map with highly variable
#' genes and visualize the result in 2D or 3D plots
#' 
#' @param Integration.l
#' Typically, it is the output from \code{InferLandmark} function
#' 
#' @param mean_gap
#' The threshold for HVG gene selection based on the mean expression
#' level. Default value is 1
#' 
#' @param sd_gap
#' The threshold for HVG gene selection based on the expression
#' level's standard deviation. Default value is 1
#' 
#' @param root
#' Define the root of the trajectory. \code{cell} indicates the
#' root is the cell with the highest SR value among all cells;
#' \code{state} indicates the root is the cell with the highest 
#' SR value within the highest inferred potency state
#' 
#' @param num_comp
#' A numeric. The dimension of diffusion map, which is related 
#' to the following plot generation. Default is 3
#' 
#' @param k
#' Number of nearest neighbors to consider. Default is 30
#' 
#' @return An updated list objects, add three elements:
#' 
#' @return DM
#' A diffusion map object
#' 
#' @return DMEigen
#' The diffusion components of the diffusion map
#' 
#' @return root
#' The index of root cell for generating diffusion map plots
#' 
#' @references 
#' Teschendorff AE, Tariq Enver. 
#' \emph{Single-cell entropy for accurate estimation of differentiation 
#' potency from a cell’s transcriptome.}
#' Nature communications 8 (2017): 15599.
#' doi:\href{https://doi.org/10.1038/ncomms15599}{
#' 10.1038/ncomms15599}.
#' 
#' Teschendorff, Andrew E., Peter Sollich, and Reimer Kuehn.
#' \emph{Signalling entropy: A novel network-theoretical framework 
#' for systems analysis and interpretation of functional omic data.}
#' Methods 67.3 (2014): 282-293.
#' doi:\href{https://doi.org/10.1016/j.ymeth.2014.03.013}{
#' 10.1016/j.ymeth.2014.03.013}.
#' 
#' @examples 
#' \dontrun{
#' data(Example.m)
#' data(net13Jun12.m)
#' Integration.l <- DoIntegPPI(exp.m = Example.m[, c(1:58,61:84,86:98,100)], ppiA.m = net13Jun12.m)
#' data(SR.v)
#' Integration.l$SR <- SR.v
#' InferPotency.o <- InferPotency(Integration.l)
#' DoDM.o <- DoDiffusionMap(InferPotency.o)}
#' 
#' @importFrom destiny DiffusionMap
#' @importFrom destiny eigenvectors
#' @importFrom cluster pam
#' @export
#'     
DoDiffusionMap <- function(Integration.l,
                           mean_gap = 1,
                           sd_gap = 1,
                           root = c("cell", "state"),
                           num_comp = 3,
                           k = 30){
  
  data.m <- Integration.l$data
  
  if (num_comp < 2) {
    stop("Dimension for diffusion maps should be at least 2!")
  }
  
  if (length(root) > 1) {
    print("Specified more than one type of root, using the first type instead!")
    root <- root[1]
  }
  
  print("Selecting highly variable genes for diffusion map construction")
  print(paste("Remove genes below ", sd_gap, " standard deviation", sep = ""))
  sd.v <- apply(data.m, 1, sd)
  print(paste("Remove genes below ", mean_gap, " mean", sep = ""))
  mean.v <- apply(data.m, 1, mean)
  selG.idx <- intersect(which(mean.v > mean_gap), 
                        which(sd.v > sd_gap))
  
  exp.m <- data.m[selG.idx ,]
  print(paste("Selected ",length(selG.idx)," genes",sep=""))
  
  print("Constructing diffusion maps")
  if (is.null(Integration.l$DM)) {
    dm <- destiny::DiffusionMap(t(exp.m), k = k)
    dms <- eigenvectors(dm)[, seq_len(num_comp)]
    dms <- data.frame(dms)
    print("Done")
  }else{
    print("Diffusion map was previously constructed already, skipped!")
    dm <- Integration.l$DM
    dms <- eigenvectors(dm)[, seq_len(num_comp)]
    dms <- data.frame(dms)
  }
  
  print("Selecting root cell")
  
  ### give warning based on cell number
  if (ncol(exp.m) <= 200) {
    warning("There are too little cells(<= 200), be careful with Plot_DiffusionMap function!")
    warning("Plot diffusion map colored by SR first befrore color it with diffusion pseudo time(DPT)!")
  }
  
  if (root == "cell") {
    print("Taking cell with the highest entropy rate as the root cell")
    root.idx <- which(Integration.l$SR == max(Integration.l$SR))
  }else if (root == "state"){
    print("Taking cell with the highest entropy rate within the highest potency cluster as the root cell")
    high.idx <- which(Integration.l$potencyState == 1)
    
    if (!is.null(Integration.l$coordinates)) {
      tsne.o <- t(Integration.l$coordinates)
      distP.o <- as.dist( 0.5*(1-cor(tsne.o[, high.idx])) )
      pam.o <- cluster::pam(distP.o, k = 3, stand = FALSE)
      clust.idx <- pam.o$cluster
    }else{
      distP.o <- as.dist( 0.5*(1-cor(exp.m[, high.idx])) )
      pam.o <- cluster::pam(distP.o, k = 3, stand = FALSE)
      clust.idx <- pam.o$cluster
    }
    ### seprate clusters
    # tmp <- as.data.frame(table(clust.idx))
    # cluster.id <- as.numeric(order(tmp$Freq)[3])
    cl1.idx <- which(clust.idx == 1)
    cl2.idx <- which(clust.idx == 2)
    cl3.idx <- which(clust.idx == 3)
    
    ### claculate the median
    median1 <- median(Integration.l$SR[cl1.idx])
    median2 <- median(Integration.l$SR[cl2.idx])
    median3 <- median(Integration.l$SR[cl3.idx])
    
    ### exclude isolated cells
    if (length(cl1.idx) < (0.1 * length(high.idx))) {
      median1 <- 0
    }else if (length(cl2.idx) < (0.1 * length(high.idx))) {
      median2 <- 0
    }else if (length(cl3.idx) < (0.1 * length(high.idx))) {
      median3 <- 0
    }
    cluster.id <- order(c(median1,
                          median2,
                          median3))[3]
    group.idx <- high.idx[which(pam.o$cluster == cluster.id)]
    maxSR <- max(Integration.l$SR[group.idx])
    root.idx <- high.idx[which(Integration.l$SR[high.idx] == maxSR)]
  }
  
  Integration.l$DM <- dm
  Integration.l$DMEigen <- dms
  Integration.l$root <- root.idx
  
  return(Integration.l)
  
}










