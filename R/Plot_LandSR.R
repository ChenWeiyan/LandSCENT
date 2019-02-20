#' @title 
#' Plot Landscape with potency states
#' 
#' @aliases Plot_LandSR
#'  
#' @description 
#' Density based visualization tool to generate figures that give the cell 
#' density distribution against potency states
#' 
#' @param Integration.l
#' Typically, it is the output from \code{PotencyInfer} function
#' 
#' @param reducedMatrix
#' The previous reduced dimension matrix, with rows lalabeling cells and two
#' colums labeling reduced dimensions. Only used when \code{reduceDim} is 
#' set to be FALSE
#' 
#' @param num_dim
#' number of dimension used for PCA
#' 
#' @param max_components
#' The dimensionality of the reduced space
#' 
#' @param num_grid
#' Number of grid points in each direction. Can be scalar or a length-2 
#' integer vector
#' 
#' @param reduceDim
#' A logical, do deminsion reduction or not, default is TRUE
#' 
#' @param mean_cutoff
#' Threshold of mean expression level for gene used to perform 
#' dimension reduction. Default is 1
#' 
#' @param sd_cutoff
#' Threshold of standard deviation for gene used to perform 
#' dimension reduction. Default is 0.5
#' 
#' @param phi
#' The angles defining the viewing direction. phi gives the colatitude
#' 
#' @param theta
#' The angles defining the viewing direction. theta gives the 
#' azimuthal direction
#' 
#' @param colpersp
#' Color palette to be used for the colvar variable, i.e. specific 
#' potency state cell density.
#' If colpersp is NULL (default) and colvar is specified, then a 
#' red-yellow-blue colorscheme will be used.
#' 
#' @param colimage
#' Color palette to be used for the colvar variable, i.e. the cell density.
#' If colimage is NULL (default) and colvar is specified, then a 
#' red-magenta-white colorscheme will be used.
#' 
#' @param colkeypersp
#' A logical, or a list (default) with parameters for the color 
#' key (legend). List parameters should be one of side, plot, 
#' length, width, dist, shift, addlines, col.clab, cex.clab, 
#' side.clab, line.clab, adj.clab, font.clab and the axis parameters 
#' at, labels, tick, line, pos, outer, font, lty, lwd, lwd.ticks, 
#' col.box, col.axis, col.ticks, hadj, padj, cex.axis, mgp, tck, 
#' tcl, las. 
#' The defaults for the parameters are length = 0.2, width = 0.4, 
#' shift = 0.15, cex.axis = 0.6, cex.clab = 0.65
#' The default is to draw the color key on side = 4, i.e. in the right 
#' margin. 
#' If colkeypersp = NULL then a color key will be added only if 
#' col is a vector. Setting colkeypersp = list(plot = FALSE) will create 
#' room for the color key without drawing it. if colkeypersp = FALSE, 
#' no color key legend will be added.
#' See more details in ?plot3D::persp3D
#' 
#' @param colkeyimage
#' A logical, or a list (default) with parameters for the color 
#' key (legend). List parameters should be one of side, plot, 
#' length, width, dist, shift, addlines, col.clab, cex.clab, 
#' side.clab, line.clab, adj.clab, font.clab and the axis parameters 
#' at, labels, tick, line, pos, outer, font, lty, lwd, lwd.ticks, 
#' col.box, col.axis, col.ticks, hadj, padj, cex.axis, mgp, tck, 
#' tcl, las. 
#' The defaults for the parameters are length = 0.2, width = 0.4, 
#' shift = -0.15, cex.axis = 0.6, cex.clab = 0.65
#' The default is to draw the color key on side = 4, i.e. in the right 
#' margin. 
#' If colkeyimage = NULL then a color key will be added only if 
#' col is a vector. Setting colkeyimage = list(plot = FALSE) will create 
#' room for the color key without drawing it. if colkeyimage = FALSE, 
#' no color key legend will be added.
#' See more details in ?plot3D::image3D
#' 
#' @param lighting
#' A logical. If TRUE, the facets will be illuminated, and colors may 
#' appear more bright. Default is FALSE
#' 
#' @param lphi
#' if finite values are specified for lphi, the surface is 
#' shaded as though it was being illuminated from the direction 
#' specified by colatitude lphi
#' 
#' @param bty
#' The type of the box ("b" or "f"), the default only 
#' drawing background panels
#' 
#' @param scale_z
#' A logical, whether to scale the density z value or not
#' 
#' @param PDF
#' A logical. Output figure via pdf file or not, default is TRUE
#' 
#' @return Integration.l
#' A list contains input information and a dimention reduced matrix
#' by tSNE
#' 
#' @return A pdf file contains the generated figures
#' 
#' @details 
#' Density based visualization tool to generate figures that give the cell 
#' density distribution against potency states
#' 
#' @references 
#' Teschendorff AE, Tariq Enver. 
#' \emph{Single-cell entropy for accurate estimation of differentiation 
#' potency from a cell’s transcriptome.}
#' Nature communications 8 (2017): 15599.
#' doi:\href{https://doi.org/10.1038/ncomms15599}{
#' 10.1038/ncomms15599}.
#' 
#' @examples 
#' data(tsne.o)
#' data(potS.v)
#' data(SR4.v)
#' scent.o <- list(potS = potS.v)
#' potencyInfer.o <- list(potencyInfer.l = scent.o, SR = SR4.v)
#' 
#' LandSR.o <- Plot_LandSR(potencyInfer.o, reducedMatrix = tsne.o, reduceDim = FALSE, PDF = FALSE)
#' 
#' @import Rtsne
#' @import MASS
#' @import plot3D
#' @import irlba
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom graphics par
#' @importFrom grDevices colorRampPalette
#' @importFrom marray maPalette
#' @export
#'    
Plot_LandSR <- function(Integration.l,
                        reducedMatrix = NULL,
                        num_dim = 50,
                        max_components = 2,
                        num_grid = 50,
                        reduceDim = TRUE,
                        mean_cutoff = 1,
                        sd_cutoff = 0.5,
                        phi = 20,
                        theta = 20,
                        colpersp = NULL,
                        colimage = NULL,
                        colkeypersp = list(length = 0.2, width = 0.4, shift = 0.15, cex.axis = 0.6, cex.clab = 0.65),
                        colkeyimage = list(length = 0.2, width = 0.4, shift = -0.15, cex.axis = 0.6, cex.clab = 0.65),
                        lighting = FALSE,
                        lphi = 90,
                        bty = "b",
                        scale_z = TRUE,
                        PDF = TRUE)
{
    ### Reduce Dimension via tSNE method
    if (reduceDim == TRUE) {
        if (is.null(Integration.l$tSNE.mat)) {
            sd.v <- apply(Integration.l$expMC, 1, sd)
            mean.v <- apply(Integration.l$expMC, 1, mean)
            selG.idx <- intersect(which(mean.v > mean_cutoff),
                                  which(sd.v > sd_cutoff));
            
            irlba_res <- irlba::prcomp_irlba(t(Integration.l$expMC[selG.idx ,]), 
                                             n = num_dim, 
                                             center = TRUE)
            irlba_pca_res <- irlba_res$x
            topDim_pca <- irlba_pca_res
            tsne_res <- Rtsne::Rtsne(as.matrix(topDim_pca), dims = max_components, 
                                     pca = FALSE)
            reducedMatrix <- tsne_res$Y[, 1:max_components]
        }else{
            reducedMatrix <- Integration.l$tSNE.mat
        }
        component1.v <- reducedMatrix[, 1]
        component2.v <- reducedMatrix[, 2]
    }else{
        if (is.null(reducedMatrix)) {
            stop("Please input reduced dimension matrix of first two components!")
        }
        temp <- dim(reducedMatrix)
        if (temp[1] < temp[2]) {
            reducedMatrix <- t(reducedMatrix)
        }
        component1.v <- reducedMatrix[, 1]
        component2.v <- reducedMatrix[, 2]
    }
    
    Integration.l$tSNE.m <- reducedMatrix
    
    ### Calculate Cell Density
    CellDensity.o <- MASS::kde2d(x = component1.v,y = component2.v,n = num_grid)
    
    inferPotency.o <- Integration.l$potencyInfer.l
    potency_state.v <- inferPotency.o$potS
    idx.p <- length(unique(potency_state.v))
    
    if (idx.p > 9) {
        warning("Too much potency states, some potency states may have little cells!")
    }
    
    num_row <- ceiling(idx.p/3)
    num_col <- min(3, idx.p)
    
    ### select color
    if (is.null(colpersp)) {
        colorpersp.v <- colorRampPalette(c("red", "orange", "yellow2", "forestgreen", "blue"))(idx.p * 2)
    }
    if (is.null(colimage)) {
        colimage <- maPalette(low="white",mid="magenta",high="red",k=10)
    }
    
    ### Output via PDF or not
    if (PDF == TRUE) {
        pdf("PotencyComparison.pdf")
        par(mar=c(2,2,2,2))
        par(mfrow=c(num_row, num_col))
        color.id1 <- 1
        color.id2 <- 2
        for (i in 0 : (idx.p - 1)) {
            potency.o <- MASS::kde2d(x = component1.v[potency_state.v == (idx.p - i)],
                                     y = component2.v[potency_state.v == (idx.p - i)],
                                     n = 2*num_grid)
            if (scale_z == TRUE) {
                potency.o$z <- potency.o$z/max(potency.o$z)
            }
            xlim.v <- c(min(range(CellDensity.o$x), range(potency.o$x)), max(range(CellDensity.o$x), range(potency.o$x)))
            ylim.v <- c(min(range(CellDensity.o$y), range(potency.o$y)), max(range(CellDensity.o$y), range(potency.o$y)))
            
            if (is.null(colpersp)) {
                persp3D(x = potency.o$x, y = potency.o$y, z = (potency.o$z + (max(potency.o$z) * 0.5)) - (max(potency.o$z) * (0.2 * i)), 
                        xlim = xlim.v, ylim = ylim.v, zlim=c(-max(potency.o$z),max(potency.o$z)),
                        phi = phi, theta = theta, col = maPalette(low="lightgray", mid=colorpersp.v[color.id2], high=colorpersp.v[color.id1], k=10), 
                        colkey = colkeypersp, lighting = lighting, lphi = lphi, clab = c("","Cell Density",paste0("PS", (i+1))), bty = bty, plot = TRUE, xlab="tSNE1",ylab="tSNE2")
            }else{
                persp3D(x = potency.o$x, y = potency.o$y, z = (potency.o$z + (max(potency.o$z) * 0.5)) - (max(potency.o$z) * (0.2 * i)), 
                        xlim = xlim.v, ylim = ylim.v, zlim=c(-max(potency.o$z),max(potency.o$z)),
                        phi = phi, theta = theta, col = colpersp, colkey = colkeypersp, lighting = lighting, lphi = lphi, clab = c("","Cell Density",paste0("PS", (i+1))), bty = bty, plot = TRUE, xlab="tSNE1",ylab="tSNE2")
            }
            image3D(x = CellDensity.o$x, y = CellDensity.o$y, z = -max(potency.o$z), xlim = xlim.v, ylim = ylim.v,
                    colvar = CellDensity.o$z, col = colimage, colkey = colkeyimage, clab = c("","Cell Density","All"), add = TRUE, plot = TRUE)
            color.id1 <- color.id1 + 2
            color.id2 <- color.id2 + 2
        }
        dev.off()
    }else{
        par(mfrow=c(num_row, num_col))
        color.id1 <- 1
        color.id2 <- 2
        for (i in 0 : (idx.p - 1)) {
            potency.o <- MASS::kde2d(x = component1.v[potency_state.v == (idx.p - i)],
                                     y = component2.v[potency_state.v == (idx.p - i)],
                                     n = 2*num_grid)
            if (scale_z == TRUE) {
                potency.o$z <- potency.o$z/max(potency.o$z)
            }
            xlim.v <- c(min(range(CellDensity.o$x), range(potency.o$x)), max(range(CellDensity.o$x), range(potency.o$x)))
            ylim.v <- c(min(range(CellDensity.o$y), range(potency.o$y)), max(range(CellDensity.o$y), range(potency.o$y)))
            if (is.null(colpersp)) {
                persp3D(x = potency.o$x, y = potency.o$y, z = (potency.o$z + (max(potency.o$z) * 0.5)) - (max(potency.o$z) * (0.2 * i)), 
                        xlim = xlim.v, ylim = ylim.v, zlim=c(-max(potency.o$z),max(potency.o$z)),
                        phi = phi, theta = theta, col = maPalette(low="lightgray", mid=colorpersp.v[color.id2], high=colorpersp.v[color.id1], k=10), 
                        colkey = colkeypersp, lighting = lighting, lphi = lphi, clab = c("","Cell Density",paste0("PS", (i+1))), bty = bty, plot = TRUE, xlab="tSNE1",ylab="tSNE2")
            }else{
                persp3D(x = potency.o$x, y = potency.o$y, z = (potency.o$z + (max(potency.o$z) * 0.5)) - (max(potency.o$z) * (0.2 * i)), 
                        xlim = xlim.v, ylim = ylim.v, zlim=c(-max(potency.o$z),max(potency.o$z)),
                        phi = phi, theta = theta, col = colpersp, colkey = colkeypersp, lighting = lighting, lphi = lphi, clab = c("","Cell Density",paste0("PS", (i+1))), bty = bty, plot = TRUE, xlab="tSNE1",ylab="tSNE2")
            }
            image3D(x = CellDensity.o$x, y = CellDensity.o$y, z = -max(potency.o$z), xlim = xlim.v, ylim = ylim.v,
                    colvar = CellDensity.o$z, col = colimage, colkey = colkeyimage, clab = c("","Cell Density","All"), add = TRUE, plot = TRUE)
            color.id1 <- color.id1 + 2
            color.id2 <- color.id2 + 2
        }
    }
    
    return(Integration.l)
}