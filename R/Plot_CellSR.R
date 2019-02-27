#' @title 
#' Plot Cell density with SR values
#' 
#' @aliases Plot_CellSR
#'  
#' @description 
#' Density based visualization tool to generate figures that give the cell 
#' density distribution with SR value distribution
#' 
#' @param Integration.l
#' Typically, it is the output from \code{InferLandmark} function
#' 
#' @param coordinates
#' Optional. The previous reduced dimension coordinates, with rows lalabeling cells 
#' and two colums labeling reduced dimensions
#' 
#' @param num_grid
#' Number of grid points in each direction. Can be scalar or a length-2 
#' integer vector
#' 
#' @param phi
#' The angles defining the viewing direction. phi gives the colatitude
#' 
#' @param theta
#' The angles defining the viewing direction. theta gives the 
#' azimuthal direction
#' 
#' @param colpersp
#' Color palette to be used for the colvar variable, i.e. the cell density.
#' If colpersp is NULL (default) and colvar is specified, then a 
#' red-magenta-lightgray colorscheme will be used
#' 
#' @param colimage
#' Color palette to be used for the colvar variable, i.e. the SR values.
#' If colimage is NULL (default) and colvar is specified, then a 
#' blue-lightblue-white colorscheme will be used.
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
#' The type of the box ("b" or "f"), the default ("b") only 
#' drawing background panels
#' 
#' @param PDF
#' A logical. Output figure via pdf file or not, default is TRUE
#' 
#' @return A pdf file contains the generated figures
#' 
#' @return Integration.l
#' If coordinates provided, it will be a list object integrates the input
#' list and coordinates
#' 
#' @details 
#' Density based visualization tool to generate figures that give the cell 
#' density distribution against SR value distribution
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
#' potS.v <- 4 - potS.v
#' InferLandmark.o <- list(potencyState = potS.v, SR = SR4.v, coordinates = tsne.o)
#' 
#' CellSR.o <- Plot_CellSR(InferLandmark.o, PDF = FALSE)
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
#' @importFrom DelayedArray isEmpty
#' @export
#'  
Plot_CellSR <- function(Integration.l,
                        coordinates = NULL,
                        num_grid = 50,
                        phi = 20,
                        theta = 20,
                        colpersp = NULL,
                        colimage = NULL,
                        colkeypersp = list(length = 0.2, width = 0.4, shift = 0.2, cex.axis = 0.6, cex.clab = 0.7),
                        colkeyimage = list(length = 0.2, width = 0.4, shift = -0.2, cex.axis = 0.6, cex.clab = 0.7),
                        lighting = FALSE,
                        lphi = 90,
                        bty = "b",
                        PDF = TRUE)
{
    ### Reduce Dimension via tSNE method
    # if (reduceDim == TRUE) {
    #     if (is.null(Integration.l$tSNE.mat)) {
    #         sd.v <- apply(Integration.l$expMC, 1, sd)
    #         mean.v <- apply(Integration.l$expMC, 1, mean)
    #         selG.idx <- intersect(which(mean.v > mean_cutoff),
    #                               which(sd.v > sd_cutoff));
    #         
    #         irlba_res <- irlba::prcomp_irlba(t(Integration.l$expMC[selG.idx ,]), 
    #                                          n = num_dim, 
    #                                          center = TRUE)
    #         irlba_pca_res <- irlba_res$x
    #         topDim_pca <- irlba_pca_res
    #         tsne_res <- Rtsne::Rtsne(as.matrix(topDim_pca), dims = max_components, 
    #                                  pca = FALSE)
    #         coordinates <- tsne_res$Y[, 1:max_components]
    #     }else{
    #         coordinates <- Integration.l$tSNE.mat
    #     }
    #     component1.v <- coordinates[, 1]
    #     component2.v <- coordinates[, 2]
    # }else{
    #     if (is.null(coordinates)) {
    #         stop("Please input reduced dimension matrix of first two components!")
    #     }
    #     temp <- dim(coordinates)
    #     if (temp[1] < temp[2]) {
    #         coordinates <- t(coordinates)
    #     }
    #     component1.v <- coordinates[, 1]
    #     component2.v <- coordinates[, 2]
    # }
    # 
    # Integration.l$tSNE.m <- coordinates
    
    ### Get inherent reduced coordinates or input coordinates
    if (is.null(coordinates)) {
        coordinates <- Integration.l$coordinates
    }
    component1.v <- coordinates[, 1]
    component2.v <- coordinates[, 2]
    Integration.l$coordinates <- coordinates
    
    ### Calculate Cell Density
    CellDensity.o <- MASS::kde2d(x = component1.v,y = component2.v,n = num_grid)
    
    ### set SR
    inc.step1 <- floor((max(component1.v) - min(component1.v))/num_grid)
    inc.step2 <- floor((max(component2.v) - min(component2.v))/num_grid)
    SR_plot <- matrix(0, nrow = num_grid, ncol = num_grid)
    for (i in seq_len(num_grid)) {
        id.1 <- which(component1.v >= (min(component1.v) + (inc.step1 * (i-1))) 
                      & component1.v <= (min(component1.v) + (inc.step1 * i)))
        for (j in seq_len(num_grid)) {
            id.2 <- which(component2.v[id.1] >= (min(component2.v) + (inc.step2 * (j-1))) 
                          & component2.v[id.1] <= (min(component2.v) + (inc.step2 * j)))
            temp.id <- id.1[id.2]
            if (DelayedArray::isEmpty(temp.id)) {
                SR_plot[i,j] <- min(Integration.l$SR) * 0.85
                next()
            }
            SR_plot[i,j] <- max(Integration.l$SR[temp.id])
        }
    }
    
    ### select color
    if (is.null(colpersp)) {
        colpersp <- maPalette(low="white",mid="magenta",high="red", k=10)
    }
    if (is.null(colimage)) {
        colimage <- maPalette(low="gray",mid="lightblue",high="blue", k=15)
    }
    
    ### Output via PDF or not
    if (PDF == TRUE) {
        pdf("CellDensitySR.pdf")
        xlim.v <- c(min(range(CellDensity.o$x), range(component1.v)), max(range(CellDensity.o$x), range(component1.v)))
        ylim.v <- c(min(range(CellDensity.o$y), range(component2.v)), max(range(CellDensity.o$y), range(component2.v)))
        
        CellDensity.o$z <- CellDensity.o$z/max(CellDensity.o$z)
        
        persp3D(x = CellDensity.o$x, y = CellDensity.o$y, z = CellDensity.o$z, xlim = xlim.v, ylim = ylim.v, zlim=c(-max(CellDensity.o$z), max(CellDensity.o$z)),
                phi = phi, theta = theta, col = colpersp, colkey = colkeypersp, lighting = lighting, lphi = lphi, clab = c("","Cell Density", "All"), bty = bty, plot = TRUE, xlab="component 1",ylab="component 2")
        
        image3D(x = CellDensity.o$x, y = CellDensity.o$y, z = -max(CellDensity.o$z), xlim = xlim.v, ylim = ylim.v,
                colvar = SR_plot, col = colimage, colkey = colkeyimage, clab = c("","Potency","SR"), add = TRUE, plot = TRUE)
        dev.off()
    }else{
        xlim.v <- c(min(range(CellDensity.o$x), range(component1.v)), max(range(CellDensity.o$x), range(component1.v)))
        ylim.v <- c(min(range(CellDensity.o$y), range(component2.v)), max(range(CellDensity.o$y), range(component2.v)))
        
        CellDensity.o$z <- CellDensity.o$z/max(CellDensity.o$z)
        
        persp3D(x = CellDensity.o$x, y = CellDensity.o$y, z = CellDensity.o$z, xlim = xlim.v, ylim = ylim.v, zlim=c(-max(CellDensity.o$z), max(CellDensity.o$z)),
                phi = phi, theta = theta, col = colpersp, colkey = colkeypersp, lighting = lighting, lphi = lphi, clab = c("","Cell Density", "All"), bty = bty, plot = TRUE, xlab="component 1",ylab="component 2")
        
        image3D(x = CellDensity.o$x, y = CellDensity.o$y, z = -max(CellDensity.o$z), xlim = xlim.v, ylim = ylim.v,
                colvar = SR_plot, col = colimage, colkey = colkeyimage, clab = c("","Potency","SR"), add = TRUE, plot = TRUE)
    }
    
    return(Integration.l)
}
