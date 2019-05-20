#' @title 
#' Visualize diffusion map
#' 
#' @aliases Plot_DiffusionMap
#'  
#' @description 
#' This function visualizes the diffusion map in 2D or 3D plots
#' 
#' @param Integration.l
#' Typically, it is the output from \code{DoDiffusionMap} function
#' 
#' @param color_by
#' Indicating the variable used for coloring. \code{SR}
#' means the plot is colored by SR values. \code{DPT}
#' means the plot is colored by diffusion pseudotime
#' 
#' @param dim
#' A numeric vector. Diffusin components order in the plot axes. And 
#' the sign of evrey entry indicates the direction of component.
#' Default is c(1, 2, 3)
#' 
#' @param phi
#' The angles defining the viewing direction. \code{phi} 
#' gives the colatitude. Default is 20
#' 
#' @param theta
#' The angles defining the viewing direction. \code{theta} 
#' gives the azimuthal direction. Default is 30
#' 
#' @param PDF
#' A logical. Output figure via pdf file or not, default is FALSE
#' 
#' @return A pdf file contains the generated figures
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
#' DoDM.o <- DoDiffusionMap(InferPotency.o)
#' Plot_DiffusionMap(DoDM.o)
#' }
#' 
#' @import ggplot2
#' @import plot3D
#' @importFrom destiny DPT
#' @export
#'     
Plot_DiffusionMap <- function(Integration.l,
                              dim = c(1, 2, 3),
                              color_by = c("SR", "DPT"),
                              phi = 20,
                              theta = 30,
                              PDF = FALSE){
  dm <- Integration.l$DM
  dms <- Integration.l$DMEigen
  DIM_dms <- dim(dms)[2]
  
  sign_dim <- sign(dim)
  dim <- abs(dim)
  
  DIMS <- max(dim)
  if (DIMS > DIM_dms) {
    stop("Diffusion map object does not contain enough dimensions, please reset dim argument!")
  }
  
  if (length(color_by) > 1) {
    print("Specified more than one type of color, using the first type instead!")
    color_by <- color_by[1]
  }
  
  if (color_by == "SR") {
    color.idx <- Integration.l$SR
    color.lab <- "SR"
    panel.text <- "Diffusion Map with SR values"
  }else{
    dpt <- destiny::DPT(dm, tips = Integration.l$root)
    color.idx <- dpt[["dpt"]]
    color.lab <- "DPT"
    panel.text <- "Diffusion Map with DPT estimation"
  }
  
  labs <- colnames(dms)[dim]
  DIMS <- length(dim)
  
  if (PDF) {
    pdf("DiffusionMap.pdf")
    par(mar = c(2,2,2,2))
    
    if (DIMS == 2) {
      g <- ggplot(dms, aes((sign_dim[1]*dms[, dim[1]]), (sign_dim[2]*dms[, dim[2]]), color = color.idx)) +
        geom_point() +
        xlab(labs[1]) +
        ylab(labs[2]) +
        labs(title = panel.text) +
        theme(legend.position = "right") +
        labs(color = color.lab)
      print(g)
    }else{
      points3D(x = (sign_dim[1]*dms[, dim[1]]), 
               y = (sign_dim[2]*dms[, dim[2]]), 
               z = (sign_dim[3]*dms[, dim[3]]), 
               colvar = color.idx, 
               phi = phi, theta = theta, 
               xlab = labs[1],
               ylab = labs[2],
               zlab = labs[3],
               clab = color.lab,
               pch = 20,
               main = panel.text)
    }
    
    dev.off()
  }else{
    if (DIMS == 2) {
      g <- ggplot(dms, aes((sign_dim[1]*dms[, dim[1]]), (sign_dim[2]*dms[, dim[2]]), color = color.idx)) +
        geom_point() +
        xlab(labs[1]) +
        ylab(labs[2]) +
        labs(title = panel.text) +
        theme(legend.position = "right") +
        labs(color = color.lab)
      print(g)
    }else{
      points3D(x = (sign_dim[1]*dms[, dim[1]]), 
               y = (sign_dim[2]*dms[, dim[2]]), 
               z = (sign_dim[3]*dms[, dim[3]]), 
               colvar = color.idx, 
               phi = phi, theta = theta, 
               xlab = labs[1],
               ylab = labs[2],
               zlab = labs[3],
               clab = color.lab,
               pch = 20,
               main = panel.text)
    }
  }
}





