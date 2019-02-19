#' Normalized example scRNA-seq data matrix
#' 
#' Normalized example scRNA-seq data matrix of 50 human 
#' embryonic stem cells(hESC) and 50 mesoderm progenitor cells(EC)
#' from Chu et al. 2016.
#'
#' \itemize{
#'   \item hESC : human embryonic stem cells(n=50)
#'   \item EC : mesoderm progenitor cells(n=50)
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name Example.m
#' @usage data(Example.m)
#' @format A matrix with 18935 rows and 100 columns
#' @references 
#' Teschendorff AE, Tariq Enver. 
#' \emph{Single-cell entropy for accurate estimation of differentiation 
#' potency from a cell’s transcriptome.}
#' Nature communications 8 (2017): 15599.
#' doi:\href{https://doi.org/10.1038/ncomms15599}{
#' 10.1038/ncomms15599}.
#' 
#' Chu Li-Fang, Ning Leng, Jue Zhang, Zhonggang Hou, Daniel Mamott, 
#' David T. Vereide, Jeea Choi, Christina Kendziorski, Ron Stewart and James A. Thomson
#' \emph{Single-cell RNA-seq reveals novel regulators of human embryonic 
#' stem cell differentiation to definitive endoderm.}
#' Genome biology 17.1 (2016): 173.
#' doi:\href{https://doi.org/10.1186/s13059-016-1033-x}{
#' 10.1186/s13059-016-1033-x}.
#' 
NULL

#' Protein-protein interaction network 13Jun12
#'
#' This protein-protein interaction network is derived from Pathway 
#' Commons (www.pathwaycommons.org) (version Jun. 2012), which is 
#' an integrated resource collating together PPIs from several 
#' distinct sources.
#'
#' \itemize{
#'   \item 121268 : Gene EntrzID 
#'   \item 2193 : Gene EntrzID 
#'   \item 79902 : Gene EntrzID
#'   \item ...
#' }
#'
#' @docType data
#' @keywords network
#' @name net13Jun12.m
#' @usage data(net13Jun12.m)
#' @format A matrix with 8434 rows and 8434 columns
#' @references 
#' Teschendorff AE, Tariq Enver. 
#' \emph{Single-cell entropy for accurate estimation of differentiation 
#' potency from a cell’s transcriptome.}
#' Nature communications 8 (2017): 15599.
#' doi:\href{https://doi.org/10.1038/ncomms15599}{
#' 10.1038/ncomms15599}.
#' 
NULL

#' Protein-protein interaction network 17Jan16
#'
#' This protein-protein interaction network is derived from Pathway 
#' Commons (www.pathwaycommons.org) (version Jan. 2016), which is 
#' an integrated resource collating together PPIs from several 
#' distinct sources.
#'
#' \itemize{
#'   \item 121268 : Gene EntrzID
#'   \item 2193 : Gene EntrzID
#'   \item 79902 : Gene EntrzID
#'   \item ...
#' }
#'
#' @docType data
#' @keywords network
#' @name net17Jan16.m
#' @usage data(net17Jan16.m)
#' @format A matrix with 11751 rows and 11751 columns
#' @references 
#' Teschendorff AE, Tariq Enver. 
#' \emph{Single-cell entropy for accurate estimation of differentiation 
#' potency from a cell’s transcriptome.}
#' Nature communications 8 (2017): 15599.
#' doi:\href{https://doi.org/10.1038/ncomms15599}{
#' 10.1038/ncomms15599}.
#' 
NULL

#' Phenotype information of the example scRNA-seq data
#' 
#' Phenotype information for the example scRNA-seq data matrix of 50 human 
#' embryonic stem cells(hESC) and 50 mesoderm progenitor cells(EC)
#' from Chu et al. 2016.
#'
#' \itemize{
#'   \item hESC : human embryonic stem cells(n=50)
#'   \item EC : mesoderm progenitor cells(n=50)
#' }
#'
#' @docType data
#' @keywords datasets
#' @name phenoExample.v
#' @usage data(phenoExample.v)
#' @format A character vector includes 100 elements
#' @references 
#' Teschendorff AE, Tariq Enver. 
#' \emph{Single-cell entropy for accurate estimation of differentiation 
#' potency from a cell’s transcriptome.}
#' Nature communications 8 (2017): 15599.
#' doi:\href{https://doi.org/10.1038/ncomms15599}{
#' 10.1038/ncomms15599}.
#' 
#' Chu Li-Fang, Ning Leng, Jue Zhang, Zhonggang Hou, Daniel Mamott, 
#' David T. Vereide, Jeea Choi, Christina Kendziorski, Ron Stewart and James A. Thomson
#' \emph{Single-cell RNA-seq reveals novel regulators of human embryonic 
#' stem cell differentiation to definitive endoderm.}
#' Genome biology 17.1 (2016): 173.
#' doi:\href{https://doi.org/10.1186/s13059-016-1033-x}{
#' 10.1186/s13059-016-1033-x}. 
#' 
NULL

#' Phenotype information for Chu et al. scRNA-seq data
#' 
#' Phenotype information for the whole scRNA-seq data matrix 
#' from Chu et al. 2016.
#'
#' \itemize{
#'   \item hESC : human embryonic stem cells(n=374)
#'   \item NPC : ectoderm progenitor(n=173)
#'   \item DEC : endoderm progenitor(n=138)
#'   \item EC : mesoderm progenitor cells(n=105)
#'   \item HFF : human foreskin fibroblasts(n=159)
#'   \item TB : trophoblasts(n=69)
#' }
#'
#' @docType data
#' @keywords datasets
#' @name phenoscChu.v
#' @usage data(phenoscChu.v)
#' @format A character vector includes 1018 elements
#' @references 
#' Teschendorff AE, Tariq Enver. 
#' \emph{Single-cell entropy for accurate estimation of differentiation 
#' potency from a cell’s transcriptome.}
#' Nature communications 8 (2017): 15599.
#' doi:\href{https://doi.org/10.1038/ncomms15599}{
#' 10.1038/ncomms15599}.
#' 
#' Chu Li-Fang, Ning Leng, Jue Zhang, Zhonggang Hou, Daniel Mamott, 
#' David T. Vereide, Jeea Choi, Christina Kendziorski, Ron Stewart and James A. Thomson
#' \emph{Single-cell RNA-seq reveals novel regulators of human embryonic 
#' stem cell differentiation to definitive endoderm.}
#' Genome biology 17.1 (2016): 173.
#' doi:\href{https://doi.org/10.1186/s13059-016-1033-x}{
#' 10.1186/s13059-016-1033-x}. 
#' 
NULL

#' Potency states for example plot scRNA-seq data
#' 
#' Potency states infered from \code{PotencyInfer} function with 
#' example scRNA-seq data.
#'
#' \itemize{
#'   \item Potency States : numbers indicates potency levels
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name potS.v
#' @usage data(potS.v)
#' @format A numeric vector includes 3473 elements
#' @references 
#' Teschendorff AE, Tariq Enver. 
#' \emph{Single-cell entropy for accurate estimation of differentiation 
#' potency from a cell’s transcriptome.}
#' Nature communications 8 (2017): 15599.
#' doi:\href{https://doi.org/10.1038/ncomms15599}{
#' 10.1038/ncomms15599}.
#' 
#' 
NULL

#' Raw example scRNA-seq data matrix
#' 
#' Raw example scRNA-seq data matrix of 50 human 
#' embryonic stem cells(hESC) and 50 mesoderm progenitor cells(EC)
#' from Chu et al. 2016.
#'
#' \itemize{
#'   \item hESC : human embryonic stem cells(n=50)
#'   \item EC : mesoderm progenitor cells(n=50)
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name rawExample.m
#' @usage data(rawExample.m)
#' @format A matrix with 19097 rows(genes) and 100 columns(cells)
#' @references 
#' Teschendorff AE, Tariq Enver. 
#' \emph{Single-cell entropy for accurate estimation of differentiation 
#' potency from a cell’s transcriptome.}
#' Nature communications 8 (2017): 15599.
#' doi:\href{https://doi.org/10.1038/ncomms15599}{
#' 10.1038/ncomms15599}.
#' 
#' Chu Li-Fang, Ning Leng, Jue Zhang, Zhonggang Hou, Daniel Mamott, 
#' David T. Vereide, Jeea Choi, Christina Kendziorski, Ron Stewart and James A. Thomson
#' \emph{Single-cell RNA-seq reveals novel regulators of human embryonic 
#' stem cell differentiation to definitive endoderm.}
#' Genome biology 17.1 (2016): 173.
#' doi:\href{https://doi.org/10.1186/s13059-016-1033-x}{
#' 10.1186/s13059-016-1033-x}.
#' 
NULL

#' SR values for example scRNA-seq data
#' 
#' SR values calculated from \code{CompSRana} function with the input 
#' of example scRNA-seq data and network information.
#'
#' \itemize{
#'   \item SR value : Normalized Signaling Entropy Rate
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name SR.v
#' @usage data(SR.v)
#' @format A vector with 96 elements
#' @references 
#' Teschendorff AE, Tariq Enver. 
#' \emph{Single-cell entropy for accurate estimation of differentiation 
#' potency from a cell’s transcriptome.}
#' Nature communications 8 (2017): 15599.
#' doi:\href{https://doi.org/10.1038/ncomms15599}{
#' 10.1038/ncomms15599}.
#' 
NULL

#' SR values for example plot scRNA-seq data
#' 
#' SR values calculated from \code{CompSRana} function with the input 
#' of example plot scRNA-seq data and network information.
#'
#' \itemize{
#'   \item SR value : Normalized Signaling Entropy Rate
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name SR4.v
#' @usage data(SR4.v)
#' @format A vector with 3473 elements
#' @references 
#' Teschendorff AE, Tariq Enver. 
#' \emph{Single-cell entropy for accurate estimation of differentiation 
#' potency from a cell’s transcriptome.}
#' Nature communications 8 (2017): 15599.
#' doi:\href{https://doi.org/10.1038/ncomms15599}{
#' 10.1038/ncomms15599}.
#' 
NULL

#' Reduced dimension matrix for example plot scRNA-seq data
#' 
#' Dimension reduced via \code{tsne} function.
#'
#' \itemize{
#'   \item data projection values in two-dimension space
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name tsne.o
#' @usage data(tsne.o)
#' @format A numeric matrix with two colums and 3473 rows
#' @references 
#' Teschendorff AE, Tariq Enver. 
#' \emph{Single-cell entropy for accurate estimation of differentiation 
#' potency from a cell’s transcriptome.}
#' Nature communications 8 (2017): 15599.
#' doi:\href{https://doi.org/10.1038/ncomms15599}{
#' 10.1038/ncomms15599}.
#' 
#' 
NULL