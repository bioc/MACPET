#-------------------Main data:
#' @title Subset sample of ChIA-PET data
#' @name SampleChIAPETData.bam
#'
#' @description A subset of ChIA-PET data:
#' \describe{
#' \item{Target: }{ESR1}
#' \item{Biosample summary: }{Homo sapiens MCF-7}
#' \item{GEO: }{GSM970212}
#' }
#'
#' @seealso \code{\link{pselfData.rda}}, \code{\link{pintraData.rda}},
#' \code{\link{pinterData.rda}},
#' \code{\link{psfitData.rda}}
#'
#' @author Yijun Ruan, GIS, 2012-05-24 (main data creators)
#' @format A BAM file.
#' @references
#' Consortium EP (2012) \emph{An integrated encyclopedia of DNA elements in the human genome.}.
#'  Nature, 489(7414), pp. 57–74. \url{http://dx.doi.org/10.1038/nature11247}.
#' @docType data
#' @source \url{https://www.encodeproject.org/experiments/ENCSR000BZZ/}
#' @keywords data
NULL

#-------------------Inter data:
#' @title Inter-chromosomal PETs from ChIA-PET data
#' @name pinterData.rda
#' @description Inter-chromosomal PETs data from
#' ESR1 ChIA-PET subset data on human MCF-7.
#'
#' @details \code{pinterData} is produced by the \code{\link{PeakCallerUlt}}
#' function at Stage=0 and it contains the Inter-chromosomal PETs of the sample data.
#'
#' @author
#' \describe{
#' \item{Main data creators}{Yijun Ruan, GIS, 2012-05-24}
#' \item{pinterData.rda creator}{Ioannis Vardaxis,
#' \email{ioannis.vardaxis@@ntnu.no}}
#' }
#'
#' @format rda object of \code{\linkS4class{PInter}} class.
#' @docType data
#' @seealso \code{\link{SampleChIAPETData.bam}}, \code{\linkS4class{PInter}}
#' @keywords data
NULL

#-------------------Intra data:
#' @title Intra-chromosomal PETs from ChIA-PET data
#' @name pintraData.rda
#' @description Intra-chromosomal PETs data from
#' ESR1 ChIA-PET subset data on human MCF-7.
#'
#' @details \code{pintraData} is produced by the \code{\link{PeakCallerUlt}}
#' function at Stage 0 and it contains the Intra-chromosomal PETs of the sample data.
#'
#' @author
#' \describe{
#' \item{Main data creators}{Yijun Ruan, GIS, 2012-05-24}
#' \item{pintraData.rda creator}{Ioannis Vardaxis,
#'  \email{ioannis.vardaxis@@ntnu.no}}
#' }
#'
#' @format rda object of \code{\linkS4class{PIntra}} class.
#' @docType data
#' @seealso \code{\link{SampleChIAPETData.bam}}, \code{\linkS4class{PIntra}}
#' @keywords data
NULL

#-------------------Self data:
#' @title Self-ligated PETs from ChIA-PET data
#' @name pselfData.rda
#' @description Self-ligated PETs data from
#' ESR1 ChIA-PET subset data on human MCF-7.
#'
#' @details \code{pselfData} is produced by the \code{\link{PeakCallerUlt}}
#' function at Stage 0 and it contains the Self-ligated PETs of the sample data.
#'
#' @author
#' \describe{
#' \item{Main data creators}{Yijun Ruan, GIS, 2012-05-24}
#' \item{pselfData.rda creator}{Ioannis Vardaxis,
#' \email{ioannis.vardaxis@@ntnu.no}}
#' }
#'
#' @format rda object of \code{\linkS4class{PSelf}} class.
#' @docType data
#' @seealso \code{\link{SampleChIAPETData.bam}}, \code{\linkS4class{PSelf}}
#' @keywords data
NULL


#-------------------Self data FIT:
#' @title Self-ligated PETs from ChIA-PET data
#' @name psfitData.rda
#' @description Self-ligated PETs data from
#' ESR1 ChIA-PET subset data on human MCF-7.
#'
#' @details \code{psfitData} is produced by the \code{\link{PeakCallerUlt}}
#' function at Stage 1 and it contains the self-ligated PETs of the sample data after
#' calling for candidate peaks.
#'
#' @author
#' \describe{
#' \item{Main data creators}{Yijun Ruan, GIS, 2012-05-24}
#' \item{psfitData.rda creator}{Ioannis Vardaxis,
#' \email{ioannis.vardaxis@@ntnu.no}}
#' }
#'
#' @format rda object of \code{\linkS4class{PSFit}} class.
#' @docType data
#' @seealso \code{\link{SampleChIAPETData.bam}},
#' \code{\linkS4class{PSFit}}
#' @keywords data
NULL
