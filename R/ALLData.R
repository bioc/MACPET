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
#' @seealso \code{\link{MACPET_pinterData.rda}}, \code{\link{MACPET_pintraData.rda}},
#' \code{\link{MACPET_pselfData.rda}},
#' \code{\link{MACPET_psfitData.rda}}, \code{\link{SampleChIAPETDataRead_1.bam}},
#'  \code{\link{SampleChIAPETDataRead_2.bam}}
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
#-------------------Main data read 1:
#' @title First reads from a subset of ChIA-PET data in \code{\link{SampleChIAPETData.bam}}
#' @name SampleChIAPETDataRead_1.bam
#'
#' @description First reads from a subset of ChIA-PET data in \code{\link{SampleChIAPETData.bam}}:
#' \describe{
#' \item{Target: }{ESR1}
#' \item{Biosample summary: }{Homo sapiens MCF-7}
#' \item{GEO: }{GSM970212}
#' }
#'
#' @seealso \code{\link{SampleChIAPETData.bam}}, \code{\link{ConvertToPE_BAM}}
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
#-------------------Main data read 2:
#' @title Second reads from a subset of ChIA-PET data in \code{\link{SampleChIAPETData.bam}}
#' @name SampleChIAPETDataRead_2.bam
#'
#' @description Second reads from a subset of ChIA-PET data in \code{\link{SampleChIAPETData.bam}}:
#' \describe{
#' \item{Target: }{ESR1}
#' \item{Biosample summary: }{Homo sapiens MCF-7}
#' \item{GEO: }{GSM970212}
#' }
#'
#' @seealso \code{\link{SampleChIAPETData.bam}},  \code{\link{ConvertToPE_BAM}}
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
#' @name MACPET_pinterData.rda
#' @description Inter-chromosomal PETs data from
#' ESR1 ChIA-PET subset data on human MCF-7.
#'
#' @details \code{MACPET_pinterData} is produced by the \code{\link{MACPETUlt}}
#' function at Stage 2 and it contains the Inter-chromosomal PETs of the sample data.
#'
#' @author
#' \describe{
#' \item{Main data creators}{Yijun Ruan, GIS, 2012-05-24}
#' \item{MACPET_pinterData.rda creator}{Ioannis Vardaxis,
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
#' @name MACPET_pintraData.rda
#' @description Intra-chromosomal PETs data from
#' ESR1 ChIA-PET subset data on human MCF-7.
#'
#' @details \code{MACPET_pintraData} is produced by the \code{\link{MACPETUlt}}
#' function at Stage 2 and it contains the Intra-chromosomal PETs of the sample data.
#'
#' @author
#' \describe{
#' \item{Main data creators}{Yijun Ruan, GIS, 2012-05-24}
#' \item{MACPET_pintraData.rda creator}{Ioannis Vardaxis,
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
#' @name MACPET_pselfData.rda
#' @description Self-ligated PETs data from
#' ESR1 ChIA-PET subset data on human MCF-7.
#'
#' @details \code{MACPET_pselfData} is produced by the \code{\link{MACPETUlt}}
#' function at Stage 2 and it contains the Self-ligated PETs of the sample data.
#'
#' @author
#' \describe{
#' \item{Main data creators}{Yijun Ruan, GIS, 2012-05-24}
#' \item{MACPET_pselfData.rda creator}{Ioannis Vardaxis,
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
#' @name MACPET_psfitData.rda
#' @description Self-ligated PETs data from
#' ESR1 ChIA-PET subset data on human MCF-7.
#'
#' @details \code{MACPET_psfitData} is produced by the \code{\link{MACPETUlt}}
#' function at Stage 3 and it contains the self-ligated PETs of the sample data after
#' calling for candidate peaks.
#'
#' @author
#' \describe{
#' \item{Main data creators}{Yijun Ruan, GIS, 2012-05-24}
#' \item{MACPET_psfitData.rda creator}{Ioannis Vardaxis,
#' \email{ioannis.vardaxis@@ntnu.no}}
#' }
#'
#' @format rda object of \code{\linkS4class{PSFit}} class.
#' @docType data
#' @seealso \code{\link{SampleChIAPETData.bam}},
#' \code{\linkS4class{PSFit}}
#' @keywords data
NULL
#-------------------GenomeMap Data:
#' @title Genomic interactions from ChIA-PET data
#' @name MACPET_GenomeMapData.rda
#' @description Genomic Interactions data from
#' ESR1 ChIA-PET subset data on human MCF-7.
#'
#' @details \code{MACPET_GenomeMapData} is produced by the \code{\link{MACPETUlt}}
#' function at Stage 4 and it contains all the genomic interactions between the peaks
#' in the data.
#'
#' @author
#' \describe{
#' \item{Main data creators}{Yijun Ruan, GIS, 2012-05-24}
#' \item{MACPET_GenomeMapData.rda creator}{Ioannis Vardaxis,
#' \email{ioannis.vardaxis@@ntnu.no}}
#' }
#'
#' @format rda object of \code{\linkS4class{GenomeMap}} class.
#' @docType data
#' @seealso  \code{\linkS4class{GenomeMap}}
#' @keywords data
NULL
