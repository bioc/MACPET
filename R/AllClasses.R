#' @importFrom methods setClass new
#' @importClassesFrom InteractionSet GInteractions
#' @title PSelf S4 Class
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#' @references
#' Vardaxis I, Drabløs F, Rye M and Lindqvist BH (2018). \emph{Model-based Analysis for ChIA-PET (MACPET)}.
#' To be published.
#' @description  \code{\linkS4class{PSelf}} class in a S4 class which
#'  inherits from the
#' \code{\link[InteractionSet]{GInteractions}} class and it contains Self-ligated
#' PETs from ChIA-PET experiment. Furthermore it also contains the following in the
#' \code{\link[S4Vectors]{metadata}} field:
#'  \describe{
#'  \item{\code{Self_info}}{ A data.frame with the count statistics for the total PETs in each chromosome.}
#'  \item{\code{SLmean}}{The mean size of the PETs in the data.}
#'  \item{\code{GenInfo}}{Information about the genome.}
#'  \item{\code{MaxSize}}{Maximum size of self-ligated PETs.}
#'  \item{\code{MinSize}}{Minimum size of self-ligated PETs.}
#'  }
#' @details \code{\linkS4class{PSelf}} class is created by
#' the \code{\link{PeakCallerUlt}}
#' function at Stage=0.
#' @export
#' @seealso \code{\link{AnalysisStatistics}}, \code{\link{plot}},
#' \code{\link{summary}},
#' \code{\link{PeakCallerUlt}}
PSelf=setClass("PSelf",contains="GInteractions")

#' @export
#' @title PSFit S4 Class
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#' @references
#' Vardaxis I, Drabløs F, Rye M and Lindqvist BH (2018). \emph{Model-based Analysis for ChIA-PET (MACPET)}.
#' To be published.
#' @description \code{\linkS4class{PSFit}} class in a S4 class which
#' inherits from the
#' \code{\link[InteractionSet]{GInteractions}} class and it contains Self-ligated
#' PETs from ChIA-PET experiment and information about the
#' genome of the data. Furthermore it also contains the following in the
#' \code{\link[S4Vectors]{metadata}} field:
#'  \describe{
#'  \item{\code{Self_info}}{Counts statistics for the total PETs,
#'  total regions and total Peaks in each chromosome.}
#'  \item{\code{SLmean}}{The mean size of the PETs in the data.}
#'  \item{\code{GenInfo}}{Information about the genome.}
#'  \item{\code{MaxSize}}{Maximum size of self-ligated PETs.}
#'  \item{\code{MinSize}}{Minimum size of self-ligated PETs.}
#'  \item{\code{Classification.Info}}{A matrix with information for the Data-row ID, region ID and Peak ID (0 represent noise)
#'  of each peak in the data.}
#'  \item{\code{Peaks.Info}}{Information for each peak found by the peak-calling algorithm:
#'   \describe{
#'   \item{\code{Chrom}}{The chromosome which the peak belongs to.}
#'   \item{\code{Region}}{The region which the peak belongs to.}
#'   \item{\code{Peak}}{The peak ID (a region might have more than one peaks).}
#'   \item{\code{Pets}}{Total PETs in the peak.}
#'   \item{\code{Peak.Summit}}{Summit of the peak.}
#'   \item{\code{Up.Summit}}{Summit of the left-stream PETs.}
#'   \item{\code{Down.Summit}}{Summit of the right-stream PETs.}
#'   \item{\code{CIQ.Up.start}}{Start of the 95 Quantile confidence interval for the left-stream PETs.}
#'   \item{\code{CIQ.Up.end}}{End of the 95 Quantile confidence interval for the left-stream PETs.}
#'   \item{\code{CIQ.Up.size}}{Size of the 95 Quantile confidence interval for the left-stream PETs.}
#'   \item{\code{CIQ.Down.start}}{Start of the 95 Quantile confidence interval for the right-stream PETs.}
#'   \item{\code{CIQ.Down.end}}{End of the 95 Quantile confidence interval for the right-stream PETs.}
#'   \item{\code{CIQ.Down.size}}{Size of the 95 Quantile confidence interval for the right-stream PETs.}
#'   \item{\code{CIQ.Peak.size}}{Size of the Peak based on the interval (CIQ.Up.start,CIQ.Down.end).}
#'   \item{\code{lambdaUp}}{The expected number of PETs in the left-stream Peak region by random chance.}
#'   \item{\code{FoldEnrichUp}}{Fold enrichment for the left-stream Peak region.}
#'   \item{\code{p.valueUp}}{p-value for the left-stream Peak region.}
#'   \item{\code{lambdaDown}}{The expected number of PETs in the right-stream Peak region by random chance.}
#'   \item{\code{FoldEnrichDown}}{Fold enrichment for the right-stream Peak region.}
#'   \item{\code{p.valueDown}}{p-value for the right-stream Peak region.}
#'   \item{\code{p.value}}{p-value for the Peak (p.valueUp*p.valueDown).}
#'   \item{\code{FDRUp}}{FDR correction for the left-stream Peak region.}
#'   \item{\code{FDRDown}}{FDR correction for the right-stream Peak region.}
#'   \item{\code{FDR}}{FDR correction for the Peak.}
#'   }
#'  }
#'  }
#' @details \code{\linkS4class{PSFit}} class is created by the
#' \code{\link{PeakCallerUlt}} function at Stage=1.
#' @seealso \code{\link{AnalysisStatistics}},\code{\link{plot}},
#' \code{\link{summary}}, \code{\link{PeakCallerUlt}},
#' \code{\link{exportPeaks}}, \code{\link{PeaksToGRanges}},
#' \code{\link{TagsToGInteractions}}, \code{\link{PeaksToNarrowPeak}}
PSFit=setClass("PSFit",contains="GInteractions")

#' @export
#' @title PInter S4 Class
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#' @references
#' Vardaxis I, Drabløs F, Rye M and Lindqvist BH (2018). \emph{Model-based Analysis for ChIA-PET (MACPET)}.
#' To be published.
#' @description \code{\linkS4class{PInter}} class in a S4 class which
#'  inherits from the
#' \code{\link[InteractionSet]{GInteractions}} class and it contains
#' Inter-chromosomal data.
#' @details \code{\linkS4class{PInter}} class is created by the
#'  \code{\link{PeakCallerUlt}}
#' function at Stage 0.
#' @seealso \code{\link{AnalysisStatistics}}, \code{\link{plot}},
#' \code{\link{summary}},
#' \code{\link{PeakCallerUlt}}
PInter=setClass("PInter",contains="GInteractions")

#' @export
#' @title PIntra Class
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#' @references
#' Vardaxis I, Drabløs F, Rye M and Lindqvist BH (2018). \emph{Model-based Analysis for ChIA-PET (MACPET)}.
#' To be published.
#' @description \code{\linkS4class{PIntra}} class in a S4 class which
#' inherits from the
#' \code{\link[InteractionSet]{GInteractions}} class and it contains
#' Intra-chromosomal data.
#' @details \code{\linkS4class{PIntra}} class is created by the
#' \code{\link{PeakCallerUlt}}
#' function at Stage 0.
#' @seealso \code{\link{AnalysisStatistics}}, \code{\link{plot}},
#' \code{\link{summary}},
#' \code{\link{PeakCallerUlt}}
PIntra=setClass("PIntra",contains="GInteractions")
