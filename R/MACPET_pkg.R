#' An R-package for binding site analysis of ChIA-PET data.
#'
#' @description
#'The MACPET package can be used for binding site analysis for ChIA-PET data.
#'MACPET reads ChIA-PET data in BAM or SAM format and separates the data into Self-ligated,
#'Intra- and Inter-chromosomal PETs. Furthermore, MACPET breaks the genome into regions and
#' applies 2D mixture models for identifying candidate peaks/binding sites using
#' skewed generalized students-t distributions (SGT). It then uses a local poisson
#'  model for finding significant binding sites. MACPET is mainly written in C++,
#'  and it supports the BiocParallel package.
#'
#' @section MACPET main function:
#'  \code{\link{PeakCallerUlt}} runs a whole binding site analysis.
#'
#' @section MACPET classes:
#' \describe{
#' \item{\code{\linkS4class{PSelf}}}{S4 class for Self-ligated PETs.}
#' \item{\code{\linkS4class{PSFit}}}{S4 class for Self-ligated PETs after peak-calling.}
#' \item{\code{\linkS4class{PInter}}}{S4 class for Inter-chromosomal PETs.}
#' \item{\code{\linkS4class{PIntra}}}{S4 class for Intra-chromosomal PETs.}
#' }
#'
#' @section MACPET methods:
#' \describe{
#' \item{\code{\link{plot}}}{Method for plotting different objects.}
#' \item{\code{\link{summary}}}{Method for summarizing different objects.}
#' \item{\code{\link{AnalysisStatistics}}}{Prints summary of multiple
#' objects.}
#' \item{\code{\link{TagsToGInteractions}}}{Method for converting Tags to
#' \code{\linkS4class{GInteractions}} class.}
#' \item{\code{\link{PeaksToGRanges}}}{Method for converting peaks to
#' \code{\linkS4class{GRanges}} class.}
#' \item{\code{\link{exportPeaks}}}{Method for exporting peaks in cvs file format.}
#' \item{\code{\link{ConvertToPSelf}}}{Method for converting a GInteractions class of
#' Self-ligated PETs to object of \code{\linkS4class{PSelf}} class.}
#' \item{\code{\link{PeaksToNarrowPeak}}}{Method for converting peaks to narrowPeak (BED) format for use in
#' interaction analysis using the MANGO algorithm.}
#'}
#' @section MACPET sample data:
#' \describe{
#' \item{\code{\link{SampleChIAPETData.bam}}}{Sample ChIA-PET data.}
#' \item{\code{\link{pinterData.rda}}}{Sample \code{\linkS4class{PInter}} data.}
#' \item{\code{\link{pintraData.rda}}}{Sample \code{\linkS4class{PIntra}} data.}
#' \item{\code{\link{pselfData.rda}}}{Sample \code{\linkS4class{PSelf}} data.}
#' \item{\code{\link{psfitData.rda}}}{Sample \code{\linkS4class{PSFit}} data.}
#'}
#'
#' @docType package
#' @name MACPET
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#--#' make reference on the article:
#' @references
#' Vardaxis I, Drabløs F, Rye M and Lindqvist BH (2018). \emph{Model-based Analysis for ChIA-PET (MACPET)}.
#' To be published.
#'
#' @import Rcpp
#' @importFrom Rcpp evalCpp sourceCpp
#' @useDynLib MACPET

NULL
