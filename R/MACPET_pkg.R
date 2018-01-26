#' An R-package for binding site analysis of ChIA-PET data.
#'
#' @description
#' The MACPET package can be used for general analysis of paired-end (PET) data like ChIA-PET.
#' MACPET currently implements the following four stages: Linker filtering (stage 0), mapping to
#' the reference genome (stage 1), PET classification (stage 2) and peak-calling (stage 3).
#' All of the MACPET stages can be run at once, or separately.
#' In stage 0, MACPET identifies the linkers in the fastq files and classifies the
#' reads as usable, chimeric or ambiguous. Usable reads are considered in the
#' subsequent stages. In stage 1, MACPET maps the usable reads to the reference
#' genome using \code{\link[Rbowtie:bowtie]{bowtie}} and produces a paired-end BAM file. This BAM file is further used in
#' stage 2 to classify the PETs as self-ligated/intra- or inter-chromosomal.
#' Self-ligated PETs are used in stage 3 for the identification of significant peaks.
#' In stage 3, MACPET segments the genome into regions and
#' applies 2D mixture models for identifying candidate peaks using
#' skewed generalized students-t distributions (SGT). It then uses a local poisson
#' model for finding significant binding sites. MACPET is mainly written in C++,
#' and it supports the BiocParallel package.
#'
#' @section MACPET main function:
#'  \code{\link{MACPETUlt}} runs the whole analysis at once.
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
#' @section MACPET supplementary functions:
#' \describe{
#' \item{\code{\link{ConvertToPE_BAM}}}{Function for converting two separate BAM files into one paired-end BAM file.}
#'  \item{\code{\link{AnalysisStatistics}}}{Prints summary of multiple
#' objects.}
#' }
#' @section MACPET sample data:
#' \describe{
#' \item{\code{\link{SampleChIAPETData.bam}}}{Sample ChIA-PET data.}
#' \item{\code{\link{SampleChIAPETDataRead_1.bam}}}{First reads from the sample ChIA-PET data.}
#' \item{\code{\link{SampleChIAPETDataRead_2.bam}}}{Second reads from the sample ChIA-PET data.}
#' \item{\code{\link{MACPET_pinterData.rda}}}{Sample \code{\linkS4class{PInter}} data.}
#' \item{\code{\link{MACPET_pintraData.rda}}}{Sample \code{\linkS4class{PIntra}} data.}
#' \item{\code{\link{MACPET_pselfData.rda}}}{Sample \code{\linkS4class{PSelf}} data.}
#' \item{\code{\link{MACPET_psfitData.rda}}}{Sample \code{\linkS4class{PSFit}} data.}
#'}
#'
#' @docType package
#' @name MACPET
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#--#' make reference on the article:
#' @references
#' Vardaxis I, Drabløs F, Rye M and Lindqvist BH (2018). \emph{MACPET: Model-based Analysis for ChIA-PET}.
#' To be published.
#'
#' @import Rcpp
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom methods is
#' @useDynLib MACPET
NULL
