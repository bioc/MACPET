#' @title Convert peaks to GRanges object
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#' @references
#' Vardaxis I, Drabløs F, Rye M and Lindqvist BH (2018). \emph{MACPET: Model-based Analysis for ChIA-PET}.
#' To be published.
#'
#' @description \code{PeaksToGRanges} converts peaks of an object of \code{\linkS4class{PSFit}}
#' class to \code{\linkS4class{GRanges}} object.
#' @details \code{PeaksToGRanges} converts peak information into a
#' \code{\linkS4class{GRanges}} object. Each row in the
#' \code{\linkS4class{GRanges}} object represents a peak with 'CIQ.Up.start' and 'CIQ.Down.end' as start
#' and end coordinates, respectively (see \code{\linkS4class{PSFit}})
#'  Metadata will also include information for the total PETs,
#' the p-value and the FDR of each peak.
#'
#' @param object An object of class \code{\linkS4class{PSFit}}.
#' @param threshold A numeric with the FDR cut-off threshold used to take a
#' subset of significant peaks. If \code{threshold=NULL}
#' then  all the peaks are returned.
#' @param ... Further arguments to be passed to \code{PeaksToGRanges} (not used).
#'
#' @seealso \code{\linkS4class{PSFit}}, \code{\link{PeaksToNarrowPeak}}
#---#'define default:
#' @name PeaksToGRanges
#' @export
#' @include AllClasses.R
#' @importFrom S4Vectors metadata
#' @importFrom GenomicRanges GRanges seqinfo
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlevelsInUse seqlevels
#'
#'
#' @examples
#' #load Self-ligated data: (class=PSFit)
#' load(system.file('extdata', 'MACPET_psfitData.rda', package = 'MACPET'))
#' class(MACPET_psfitData)
#' PeaksToGRanges(object=MACPET_psfitData,threshold=1e-5)
#'
#'
# default:
PeaksToGRanges = function(object, ...) {
    UseMethod("PeaksToGRanges", object = object)
}
#' @rdname PeaksToGRanges
#' @method PeaksToGRanges default
#' @export
PeaksToGRanges.default = function(object, ...) {
    stop("No PeaksToGRanges method for class ", class(object), ".", call. = FALSE)
}
#' @rdname PeaksToGRanges
#' @method PeaksToGRanges PSFit
#' @return For \code{\linkS4class{PSFit}} class, a
#' \code{\linkS4class{GRanges}}
#' object created by the estimated
#' peak information including metadata columns for the total PETs,
#' the p-value and the FDR of each peak.
#' @export
PeaksToGRanges.PSFit = function(object, threshold = NULL, ...) {
    # global variables for Rcheck:
    FDR = NULL
    #--------------------------
    # keep significant peaks:
    Peaks = S4Vectors::metadata(object)$Peaks.Info
    if (is.numeric(threshold)) {
        Peaks = subset(Peaks, FDR < threshold)
        if (nrow(Peaks) == 0) {
            stop("No significant peaks in the data!Try a higher threshold.", call. = FALSE)
        }
    } else {
        warning("No threshold given, all the peaks are returned.")
    }
    # keep seqinfo:
    Seqinformation = GenomicRanges::seqinfo(object)
    # make GRanges:
    RES = GenomicRanges::GRanges(seqnames = Peaks$Chrom, seqinfo = Seqinformation,
        ranges = IRanges::IRanges(start = round(Peaks$CIQ.Up.start), end = round(Peaks$CIQ.Down.end)))
    # reduce unused levels:
    LevelsUsed = GenomeInfoDb::seqlevelsInUse(RES)
    GenomeInfoDb::seqlevels(RES) = LevelsUsed
    # add metadata:
    RES$TotPETs = Peaks$Pets
    RES$p.value = Peaks$p.value
    RES$FDR = Peaks$FDR
    return(RES)
}
