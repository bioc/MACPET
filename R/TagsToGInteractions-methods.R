#' @title Convert PETs to GInteractions object
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#' @references
#' Vardaxis I, Drabløs F, Rye M and Lindqvist BH (2018). \emph{MACPET: Model-based Analysis for ChIA-PET}.
#' To be published.
#'
#' @description \code{TagsToGInteractions} converts the PETs of an object of \code{\linkS4class{PSFit}}
#' class to \code{\link[InteractionSet:GInteractions-class]{GInteractions}} object.
#'
#'
#' @param object An object of class \code{\linkS4class{PSFit}}.
#' @param threshold A numeric for the FDR threshold used to take a subset of
#' significant peaks/binding sites. If \code{threshold=NULL} then
#' all the peaks are returned.
#' @param ... (not used).
#'
#' @seealso \code{\linkS4class{PSFit}}
#---#'define default:
#' @name TagsToGInteractions
#' @export
#' @include AllClasses.R
#' @importFrom S4Vectors metadata
#' @importFrom GenomicRanges seqinfo
#'
#'
#' @examples
#' #load Self-ligated data: (class=PSFit)
#' load(system.file('extdata', 'MACPET_psfitData.rda', package = 'MACPET'))
#' class(MACPET_psfitData)
#' object=TagsToGInteractions(object=MACPET_psfitData,threshold=1e-5)
#' object
#' S4Vectors::metadata(object)$Peaks.Info #peak/binding site information
#'
# default:
TagsToGInteractions = function(object, ...) {
    UseMethod("TagsToGInteractions", object = object)
}
#' @rdname TagsToGInteractions
#' @method TagsToGInteractions default
#' @export
TagsToGInteractions.default = function(object, ...) {
    stop("No TagsToGInteractions method for class ", class(object), ".",
         call. = FALSE)
}
#' @rdname TagsToGInteractions
#' @method TagsToGInteractions PSFit
#' @return For \code{\linkS4class{PSFit}} class:
#' A \code{\link[InteractionSet:GInteractions-class]{GInteractions}} object containing PETs from all
#'  the peaks found in the data (removing noisy and insignificant PETs).
#'  Furthermore, it also includes information about
#'  the binding sites which can be accessed via the
#'  \code{\link[S4Vectors:Annotated-class]{metadata}} function.
#' @export
TagsToGInteractions.PSFit = function(object, threshold = NULL, ...) {
    # global variables for Rcheck:
    FDR = NULL
    #--------------------------
    # keep significant peaks:
    Peaks.info = S4Vectors::metadata(object)$Peaks.Info
    if (nrow(Peaks.info) == 0) {
        stop("No peaks in the data!", call. = FALSE)
    }
    if (is.numeric(threshold)) {
        Peaks.info = subset(Peaks.info, FDR < threshold)
        if (nrow(Peaks.info) == 0) {
            stop("No significant peaks in the data!Try a higher threshold.",
                 call. = FALSE)
        }
    } else {
        warning("No threshold given, all the peaks are returned.")
    }
    Peaks.infoSave = Peaks.info
    Peaks.info = Peaks.info[, c("Chrom", "Region", "Peak")]
    # make data to data.frame to find the tags. keep seqinfo:
    Seqinformation = GenomicRanges::seqinfo(object)
    # make df of the data:
    objectdf = data.frame(object)
    # take chromosomes:
    ChromInfo = objectdf$seqnames1
    ChromInfo = as.character(ChromInfo)
    # Make stings and find overlaps:
    From.significant = paste(Peaks.info$Chrom, "-", Peaks.info$Region, "-",
                             Peaks.info$Peak, sep = "")
    # take EMinformation:
    Classification.Info = S4Vectors::metadata(object)$Classification.Info
    To.data = paste(ChromInfo[Classification.Info$MainIndex], "-",
                    Classification.Info$Region, "-", Classification.Info$Peak.ID,
                    sep = "")
    Keep.data = which(To.data %in% From.significant)
    # take MainIndex:
    MainIndex = Classification.Info$MainIndex[Keep.data]
    # keep those in object:
    object = object[MainIndex]  #keep those in peaks only
    S4Vectors::metadata(object)$Peaks.Info = Peaks.infoSave
    class(object) = "GInteractions"
    return(object)
}
#defining the method as S4 method:
#' @rdname TagsToGInteractions
#' @aliases TagsToGInteractions,PSFit,TagsToGInteractions-method
#' @export
setMethod("TagsToGInteractions", "PSFit", TagsToGInteractions.PSFit)
