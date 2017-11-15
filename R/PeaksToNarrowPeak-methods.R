#' @title Convert Peaks to narrowPeak (BED) object.
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#' @references
#' Vardaxis I, Drabløs F, Rye M and Lindqvist BH (2018). \emph{Model-based Analysis for ChIA-PET (MACPET)}.
#' To be published.
#'
#' @description \code{PeaksToNarrowPeak} converts peaks of an object of \code{\linkS4class{PSFit}}
#' class to \code{narrowPeak} object. The object is saved in a user specified directory and can be used in the
#' MANGO or MICC algorithms for interaction analysis.
#' @details Each Peak in the narrowPeak object is represented by an interval starting from the 'CIQ.Up.start'
#' estimated variable to its 'CIQ.Down.end' (see \code{\linkS4class{PSFit}}). Close Peaks in genomic distance are NOT
#' merged by the \code{PeaksToNarrowPeak} function. However the user can specify a distance window for merging in the
#' MANGO or MICC algorithms. Note also that MANGO and MICC find a self-ligated cut-off by itself which is usually very different than
#' that found by MACPET. We suggest that the user overwrites MANGOS's or MICC's cut-off with that of MACPET.
#'
#'
#' @param object An object of class \code{\linkS4class{PSFit}}.
#' @param threshold A numeric with the FDR cut-off threshold used to take a
#' subset of significant peaks. If \code{threshold=NULL}
#' then  all the peaks are returned.
#' @param ... Further arguments to be passed to \code{PeaksToNarrowPeak} (not used).
#' @param savedir A string with the directory to save the ouput file.
#' @param file.out A string with the name of the output to be saved to
#' \code{savedir}.
#'
#' @seealso \code{\linkS4class{PSFit}}
#---#'define default:
#' @name PeaksToNarrowPeak
#' @export
#' @include AllClasses.R
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom S4Vectors sort
#' @importFrom rtracklayer export
#'
#'
#' @examples
#' #Create a temporary forder, or anywhere you want:
#' savedir=file.path(tempdir(),'MACPETtest')
#' dir.create(savedir)#where you will save the results
#' file.out='MACPET_peaks.narrowPeak'
#'
#' #load Self-ligated data: (class=PSFit)
#' load(system.file('extdata', 'psfitData.rda', package = 'MACPET'))
#' class(psfitData)
#' PeaksToNarrowPeak(object=psfitData,threshold=1e-5,file.out=file.out,savedir=savedir)
#'
#' #-----delete test directory:
#' unlink(savedir,recursive=TRUE)
# default:
PeaksToNarrowPeak = function(object, ...) {
    UseMethod("PeaksToNarrowPeak", object = object)
}
#' @rdname PeaksToNarrowPeak
#' @method PeaksToNarrowPeak default
#' @export
PeaksToNarrowPeak.default = function(object, ...) {
    stop("No PeaksToNarrowPeak method for class ", class(object), ".", call. = FALSE)
}
#' @rdname PeaksToNarrowPeak
#' @method PeaksToNarrowPeak PSFit
#' @return A \code{narrowPeak} object named after the value of \code{file.out} and saved in the \code{savedir.}
#' @export
PeaksToNarrowPeak.PSFit = function(object, threshold = NULL, savedir, file.out, ...) {
    # global variables for Rcheck:
    FDR = NULL
    #--------------------------
    if (any(!is.character(c(savedir, file.out)))) {
        stop("savedir or file.out are not of character class!", call. = FALSE)
    }
    if (!dir.exists(savedir)) {
        stop("savedir does not exist!", call. = FALSE)
    }
    # ----convert to GRanges:
    object = MACPET::PeaksToGRanges(object = object, threshold = threshold)
    # remove metadata:
    object$TotPETs = NULL
    object$p.value = NULL
    object$FDR = NULL
    #----sort by chrom and ranges:
    object = GenomeInfoDb::sortSeqlevels(object)
    object = S4Vectors::sort(object)
    #-----add names:
    object$name = paste("peak_", 1:length(object), sep = "")
    # save:
    connection = file.path(savedir, file.out)
    rtracklayer::export(object, con = connection, format = "BED")
    return("Done! Check savedir!")
}
