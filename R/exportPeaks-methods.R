#' @title Exports peaks to csv file
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#'
#' @references
#' Vardaxis I, Drabløs F, Rye M and Lindqvist BH (2018). \emph{MACPET: Model-based Analysis for ChIA-PET}.
#' To be published.
#' @description \code{exportPeaks} is an S3 method for the  \linkS4class{PSFit}
#' class. It exports peak information to a csv file in a given directory.
#'
#' @param object An object of \code{\linkS4class{PSFit}}  class.
#' @param savedir A string with the directory to save the output.
#' @param file.out A string with the name of the output to be saved to
#' \code{savedir}.
#' @param threshold A numeric indicating the FDR cut-off used for subseting
#' significant peaks. If NULL all the peaks are returned.
#' @param ... (not used).
#'
#' @seealso \linkS4class{PSFit}
#---#'define default:
#' @name exportPeaks
#' @export
#' @include AllClasses.R
#' @importFrom S4Vectors metadata
#' @importFrom utils write.table
#' @examples
#'
#' #Create a temporary forder, or anywhere you want:
#' savedir=file.path(tempdir(),'MACPETtest')
#' dir.create(savedir)#where you will save the results
#'
#' #load Self-ligated data: (class=PSFit)
#' load(system.file('extdata', 'MACPET_psfitData.rda', package = 'MACPET'))
#' class(MACPET_psfitData)
#' exportPeaks(object=MACPET_psfitData,file.out='Peaks',threshold=1e-5,savedir=savedir)
#'
#' #-----delete test directory:
#' unlink(savedir,recursive=TRUE)
# default:
exportPeaks = function(object, ...) {
    UseMethod("exportPeaks", object = object)
}
#' @rdname exportPeaks
#' @method exportPeaks default
#' @export
exportPeaks.default = function(object, ...) {
    stop("No exportPeaks method for class ", class(object), ".", call. = FALSE)
}
#' @rdname exportPeaks
#' @method exportPeaks PSFit
#' @return For \code{\linkS4class{PSFit}} class: a csv file named after the
#' value of \code{file.out} with all the information about the peaks found by the
#' \code{\link{MACPETUlt}} function, plus comments which
#'  explain the column names.
#' @export
exportPeaks.PSFit = function(object, file.out, savedir, threshold = NULL, ...) {
    if (!dir.exists(savedir)) {
        stop("savedir does not exist!", call. = FALSE)
    }
    x = S4Vectors::metadata(object)
    x = x$Peaks.Info
    if (is.numeric(threshold)) {
        x = subset(x, FDR < threshold)
    } else {
        message("No threshold given, all the peaks are returned.")
    }
    if (nrow(x) == 0) {
        stop("Threshold too low: No significant peaks in data. Try a higher threshold.",
            call. = FALSE)
    }
    x = x[, -which(colnames(x) %in% c("Region", "Peak"))]
    # add comments to CSV file:
    Desc1 = "#Output of peaks found by the algorithm."
    Desc2 = "#Each row represents a peak with the following information:"
    Chrom = "#Chrom: The name of the chromosome"
    Pets = "#Pets: Total PETs in the peak."
    Peak.Summit = "#Peak.Summit: Summit of the peak."
    Up.Summit = "#Up.Summit: Summit of the left-stream PETs."
    Down.Summit = "#Down.Summit: Summit of the right-stream PETs."
    CIQ.Up.start = "#CIQ.Up.start: Start of 95 Quantile confidence interval for the left-stream PETs."
    CIQ.Up.end = "#CIQ.Up.end: End of 95 Quantile confidence interval for the left-stream PETs."
    CIQ.Up.size = "#CIQ.Up.size: Size of 95 Quantile confidence interval for the left-stream PETs."
    CIQ.Down.start = "#CIQ.Down.start: Start of 95 Quantile confidence interval for the right-stream PETs."
    CIQ.Down.end = "#CIQ.Down.end: End of 95 Quantile confidence interval for the right-stream PETs."
    CIQ.Down.size = "#CIQ.Down.size: Size of 95 Quantile confidence interval for the right-stream PETs."
    CIQ.Peak.size = "#CIQ.Peak.size: Size of the Peak based on the interval (CIQ.Up.start,CIQ.Down.end)."
    sdx = "#sdx: The standard deviation of the upstream PETs."
    lambdax = "#lambdax: The skewness of the upstream PETs."
    sdy = "#sdy: The standard deviation of the downstream PETs."
    lambday = "#lambday: The skewness of the downstream PETs."
    lambdaUp = "#lambdaUp: The expected number of PETs in the left-stream Peak region by random chance."
    FoldEnrichUp = "#FoldEnrichUp: Fold enrichment for the left-stream Peak region."
    p.valueUp = "#p.valueUp: p-value for the left-stream Peak region."
    lambdaDown = "#lambdaDown: The expected number of PETs in the right-stream Peak region by random chance."
    FoldEnrichDown = "#FoldEnrichDown: Fold enrichment for the right-stream Peak region."
    p.valueDown = "#p.valueDown: p-value for the right-stream Peak region."
    p.value = "#p.value: p-value for the Peak (p.valueUp*p.valueDown)."
    FDRUp = "#FDRUp: FDR correction for the left-stream Peak region."
    FDRDown = "#FDRDown: FDR correction for the right-stream Peak region."
    FDR = "#FDR:  FDR correction for the Peak."
    skip = "\n"
    writeLines(paste(Desc1, Desc2, skip, Chrom, Pets, Peak.Summit, Up.Summit, Down.Summit,
        CIQ.Up.start, CIQ.Up.end, CIQ.Up.size, CIQ.Down.start, CIQ.Down.end, CIQ.Down.size,
        CIQ.Peak.size, sdx, lambdax, sdy, lambday, lambdaUp, FoldEnrichUp, p.valueUp,
        lambdaDown, FoldEnrichDown, p.valueDown, p.value, FDRUp, FDRDown, FDR,
        skip, sep = "\n"), file.path(savedir, paste(file.out, ".csv", sep = "")))
    suppressWarnings(utils::write.table(x, quote = FALSE, file = file.path(savedir,
        paste(file.out, ".csv", sep = "")), sep = ";", col.names = TRUE, row.names = FALSE,
        qmethod = "double", append = TRUE))
    return("The output is saved at savedir")
}
