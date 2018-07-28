#' @title Get the significant interactions of a \code{\linkS4class{GenomeMap}} object
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#' @references
#' Vardaxis I, Drabløs F, Rye M and Lindqvist BH (2018). \emph{MACPET: Complete pipeline for ChIA-PET}.
#' To be published.
#'
#' @description \code{GetSignInteractions} subsets the significant interactions of a \code{\linkS4class{GenomeMap}} object
#' given a user specified FDR cut-off. It returns/saves different kind of objects based on the \code{ReturnedAs} value.
#'
#' @details  \code{\link{MACPETUlt}} at Stage 4 will run the interaction analysis for the ChIA-PET data.
#' This stage will proccess all the interactions in the data. After Stage 4 is run, the user can use the
#' \code{GetSignInteractions} function on the saved object of class  \code{\linkS4class{GenomeMap}},
#'  to subset the significant interactions given a FDR cut-off.
#'
#'
#' @param object An object of class \code{\linkS4class{GenomeMap}}.
#' @param threshold A numeric with the FDR cut-off threshold used to take a
#' subset of significant interactions If \code{threshold=NULL}
#' then  all the interactions are returned.
#' @param ReturnedAs What type of object to return. Possible values are:
#' \describe{
#' \item{\code{'GInteractions'}}{ Will return a \code{\link[InteractionSet:GInteractions-class]{GInteractions}} object.
#' Each line corresponds to an interaction between two anchors (peaks). It provides information about the intervals of the
#' peaks, the anchor/peak summits, the interaction p-value/FDR and order (see \code{\linkS4class{GenomeMap}} for further details).}
#' \item{\code{'csv'}}{ Will save a .csv file in the \code{SaveTo}. The file will
#' be named after the \code{Prefix} parameter \code{Prefix_InteractionResults.csv}}
#' }
#' @param SaveTo The directory to save the csv file in case of \code{ReturnedAs=='csv'}.
#' @param Prefix A string which will be used a prefix on the object saved in \code{SaveTo}.
#' @param ... Further arguments to be passed to \code{GetSignInteractions} (not used).
#'
#' @seealso \code{\linkS4class{GenomeMap}}
#---#'define default:
#' @name GetSignInteractions
#' @export
#' @include AllClasses.R
#' @importFrom GenomicRanges seqinfo
#' @importFrom S4Vectors metadata
#' @importFrom methods as is
#' @importFrom InteractionSet GInteractions
#' @importFrom utils write.table
#'
#'
#'
#' @examples
#'
#' #load Interaction data: (class=GenomeMap)
#' load(system.file('extdata', 'MACPET_GenomeMapData.rda', package = 'MACPET'))
#' class(MACPET_GenomeMapData)
#' GetSignInteractions(object=MACPET_GenomeMapData,
#'                     threshold = NULL,
#'                     ReturnedAs='GInteractions')
#'
# default:
GetSignInteractions = function(object, ...) {
    UseMethod("GetSignInteractions", object = object)
}
#' @rdname GetSignInteractions
#' @method GetSignInteractions default
#' @export
GetSignInteractions.default = function(object, ...) {
    stop("No GetSignInteractions method for class ", class(object), ".", call. = FALSE)
}
#' @rdname GetSignInteractions
#' @method GetSignInteractions GenomeMap
#' @return For \code{\linkS4class{GenomeMap}} class, different outputs based on the \code{ReturnedAs} parameter.
#' @export
GetSignInteractions.GenomeMap = function(object, threshold = NULL, ReturnedAs = "GInteractions",
    SaveTo = "", Prefix = "MACPET", ...) {
    #--------------------------
    # Check input here:
    #--------------------------
    if ((!is.null(threshold)) && (!methods::is(threshold, "numeric"))) {
        stop("The threshold parameter has to be either NULL or numeric!", call. = FALSE)
    }
    if (!ReturnedAs %in% c("GInteractions", "csv")) {
        stop("The ReturnedAs parameter has to be a string with values: 'GInteractions' or 'csv' !",
            call. = FALSE)
    } else if ((ReturnedAs == "csv") && (!methods::is(SaveTo, "character"))) {
        stop("SaveTo: ", SaveTo, " is not a directory!", call. = FALSE)
    } else if ((ReturnedAs == "csv") && (!dir.exists(SaveTo))) {
        stop("SaveTo: ", SaveTo, " directory does not exist!", call. = FALSE)
    } else if ((ReturnedAs == "csv") && (!methods::is(Prefix, "character"))) {
        stop("Prefix: ", Prefix, " should be string!", call. = FALSE)
    }
    #--------------------------
    # break the object and subset to threshold
    #--------------------------
    # get mapping data/Those are shorted by the order:
    InteractionInfo = S4Vectors::metadata(object)$InteractionInfo
    # subset by threshold:
    if (!is.null(threshold)) {
        # count the line
        CutLine = 0
        for (i in seq_len(nrow(InteractionInfo))) {
            if (InteractionInfo$FDR[i] >= threshold) {
                break
            }
            CutLine = CutLine + 1
        }
        # subset based CutLine
        if (CutLine == 0) {
            stop("The FDR threshold ", threshold, " is too low! Use a higher one!",
                call. = FALSE)
        } else if (CutLine == nrow(InteractionInfo)) {
            message("Threshold ", threshold, " is reached, all the interactions are returned.")
        }
        # subset
        CutSeq = seq_len(CutLine)
        InteractionInfo = InteractionInfo[CutSeq, ]
        object = object[CutSeq]
    } else {
        message("No threshold given, all the interactions are returned.")
    }
    #--------------------------
    # Finalize object:
    #--------------------------
    object$pvalue = InteractionInfo$pvalue
    object$FDR = InteractionInfo$FDR
    object$Order = InteractionInfo$Order
    object$TotalInterPETs = InteractionInfo$TotalInterPETs
    S4Vectors::metadata(object)$InteractionInfo = NULL
    #--------------------------
    # take cases of return:
    #--------------------------
    if (ReturnedAs == "GInteractions") {
        # return object
        object=methods::as(object,"GInteractions")
        return(object)
    } else {
        # create csv and save:
        object = as.data.frame(object)
        object = object[, c(1, 2, 3, 4, 11, 6, 7, 8, 9, 12,
            13, 14, 15, 16)]
        colnames(object) = c("ChromFrom", "StartFrom", "EndFrom", "SizeFrom",
            "SummitFrom", "ChromTo", "StartTo", "EndTo", "SizeTo", "SummitTo", "pvalue",
            "FDR", "Order", "TotalInterPETs")
        # add comments to CSV file:
        Desc1 = "#Output of the significant interactions found by MACPET."
        Desc2 = "#Each row represents an interaction between two peaks with the following information:"
        ChromFrom = "#ChromFrom: The name of the chromosome of the one peak."
        StartFrom = "#StartFrom: The leftmost coordinate of the one peak (CIQ.Up.start)."
        EndFrom = "#EndFrom: The rightmost coordinate of the one peak (CIQ.Down.end)."
        SizeFrom = "#SizeFrom: The size of the one peak (CIQ.Down.end-CIQ.Up.start)."
        SummitFrom = "#SummitFrom: The summit of the one peak."
        ChromTo = "#ChromTo: The name of the chromosome of the other peak."
        StartTo = "#StartTo: The leftmost coordinate of the other peak (CIQ.Up.start)."
        EndTo = "#EndTo: The rightmost coordinate of the other peak (CIQ.Down.end)."
        SizeTo = "#SizeTo: The size of the other peak (CIQ.Down.end-CIQ.Up.start)."
        SummitTo = "#SummitTo: The summit of the other peak."
        pvalue = "#pvalue: The p-value of the interaction."
        FDR = "#FDR: The FDR of the interaction."
        Order = "#Order: The Order of the interaction for which it was added to the model."
        TotalInterPETs = "#TotalInterPETs: The total PETs connecting the two peaks."
        skip = "\n"
        writeLines(paste(Desc1, Desc2, skip, ChromFrom, StartFrom, EndFrom, SizeFrom,
            SummitFrom, ChromTo, StartTo, EndTo, SizeTo, SummitTo, pvalue, FDR, Order,
            TotalInterPETs, skip, sep = "\n"), file.path(SaveTo, paste(Prefix, "_InteractionResults.csv",
            sep = "")))
        suppressWarnings(utils::write.table(object, quote = FALSE, file = file.path(SaveTo,
            paste(Prefix, "_InteractionResults.csv", sep = "")), sep = ";", col.names = TRUE,
            row.names = FALSE, qmethod = "double", append = TRUE))
        # return:
        cat("The output is saved at:\n", SaveTo)
    }
}  #done
