#' @title Convert GInteraction object to PSelf object
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#' @references
#' Vardaxis I, Drabløs F, Rye M and Lindqvist BH (2018). \emph{MACPET: Model-based Analysis for ChIA-PET}.
#' To be published.
#'
#' @description \code{ConvertToPSelf} converts a \code{\link[InteractionSet:GInteractions-class]{GInteractions}} object to
#' class to \code{\linkS4class{PSelf}} object.
#'
#' @details \code{\link{MACPETUlt}} at State 2 separates the Inter-chromosomal,
#'  Intra-chromosomal and Self-ligated PETs by taking
#' the paired-end BAM/SAM file as input. However the user might only have
#' Self-ligated data available and already separated from the Inter/Intra-chromosomal
#' PETs. \code{ConvertToPSelf} can then be used in the Self-ligated data to convert
#' a \code{\link[InteractionSet:GInteractions-class]{GInteractions}} object containing only the Self-ligated
#' PETs to a \code{\linkS4class{PSelf}} class for further analysis in Stage 3.
#' The object will be saved in the \code{S2_AnalysisDir} directory with the
#'  name \code{SA_prefix_pselfData}.
#' Note that if \code{S2_BlackList==TRUE} then the \code{\link[InteractionSet:GInteractions-class]{GInteractions}}
#' object given as input has to include the genome name in the \code{seqinfo} slot.
#' Also, the sequences lengths are mandatory in the \code{seqinfo} slot since they
#' are used in stage 3 of the analysis.
#'
#' @param object An object of \code{\link[InteractionSet:GInteractions-class]{GInteractions}} class.
#' @param ... not used.
#' @param S2_BlackList See  \code{\link{MACPETUlt}}.
#' @param SA_prefix See  \code{\link{MACPETUlt}}.
#' @param S2_AnalysisDir The directory in which the object will be saved.
#'
#' @seealso \code{\linkS4class{PSelf}}
#---#'define default:
#' @name ConvertToPSelf
#' @export
#' @include AllClasses.R
#' @importFrom GenomeInfoDb genome seqlengths seqinfo seqnames
#' @importFrom methods is
#'
#' @examples
#' #load Self-ligated data: (class=PSelf)
#' load(system.file('extdata', 'MACPET_pselfData.rda', package = 'MACPET'))
#' class(MACPET_pselfData)
#'
#' object=MACPET_pselfData
#' #--remove information and convert to GInteractions:
#' S4Vectors::metadata(object)=list(NULL)
#' class(object)='GInteractions'
#' #----input parameters
#' S2_BlackList=TRUE
#' SA_prefix='MACPET'
#' S2_AnalysisDir=file.path(tempdir(),'MACPETtest')
#' if(!dir.exists(S2_AnalysisDir)) dir.create(S2_AnalysisDir)
#'
#' ConvertToPSelf(object=object,
#'                       S2_BlackList=S2_BlackList,
#'                       SA_prefix=SA_prefix,
#'                       S2_AnalysisDir=S2_AnalysisDir)
#' #load object:
#' rm(MACPET_pselfData)#old object
#' load(file.path(S2_AnalysisDir,'MACPET_pselfData'))
#' class(MACPET_pselfData)
#' #-----delete test directory:
#' unlink(S2_AnalysisDir,recursive=TRUE)
#'
# default:
ConvertToPSelf = function(object, ...) {
    UseMethod("ConvertToPSelf", object = object)
}
#' @rdname ConvertToPSelf
#' @method ConvertToPSelf default
#' @export
ConvertToPSelf.default = function(object, ...) {
    stop("No ConvertToPSelf method for class ", class(object), ".", call. = FALSE)
}
#' @rdname ConvertToPSelf
#' @method ConvertToPSelf GInteractions
#' @return An object of class \code{\linkS4class{PSelf}}.
#' @export
ConvertToPSelf.GInteractions = function(object, S2_BlackList, SA_prefix, S2_AnalysisDir,
    ...) {
    #----R-check:
    #----
    #------------
    # check directory:
    #------------
    if (!methods::is(S2_AnalysisDir, "character")) {
        stop("S2_AnalysisDir:", S2_AnalysisDir, " is not a directory!", call. = FALSE)
    } else if (!dir.exists(S2_AnalysisDir)) {
        stop("S2_AnalysisDir:", S2_AnalysisDir, " directory does not exist!", call. = FALSE)
    }
    #------------
    # check prefix:
    #------------
    if (!methods::is(SA_prefix, "character")) {
        stop("SA_prefix: ", SA_prefix, " variable has to be a string!", call. = FALSE)
    } else if (nchar(SA_prefix) == 0) {
        stop("SA_prefix: ", SA_prefix, " variable has to be a non-empty string!",
            call. = FALSE)
    }
    #------------
    # check seqinfo in object:
    #------------
    SeqInfo = GenomeInfoDb::seqinfo(object)
    # names:
    Names = GenomeInfoDb::seqnames(SeqInfo)
    if (any(is.na(Names))) {
        WhichNA = which(is.na(Names))
        stop("GenomeInfoDb::seqnames(object) has NA names at positions: ", paste(WhichNA,
            collapse = "/"), call. = FALSE)
    }
    # sizes:
    Lengths = GenomeInfoDb::seqlengths(SeqInfo)
    if (any(is.na(Lengths))) {
        WhichNA = which(is.na(Lengths))
        stop("GenomeInfoDb::seqlengths(object) has NA lengths at positions: ", paste(names(Lengths[WhichNA]),
            collapse = "/"), call. = FALSE)
    }
    # get genome:
    Genome = GenomeInfoDb::genome(SeqInfo)
    Genome = unique(Genome)
    if (length(Genome) > 1) {
        stop("There are more than one genomes defined in the object. Those are: ",
            paste(Genome, collapse = "/"), call. = FALSE)
    }
    #------------
    # check S2_BlackList:
    #------------
    if (!methods::is(S2_BlackList, "logical") && !methods::is(S2_BlackList, "GRanges")) {
        stop("S2_BlackList: has to be logical or a GRanges object!", call. = FALSE)
    } else if (isTRUE(S2_BlackList) & !Genome %in% names(sysdata)) {
        LogFile = paste("The genome: ", Genome, " is not one of the following: ",
            paste(names(sysdata), collapse = "/"), ". No black listed regions will be removed!")
        warning(LogFile)
    } else if (isTRUE(S2_BlackList) & Genome %in% names(sysdata)) {
        # get black list:
        S2_BlackList = sysdata[[Genome]]
    } else {
        S2_BlackList = NULL
    }
    #-------------------------
    #---------Get correct GInteractions object:
    #-------------------------
    object = GInteractionsCovnert_fun(S2_PairedData = object, S2_BlackList = S2_BlackList,
        PselfConvert = TRUE)
    if (length(object) == 0) {
        stop("The object contained only black-listed regions!", call. = FALSE)
    }
    #-------------------------
    #---------Check if Inter-chromosomal:
    #-------------------------
    Nobject = length(object)
    object = subset(object, !is.na(object$Dist))
    Nreduced = length(object)
    if (Nobject != Nreduced) {
        cat("Total Inter-Chromsomal removed from data:", Nobject - Nreduced, "\n")
    }
    if (Nreduced == 0) {
        # only Inter-chromosomal:
        stop("object contained only Inter-chromosomal PETs!")
    }
    #-------------------------
    #-----Convert to PSelf class:
    #-------------------------
    SelfIndicator = seq_len(Nreduced)  #all the data
    Nself = FindSelf_fun(S2_PairedData = object, SelfIndicator = SelfIndicator, S2_AnalysisDir = S2_AnalysisDir,
        SA_prefix = SA_prefix)
    cat("The PSelf object is saved in: \n", S2_AnalysisDir)
}
