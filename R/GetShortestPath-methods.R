#' @title Find shortest path between pairs of peaks given a set of significant interactions.
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#' @references
#' Vardaxis I, Drabløs F, Rye M and Lindqvist BH (2018). \emph{MACPET: Complete pipeline for ChIA-PET}.
#' To be published.
#'
#' @description \code{GetShortestPath} function subsets the significant interactions of a \code{\linkS4class{GenomeMap}} object
#' given a user specified FDR cut-off. Based on the leftover interactions, it creates a network and finds the length of the
#' shortest path between two user-specified peaks. Currently it only finds the shortest paths between intra-chromosomal peaks.
#' Therefore, the peaks have to be on the same chromosome.
#'
#'
#'
#' @param object An object of class \code{\linkS4class{GenomeMap}}.
#' @param threshold A numeric with the FDR cut-off threshold used to take a
#' subset of significant interactions If \code{threshold=NULL}
#' then  all the interactions will be used
#' @param ChrFrom Character specifying the chromosome of the 'From' peak, for example 'chr1'.
#' @param ChrTo Character specifying the chromosome of the 'To' peak.
#' @param SummitFrom Numeric specifying the peak summit of the 'From' peak.
#' @param SummitTo Numeric specifying the peak summit of the 'To' peak.
#' @param ... Further arguments to be passed to \code{GetShortestPath} (not used).
#'
#' @seealso \code{\linkS4class{GenomeMap}}
#---#'define default:
#' @name GetShortestPath
#' @export
#' @include AllClasses.R
#' @importFrom methods as is
#'
#'
#'
#' @examples
#'
#' #load Interaction data: (class=GenomeMap)
#' load(system.file('extdata', 'MACPET_GenomeMapData.rda', package = 'MACPET'))
#' class(MACPET_GenomeMapData)
#' GetShortestPath(object=MACPET_GenomeMapData,
#'                     threshold = NULL,
#'                     ChrFrom='chr1',
#'                     ChrTo='chr1',
#'                     SummitFrom=10000,
#'                     SummitTo=1000000)
#'
# default:
GetShortestPath = function(object, ...) {
    UseMethod("GetShortestPath", object = object)
}
#' @rdname GetShortestPath
#' @method GetShortestPath default
#' @export
GetShortestPath.default = function(object, ...) {
    stop("No GetShortestPath method for class ", class(object), ".", call. = FALSE)
}
#' @rdname GetShortestPath
#' @method GetShortestPath GenomeMap
#' @return A two-element list with the first element named \code{LinearPathLength}
#' for the linear length of the path between \code{SummitFrom} and \code{SummitTo},
#' and the second element named \code{ThreeDPathLength} for the 3D length of the shortest path
#' between \code{SummitFrom} and \code{SummitTo}.
#' @export
GetShortestPath.GenomeMap = function(object, threshold = NULL, ChrFrom, ChrTo, SummitFrom, 
    SummitTo, ...) {
    # R check:
    seqnames1 = NULL
    seqnames2 = NULL
    #--------------------------
    # Check input here:
    #--------------------------
    if (!methods::is(ChrFrom, "character")) {
        stop("The ChrFrom parameter has to be a character!", call. = FALSE)
    }
    if (!methods::is(ChrTo, "character")) {
        stop("The ChrTo parameter has to be a character!", call. = FALSE)
    }
    if (ChrFrom != ChrTo) {
        stop("ChrFrom: ", ChrFrom, " has to be the same as ChrTo: ", ChrTo, call. = FALSE)
    }
    if (!methods::is(SummitFrom, "numeric")) {
        stop("The SummitFrom parameter has to be numeric!", call. = FALSE)
    } else if (SummitFrom < 0) {
        stop("The SummitFrom parameter has to be positive!", call. = FALSE)
    }
    if (!methods::is(SummitTo, "numeric")) {
        stop("The SummitTo parameter has to be numeric!", call. = FALSE)
    } else if (SummitTo < 0) {
        stop("The SummitTo parameter has to be positive!", call. = FALSE)
    }
    #--------------------------
    # Subset the significant interactions:
    #--------------------------
    SignInteractions = GetSignInteractions(object = object, threshold = threshold, 
        ReturnedAs = "GInteractions")
    # make data frame, keep what you need:
    SignInteractions = as.data.frame(SignInteractions)
    SignInteractions = subset(SignInteractions, seqnames1 == ChrFrom & seqnames2 == 
        ChrTo)
    Ninteractions = nrow(SignInteractions)
    if (Ninteractions == 0) {
        stop("No significant interactions were found between chromosomes: ", ChrFrom, 
            " and ", ChrTo, " at the given threshold!", call. = FALSE)
    }
    SignInteractions = SignInteractions[, c("Anchor1Summit", "Anchor2Summit")]
    #--------------------------
    # Make unique peaks data:
    #--------------------------
    Peaks = unique(c(SignInteractions$Anchor1Summit, SignInteractions$Anchor2Summit, 
        SummitFrom, SummitTo))
    Peaks = sort(Peaks, decreasing = FALSE)
    NPeaksInvolved = length(Peaks)
    Peaks = data.frame(Summit = Peaks, PeakSummitIDS = seq_len(NPeaksInvolved))
    #--------------------------
    # Change the names of the summits to those in peaks ids:
    #--------------------------
    LinearPathLength = abs(SummitFrom - SummitTo)  #the linear path length
    SummitFrom = Peaks$PeakSummitIDS[which(Peaks$Summit == SummitFrom)]
    SummitTo = Peaks$PeakSummitIDS[which(Peaks$Summit == SummitTo)]
    Match1 = match(SignInteractions$Anchor1Summit, Peaks$Summit)
    SignInteractions$Anchor1Summit = Peaks$PeakSummitIDS[Match1]
    Match2 = match(SignInteractions$Anchor2Summit, Peaks$Summit)
    SignInteractions$Anchor2Summit = Peaks$PeakSummitIDS[Match2]
    #--------------------------
    # Create the network structure:
    #--------------------------
    # The edges are in R index
    Network = Initiate_GenomeMap_fun_Rcpp(NPeaksInvolved, Peaks$PeakSummitIDS, Peaks$PeakSummitIDS, 
        Peaks$Summit, 0, FALSE)
    #--------------------------
    # Merge the interactions:
    #--------------------------
    for (i in seq_len(Ninteractions)) {
        #--------------------------
        # take interaction peaks:
        #--------------------------
        Node_kh = sort(c(SignInteractions[i, c(1)], SignInteractions[i, c(2)]), decreasing = FALSE)
        if (Node_kh[1] == Node_kh[2]) 
            next  #everything is done so go to next
        #--------------------------
        # update the leftover interactions:
        #--------------------------
        k = Node_kh[1]
        h = Node_kh[2]
        Math1h = which(SignInteractions$Anchor1Summit == h)
        if (length(Math1h) != 0) 
            SignInteractions$Anchor1Summit[Math1h] = k
        Math2h = which(SignInteractions$Anchor2Summit == h)
        if (length(Math2h) != 0) 
            SignInteractions$Anchor2Summit[Math2h] = k
        #--------------------------
        # update the network:
        #--------------------------
        Network = Update_Network_kh_fun(Network = Network, k = k, h = h)
    }
    #--------------------------
    # Find new shortest path:
    #--------------------------
    if (SummitFrom > SummitTo) {
        Change = SummitFrom
        SummitFrom = SummitTo
        SummitTo = Change
    }
    SPDistances = Dijkstra_GSP_fun_Rcpp(SummitFrom, Network, NPeaksInvolved)
    ThreeDPathLength = SPDistances[SummitTo]
    #--------------------------
    # return
    #--------------------------
    return(list(LinearPathLength = LinearPathLength, ThreeDPathLength = ThreeDPathLength))
}  #done
