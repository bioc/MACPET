#' @title plot methods for MACPET classes
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#'
#' @references
#' Vardaxis I, Drabløs F, Rye M and Lindqvist BH (2018). \emph{MACPET: Model-based Analysis for ChIA-PET}.
#' To be published.
#'
#' Vardaxis I, Drabløs F, Rye M and Lindqvist BH (2018). \emph{MACPET: Complete pipeline for ChIA-PET}.
#' To be published.
#'
#' @description Different plot methods for the classes in the
#' \code{\link{MACPET}} package.
#'
#' @param x An object of correct class used to create different plots.
#' @param ... further arguments to be passed in the plot functions.
#'
#' @seealso \code{\linkS4class{PSelf}},
#' \code{\linkS4class{PSFit}}
#' \code{\linkS4class{PInter}}, \code{\linkS4class{PIntra}},
#' \code{\linkS4class{GenomeMap}}
#'
#' @name plot
#' @include AllClasses.R
#' @importFrom plyr ddply . ldply
#' @importFrom utils methods
#' @importFrom gtools mixedsort
#' @importFrom S4Vectors metadata
NULL
# > NULL
#---PInter:
#' @rdname plot
#' @method plot PInter
#' @return For the \code{\linkS4class{PInter}} class:
#' A network plot.
#'  Each node is a chromosome with size proportional to the total PETs of the
#'  corresponding chromosome. Edges connect chromosomes which have common PETs,
#'   where the thickness of an edge is proportional on the total number of PETs
#'   connecting the two chromosomes.
#' @export
#'
#' @examples
#'
#' #load Inter-chromosomal data:
#' load(system.file('extdata', 'MACPET_pinterData.rda', package = 'MACPET'))
#' class(MACPET_pinterData)
#' requireNamespace('igraph')
#' plot(MACPET_pinterData)
plot.PInter = function(x, ...) {
    # global variables for Rcheck:
    V1 = NULL
    #--------------------------
    #-------check package:
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("igraph needed for this function to work. Please install it.", call. = FALSE)
    }
    PETcounts = S4Vectors::metadata(x)$InteractionCounts
    PETcounts = as.matrix(PETcounts)
    nodes = PETcounts
    nodes = as.data.frame(colSums(nodes))
    nodes$sepnames1 = rownames(nodes)
    colnames(nodes) = c("V1", "from")
    nodes = nodes[, c("from", "V1")]
    nodes$V1 = nodes$V1/max(nodes$V1)
    nodes$from = as.character(nodes$from)
    # network plot:
    name.edges = colnames(PETcounts)
    edges = split(PETcounts, rep(seq_len(ncol(PETcounts)), each = nrow(PETcounts)))
    edges = plyr::ldply(seq_len(length(edges)), function(i, name.edges, edges) {
        data.frame(from = name.edges[i], to = name.edges, V1 = edges[[i]])
    }, name.edges = name.edges, edges = edges)
    edges$V1 = edges$V1/max(edges$V1)
    edges = subset(edges, V1 != 0)
    edges$from = as.character(edges$from)
    edges$to = as.character(edges$to)
    net = igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
    res = igraph::plot.igraph(net, edge.color = "blue", vertex.color = "red", edge.arrow.size = 0,
        vertex.size = igraph::V(net)$V1 * 10, edge.width = igraph::E(net)$V1 * 2,
        main = "Inter Interaction Network Plot")
    return(res)
}
#---Intra PETs:
#' @rdname plot
#' @method plot PIntra
#' @return For the \code{\linkS4class{PIntra}} class:
#' A bar-plot. Each bar
#' represents the total number of Intra-chromosomal PETs for each chromosome in
#' the data.
#' @export
#'
#' @examples
#'
#' #load Intra-chromosomal data:
#' load(system.file('extdata', 'MACPET_pintraData.rda', package = 'MACPET'))
#' class(MACPET_pintraData)
#' requireNamespace('ggplot2')
#' plot(MACPET_pintraData)
plot.PIntra = function(x, ...) {
    # global variables for Rcheck:
    Chrom = NULL
    #--------------------------
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 needed for this function to work. Please install it.", call. = FALSE)
    }
    x = S4Vectors::metadata(x)$InteractionCounts
    x = data.frame(Chrom = rep(x$Chrom, x$Counts))
    res = ggplot2::ggplot(x, ggplot2::aes(x = Chrom, fill = Chrom)) + ggplot2::geom_bar() +
        ggplot2::xlab("Chromosome") + ggplot2::ylab("Intra count") + ggplot2::ggtitle("Intra-chromosomal PET counts by chromosome") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    return(res)
}
#---Self PETs:
#' @rdname plot
#' @method plot PSelf
#' @return For the \code{\linkS4class{PSelf}} class:
#' A bar-plot. Each bar
#' represents the total number of Self-ligated PETs for each chromosome in
#' the data.
#' @export
#'
#' @examples
#'
#' #load Self-ligated data:
#' load(system.file('extdata', 'MACPET_pselfData.rda', package = 'MACPET'))
#' class(MACPET_pselfData)
#' requireNamespace('ggplot2')
#' plot(MACPET_pselfData)
plot.PSelf = function(x, ...) {
    # global variables for Rcheck:
    Chrom = NULL
    #--------------------------
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 needed for this function to work. Please install it.", call. = FALSE)
    }
    x = S4Vectors::metadata(x)$Self_info
    x = data.frame(Chrom = rep(x$Chrom, x$PET.counts))
    res = ggplot2::ggplot(x, ggplot2::aes(x = Chrom, fill = Chrom)) + ggplot2::geom_bar() +
        ggplot2::xlab("Chromosome") + ggplot2::ylab("PET counts") + ggplot2::ggtitle("Self-ligated PET counts by chromosome") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    return(res)
}
#----PSFit:
#' @rdname plot
#' @method plot PSFit
#' @return For the \code{\linkS4class{PSFit}} class:
#'  Different plots depenting on the \code{kind} argument.
#' @param RegIndex an integer indicating which region to plot (1 means the biggest
#' in terms of total PETs.)
#' @param threshold The FDR cut-off when plotting the total significant peaks/interactions for each chromosome in the data.
#' @param kind A string with one of the following arguments. Note that if a region visualization is plotted, the vertical lines represent peak-summits.
#'  \describe{
#'   \item{\code{PETcounts}}{ For a bar-plot of the PET-counts in each chromosome.}
#'   \item{\code{RegionCounts}}{ For a bar-plot for the region counts in each chromosome.}
#'   \item{\code{PeakCounts}}{ For a bar-plot for the Peak-counts in each chromosome.}
#'   \item{\code{RegionPETs}}{ For a ggplot for a visualization of the PETs in a region.}
#'   \item{\code{RegionTags}}{ For a ggplot for a visualization of the tags in a region. The tags are classified by stream (upper/lower)}
#'   \item{\code{PeakPETs}}{ For a ggplot for a visualization of the PETs in a region. The PETs are classified by the peak they belong to.}
#'   \item{\code{PeakTags}}{ For a ggplot for a visualization of the tags in a region. The tags are classified by the peak they belong to.}
#'   \item{\code{SigPETCounts}}{ For a bar-plot with the significant PET-counts in each chromosome.}
#'   \item{\code{SigRegionCounts}}{For a bar-plot with the significant region-counts in each chromosome.}
#'   \item{\code{SigPeakCounts}}{For a bar-plot with the significant peak-counts in each chromosome.}
#'
#'  }
#'
#' @export
#'
#' @examples
#'
#' #load Self-ligated data:
#' load(system.file('extdata', 'MACPET_psfitData.rda', package = 'MACPET'))
#' class(MACPET_psfitData)
#' requireNamespace('ggplot2')
#' plot(MACPET_psfitData,kind='PETcounts')
#' plot(MACPET_psfitData,kind='PeakCounts')
#' plot(MACPET_psfitData,kind='PeakPETs',RegIndex=1)
#' plot(MACPET_psfitData,kind='PeakTags',RegIndex=1)
plot.PSFit = function(x, kind, RegIndex = NULL, threshold = NULL, ...) {
    # global variables for Rcheck:
    FDR = RegCount = X = Y = ymin = ymax = Dist = Tag = Stream = PeakID = NULL
    #--------------------------
    # check that package exists:
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 needed for this function to work if
             create.self_intra.image==T. Please install it.",
            call. = FALSE)
    }
    if (!kind %in% c("PETcounts", "RegionCounts", "PeakCounts", "RegionPETs", "RegionTags",
        "PeakPETs", "PeakTags", "SigPETCounts", "SigRegionCounts", "SigPeakCounts")) {
        stop("kind has been given wrong value!", call. = FALSE)
    }
    if (!methods::is(threshold, "numeric"))
        threshold = NULL
    Peaks.Info = S4Vectors::metadata(x)$Peaks.Info
    if (!is.null(threshold))
        Peaks.Info = subset(Peaks.Info, FDR < threshold)
    if (nrow(Peaks.Info) == 0)
        stop("threshold too low, try a lower one!", call. = FALSE)
    if (kind == "PETcounts") {
        # plot PET count for chromosomes:
        class(x) = "PSelf"
        res = plot.PSelf(x = x)
    } else if (kind == "SigPETCounts") {
        # plot significant pet counts:
        SigPETCounts = plyr::ddply(Peaks.Info, plyr::.(Chrom), function(y) sum(y$Pets))
        SigPETCounts = data.frame(Chrom = rep(SigPETCounts$Chrom, SigPETCounts$V1))
        res = ggplot2::ggplot(SigPETCounts, ggplot2::aes(x = Chrom, fill = Chrom)) +
            ggplot2::geom_bar() + ggplot2::xlab("Chromosome") + ggplot2::ylab("Significant PET counts") +
            ggplot2::ggtitle("significant Self-ligated PET counts by chromosome") +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    } else if (kind == "RegionCounts") {
        # plot region count for chromosomes:
        RegionCounts = S4Vectors::metadata(x)$Self_info
        RegionCounts = data.frame(Chrom = rep(RegionCounts$Chrom, RegionCounts$Region.counts))
        res = ggplot2::ggplot(RegionCounts, ggplot2::aes(x = Chrom, fill = Chrom)) +
            ggplot2::geom_bar() + ggplot2::xlab("Chromosome") + ggplot2::ylab("Region counts") +
            ggplot2::ggtitle("Self-ligated region counts by chromosome") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    } else if (kind == "SigRegionCounts") {
        # plot significant region count for chromosomes:
        SigRegionCounts = plyr::ddply(Peaks.Info, plyr::.(Chrom), function(y) nrow(y))
        SigRegionCounts = data.frame(Chrom = rep(SigRegionCounts$Chrom, SigRegionCounts$V1))
        res = ggplot2::ggplot(SigRegionCounts, ggplot2::aes(x = Chrom, fill = Chrom)) +
            ggplot2::geom_bar() + ggplot2::xlab("Chromosome") + ggplot2::ylab("Region counts") +
            ggplot2::ggtitle("Significant Self-ligated region counts by chromosome") +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    } else if (kind == "PeakCounts") {
        PeakCounts = S4Vectors::metadata(x)$Self_info
        PeakCounts = data.frame(Chrom = rep(PeakCounts$Chrom, PeakCounts$Peak.counts))
        res = ggplot2::ggplot(PeakCounts, ggplot2::aes(x = Chrom, fill = Chrom)) +
            ggplot2::geom_bar() + ggplot2::xlab("Chromosome") + ggplot2::ylab("Peak counts") +
            ggplot2::ggtitle("Self-ligated peak counts by chromosome") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    } else if (kind == "SigPeakCounts") {
        SigPeakCounts = plyr::ddply(Peaks.Info, plyr::.(Chrom), function(y) nrow(y))
        SigPeakCounts = data.frame(Chrom = rep(SigPeakCounts$Chrom, SigPeakCounts$V1))
        res = ggplot2::ggplot(SigPeakCounts, ggplot2::aes(x = Chrom, fill = Chrom)) +
            ggplot2::geom_bar() + ggplot2::xlab("Chromosome") + ggplot2::ylab("Peak counts") +
            ggplot2::ggtitle("Significant Self-ligated peak counts by chromosome") +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    } else if (kind %in% c("RegionPETs", "RegionTags", "PeakPETs", "PeakTags")) {
        # (with summits):
        #----choose the region to plot:
        Classification.Info = S4Vectors::metadata(x)$Classification.Info
        xdf = as.data.frame(x)
        Chrom = xdf$seqnames1
        Chrom = as.character(Chrom)
        Classification.Info$Chrom = Chrom[Classification.Info$MainIndex]
        Classification.Info$RegCount = paste(Classification.Info$Region, "-", Classification.Info$Chrom,
            sep = "")
        MaxReg = table(Classification.Info$RegCount)
        MaxReg = sort(MaxReg, decreasing = TRUE)
        if (is.null(RegIndex) | !is.numeric(RegIndex)) {
            # take max:
            which.id.region = max(which(MaxReg == max(MaxReg)))
            which.id.region = names(MaxReg)[which.id.region]
        } else {
            which.id.region = names(MaxReg)[min(RegIndex, length(MaxReg))]
        }
        # take subset which will be plotted:
        Classification.Info = subset(Classification.Info, RegCount == which.id.region)
        xsub = xdf[Classification.Info$MainIndex, ]
        xsub$Peak.ID = 0
        xsub$Peak.ID = Classification.Info$Peak.ID
        xsub = xsub[, c("start1", "end1", "start2", "end2", "Peak.ID")]
        # sort:
        Tosort = which(xsub$start1 > xsub$start2)
        if (length(Tosort) > 0) {
            Start2new = xsub$start1[Tosort]
            End2new = xsub$end1[Tosort]
            strand2new = xsub$strand1[Tosort]
            xsub$start1[Tosort] = xsub$start2[Tosort]
            xsub$end1[Tosort] = xsub$end2[Tosort]
            xsub$strand1[Tosort] = xsub$strand2[Tosort]
            xsub$start2[Tosort] = Start2new
            xsub$end2[Tosort] = End2new
            xsub$strand2[Tosort] = strand2new
        }
        xsub$Dist = xsub$end2 - xsub$start1 + 1
        # take the peak summits of the region:
        Peaks.Info = S4Vectors::metadata(x)$Peaks.Info
        Peaks.Info$RegCount = paste(Peaks.Info$Region, "-", Peaks.Info$Chrom, sep = "")
        Peaks.Info = subset(Peaks.Info, RegCount == which.id.region)
        if (nrow(Peaks.Info) == 0) {
            Peaks.Info = NULL
        }
        #-----plot according to what is asked:
        if (kind == "RegionPETs") {
            # plot region PETs, no strand info.
            Dens = data.frame(Y = (xsub$start1 + xsub$end2)/2, ymin = xsub$start1,
                ymax = xsub$end2, X = xsub$Dist)
            res = ggplot2::ggplot(Dens, ggplot2::aes(x = X, y = Y, ymin = ymin, ymax = ymax)) +
                ggplot2::geom_errorbar(width = 35) + ggplot2::coord_flip() + ggplot2::ggtitle("Visualization of PETs in region") +
                ggplot2::ylab("Midpoints of PETs") + ggplot2::xlab("PET sizes") +
                ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
            if (!is.null(Peaks.Info)) {
                # add peaks
                res = res + ggplot2::geom_hline(yintercept = Peaks.Info$Peak.Summit,
                  color = "red", linetype = "dashed") + ggplot2::ggtitle("Visualization of PETs in region. Vertical lines correspond to peak-summits.")
            }
        } else if (kind == "RegionTags") {
            # plot Region tags(with summits):
            Dens = data.frame(Tag = c((xsub$start1 + xsub$end1)/2, (xsub$start2 +
                xsub$end2)/2), ymin = c(xsub$start1, xsub$start2), ymax = c(xsub$end1,
                xsub$end2), Dist = c(xsub$Dist, xsub$Dist), Stream = c(rep("Upper",
                nrow(xsub)), rep("Lower", nrow(xsub))))
            res = ggplot2::ggplot(Dens, ggplot2::aes(x = Dist, y = Tag, ymin = ymin,
                ymax = ymax, color = factor(Stream))) + ggplot2::geom_errorbar(width = 10) +
                ggplot2::coord_flip() + ggplot2::ggtitle("Visualization of Tags in region") +
                ggplot2::ylab("Midpoints of Tags") + ggplot2::xlab("PET sizes") +
                ggplot2::labs(color = "Stream") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
            if (!is.null(Peaks.Info)) {
                res = res + ggplot2::geom_hline(yintercept = Peaks.Info$Peak.Summit,
                  color = "red", linetype = "dashed") + ggplot2::ggtitle("Visualization of Tags in region. Vertical lines correspond to peak-summits.")
            }
        } else if (kind == "PeakPETs") {
            # (with summits) plot region PETs, no strand info.
            Dens = data.frame(Y = (xsub$start1 + xsub$end2)/2, ymin = xsub$start1,
                ymax = xsub$end2, X = xsub$Dist, PeakID = xsub$Peak.ID)
            res = ggplot2::ggplot(Dens, ggplot2::aes(x = X, y = Y, ymin = ymin, ymax = ymax,
                color = factor(PeakID))) + ggplot2::geom_errorbar(width = 35) + ggplot2::coord_flip() +
                ggplot2::ggtitle("Visualization of PETs in region") + ggplot2::ylab("Midpoints of PETs") +
                ggplot2::xlab("PET sizes") + ggplot2::labs(color = "PeakID") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
            if (!is.null(Peaks.Info)) {
                res = res + ggplot2::geom_hline(yintercept = Peaks.Info$Peak.Summit,
                  color = "red", linetype = "dashed") + ggplot2::ggtitle("Visualization of PETs in region. Vertical lines correspond to peak-summits.")
            }
        } else if (kind == "PeakTags") {
            # (with summits)
            Dens = data.frame(Tag = c((xsub$start1 + xsub$end1)/2, (xsub$start2 +
                xsub$end2)/2), ymin = c(xsub$start1, xsub$start2), ymax = c(xsub$end1,
                xsub$end2), Dist = c(xsub$Dist, xsub$Dist), PeakID = c(xsub$Peak.ID,
                xsub$Peak.ID))
            res = ggplot2::ggplot(Dens, ggplot2::aes(x = Dist, y = Tag, ymin = ymin,
                ymax = ymax, color = factor(PeakID))) + ggplot2::geom_errorbar(width = 10) +
                ggplot2::coord_flip() + ggplot2::ggtitle("Visualization of Tags in region") +
                ggplot2::ylab("Midpoints of Tags") + ggplot2::xlab("PET sizes") +
                ggplot2::labs(color = "PeakID") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
            if (!is.null(Peaks.Info)) {
                res = res + ggplot2::geom_hline(yintercept = Peaks.Info$Peak.Summit,
                  color = "red", linetype = "dashed") + ggplot2::ggtitle("Visualization of Tags in region. Vertical lines correspond to peak-summits.")
            }
        }
    }
    return(res)
}
#---GenomeMap:
#' @rdname plot
#' @method plot GenomeMap
#' @return For the \code{\linkS4class{GenomeMap}} class:
#' A heatmap plot or different kinds of network plots showing the relations between
#' different chromosomes based on the interactions between them.
#' @param Type A string with one of the following values:
#'  \code{'heatmap'}, \code{'network-random'},\code{'network-circle'},\code{'network-sphere'}.
#'  Each one shows the relations between the chromosomes based on their interactions.
#' @export
#'
#' @examples
#'
#' #load Interactions data:
#' load(system.file('extdata', 'MACPET_GenomeMapData.rda', package = 'MACPET'))
#' class(MACPET_GenomeMapData)
#' requireNamespace('igraph')
#'  requireNamespace('reshape2')
#' plot(MACPET_GenomeMapData,Type='network-circle')
plot.GenomeMap = function(x, Type, threshold = NULL, ...) {
    # global variables for Rcheck:
    seqnames1 = NULL
    seqnames2 = NULL
    value = NULL
    V1 = NULL
    #--------------------------
    #-------check package:
    if (!methods::is(Type, "character")) {
        stop("Type variable has to be a character", call. = FALSE)
    } else if (!Type %in% c("heatmap", "network-random", "network-circle", "network-sphere")) {
        stop("Type variable has to be: 'heatmap', 'network-random', 'network-circle' or 'network-sphere'.",
            call. = FALSE)
    }
    # take the data:
    object = GetSignInteractions(object = x, threshold = threshold, ReturnedAs = "GInteractions")
    #--------------------------
    # take cases:
    #--------------------------
    # make the CountsMat:
    object = as.data.frame(object)[, c("seqnames1", "seqnames2")]
    CountsMat = table(object)
    if (Type == "heatmap") {
        # check package:
        if (!requireNamespace("ggplot2", quietly = TRUE) | !requireNamespace("reshape2",
            quietly = TRUE)) {
            stop("ggplot2 and reshape2 are needed for this
                 function to work if heatmap==TRUE.
                 Please install it.",
                call. = FALSE)
        }
        plot.data = reshape2::melt(as.matrix(CountsMat))
        plot.data$value = plot.data$value/max(plot.data$value)
        RES = ggplot2::ggplot(data = plot.data, ggplot2::aes(x = seqnames1, y = seqnames2,
            fill = value)) + ggplot2::geom_tile(color = "red") + ggplot2::ggtitle("Heat-Map plot for Interactions") +
            ggplot2::xlab("Chromosome") + ggplot2::ylab("Chromosome") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        return(RES)
    } else if (Type %in% c("network-random", "network-circle", "network-sphere")) {
        # check package
        if (!requireNamespace("igraph", quietly = TRUE)) {
            stop("igraph needed for this function to work. Please install it.", call. = FALSE)
        }
        # make data:
        nodes = as.data.frame(colSums(CountsMat))
        nodes$sepnames1 = rownames(nodes)
        colnames(nodes) = c("V1", "from")
        nodes = nodes[, c("from", "V1")]
        nodes$V1 = nodes$V1/max(nodes$V1)
        nodes$from = as.character(nodes$from)
        # network plot:
        name.edges = colnames(CountsMat)
        edges = split(CountsMat, rep(seq_len(ncol(CountsMat)), each = nrow(CountsMat)))
        edges = plyr::ldply(seq_len(length(edges)), function(i, name.edges, edges) {
            data.frame(from = name.edges[i], to = name.edges, V1 = edges[[i]])
        }, name.edges = name.edges, edges = edges)
        edges$V1 = edges$V1/max(edges$V1)
        edges = subset(edges, V1 != 0)
        edges$from = as.character(edges$from)
        edges$to = as.character(edges$to)
        net = igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
        if (Type == "network-random") {
            igraph::plot.igraph(net, edge.color = "blue", vertex.color = "red", main = "Interactions Network Plot")
        } else if (Type == "network-circle") {
            igraph::plot.igraph(net, layout = igraph::layout_in_circle, edge.curved = 0.3,
                edge.color = "blue", vertex.color = "red")
        } else if (Type == "network-sphere") {
            igraph::plot.igraph(net, layout = igraph::layout_on_sphere(net), edge.color = "blue",
                vertex.color = "red")
        }
    }
}
