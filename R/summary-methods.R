#' @title summary methods for the MACPET classes.
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#'
#' @references
#' Vardaxis I, Drabløs F, Rye M and Lindqvist BH (2018). \emph{MACPET: Model-based Analysis for ChIA-PET}.
#' To be published.
#'
#' @description Different summary methods for the classes in the
#' \code{\link{MACPET}} package.
#'
#' @param object An object of correct class used to create different summaries.
#' @param threshold A numeric representing the FDR cut-off or summarizing
#' singificant peaks (or interactions), if NULL the summary is based on all the peaks found.
#' @param ... Further arguments to be passed to the summary function.
#' @param heatmap TRUE or FALSE indicating whether the user wants to plot a
#' heat-map plot for the Intra/Inter-chromosomal PET counts within chromosomes or
#'  between different chromosomes.
#'
#' @return A summary of the \code{object} and a heat-map plot depending on the
#' class of the input.
#'
#' @seealso \code{\linkS4class{PSelf}},
#' \code{\linkS4class{PSFit}},
#' \code{\linkS4class{PInter}},\code{\linkS4class{PIntra}},
#' \code{\linkS4class{GenomeMap}}
#'
#' @name summary
#' @include AllClasses.R
#'
#' @importFrom utils methods
#' @importFrom S4Vectors metadata
#' @importFrom GenomeInfoDb genome
#' @importFrom knitr kable
NULL
# > NULL
#' @rdname summary
#' @export
#' @method summary PSelf
#'
#' @examples
#'
#' #load Self-ligated data: (class=PSelf)
#' load(system.file('extdata', 'MACPET_pselfData.rda', package = 'MACPET'))
#' class(MACPET_pselfData)
#' summary(MACPET_pselfData)
summary.PSelf = function(object, ...) {
    # prepare data:
    N = length(object)
    SLmean = S4Vectors::metadata(object)$SLmean
    hg = unique(GenomeInfoDb::genome(GenomeInfoDb::seqinfo(object)))
    MaxSize = S4Vectors::metadata(object)$MaxSize
    MinSize = S4Vectors::metadata(object)$MinSize
    PET.counts = S4Vectors::metadata(object)$Self_info
    colnames(PET.counts) = c("Chrom", "Self-lig.")
    cat("|-Self-ligatead PETs|\n")
    cat("|------Summary------|")
    Info = data.frame(Total.PETs = N, SLmean = SLmean, hg = hg, MIN = paste(MinSize,
        "bp"), MAX = paste(MaxSize, "bp"), class = class(object))
    colnames(Info) = c("Tot. Self-lig.", "Self-lig. mean size", "Genome", "Sortest Self-PET",
        "Longest Self-PET", "class")
    print(knitr::kable(PET.counts, row.names = FALSE, align = c("c"), format = "markdown",
        padding = 1, output = TRUE))
    print(knitr::kable(Info[, seq_len(3)], row.names = FALSE, align = c("c"), format = "rst",
        padding = 1, output = TRUE))
    print(knitr::kable(Info[, c(4, 5)], row.names = FALSE, align = c("c"), format = "rst",
        padding = 1, output = TRUE))
}
#' @rdname summary
#' @export
#' @method summary PSFit
#' @examples
#'
#' #load Self-ligated data: (class=PSFit)
#' load(system.file('extdata', 'MACPET_psfitData.rda', package = 'MACPET'))
#' class(MACPET_psfitData)
#' summary(MACPET_psfitData)
#' summary(MACPET_psfitData,threshold=1e-5)
summary.PSFit = function(object, threshold = NULL, ...) {
    # for R check:
    FDR = NULL
    Chrom = NULL
    # prepare data:
    N = length(object)
    SLmean = S4Vectors::metadata(object)$SLmean
    hg = unique(GenomeInfoDb::genome(GenomeInfoDb::seqinfo(object)))
    MaxSize = S4Vectors::metadata(object)$MaxSize
    MinSize = S4Vectors::metadata(object)$MinSize
    PET.counts = S4Vectors::metadata(object)$Self_info
    if (!is.numeric(threshold)) {
        colnames(PET.counts) = c("Chrom", "Self-lig.", "Regions", "Peaks")
        cat("|", rep("-", 6), "Self-ligated PETs Summary", rep("-", 6), "|", sep = "")
        Info = data.frame(Total.PETs = N, Regions = sum(PET.counts$Regions), BS = sum(PET.counts$Peaks),
            SLmean = SLmean, hg = hg, MIN = paste(MinSize, "bp"), MAX = paste(MaxSize,
                "bp"), class = class(object))
        colnames(Info) = c("Tot. Self-lig.", "Regions", "Peaks", "Self-lig. mean size",
            "Genome", "Sortest Self-PET", "Longest Self-PET", "class")
        print(knitr::kable(PET.counts, row.names = FALSE, align = c("c"), format = "markdown",
            padding = 1, output = TRUE))
        print(knitr::kable(Info[, seq_len(4)], row.names = FALSE, align = c("c"),
            format = "rst", padding = 1, output = TRUE))
        print(knitr::kable(Info[, seq_len(4) + 4], row.names = FALSE, align = c("c"),
            format = "rst", padding = 1, output = TRUE))
    } else {
        Sign.peaks = S4Vectors::metadata(object)$Peaks.Info
        Sign.peaks = subset(Sign.peaks, FDR < threshold)
        Sign.peaks = plyr::ddply(Sign.peaks, plyr::.(Chrom), nrow)
        PET.counts$Sign.BS.counts = 0
        PET.counts$Sign.BS.counts[match(Sign.peaks$Chrom, PET.counts$Chrom)] = Sign.peaks$V1
        colnames(PET.counts) = c("Chrom", "Self-lig.", "Regions", "Peaks", "Sign. Peaks")
        cat("|", rep("-", 13), "Self-ligated PETs Summary", rep("-", 13), "|", sep = "")
        Info = data.frame(Total.PETs = N, Regions = sum(PET.counts$Regions), BS = sum(PET.counts$Peaks),
            Sign.BS = sum(PET.counts$`Sign. Peaks`), SLmean = SLmean, hg = hg, MIN = paste(MinSize,
                "bp"), MAX = paste(MaxSize, "bp"), class = class(object))
        colnames(Info) = c("Tot. Self-lig.", "Regions", "Peaks", "Tot. Sign. Peaks",
            "Self-lig. mean size", "Genome", "Sortest Self-PET", "Longest Self-PET",
            "class")
        print(knitr::kable(PET.counts, row.names = FALSE, align = c("c"), format = "markdown",
            padding = 1, output = TRUE))
        print(knitr::kable(Info[, seq_len(4)], row.names = FALSE, align = c("c"),
            format = "rst", padding = 1, output = TRUE))
        print(knitr::kable(Info[, seq_len(5) + 4], row.names = FALSE, align = c("c"),
            format = "rst", padding = 1, output = TRUE))
    }
}
#' @rdname summary
#' @export
#' @method summary PIntra
#' @examples
#'
#' #load Intra-chromosomal data: (class=PIntra)
#' load(system.file('extdata', 'MACPET_pintraData.rda', package = 'MACPET'))
#' class(MACPET_pintraData)
#' summary(MACPET_pintraData)
#' requireNamespace('ggplot2')
#' requireNamespace('reshape2')
#' summary(MACPET_pintraData,heatmap=TRUE)#sample data, not good heatmap plot.
summary.PIntra = function(object, heatmap = FALSE, ...) {
    # for R check:
    Var1 = NULL
    Var2 = NULL
    value = NULL
    # prepare data:
    N = length(object)
    genome = unique(GenomeInfoDb::genome(GenomeInfoDb::seqinfo(object)))
    PET.counts = S4Vectors::metadata(object)$InteractionCounts
    colnames(PET.counts) = c("Chrom", "Intra-chrom.")
    PET.counts$Chrom = as.character(PET.counts$Chrom)
    cat("|--Intra-chrom. PETs--|\n")
    cat("|-------Summary-------|")
    Info = data.frame(Total.PETs = N, genome = genome, class = class(object))
    colnames(Info) = c("Tot. Intra-chrom.", "Genome", "class")
    print(knitr::kable(PET.counts, row.names = FALSE, align = c("l", "c"), format = "markdown",
        padding = 1, output = TRUE))
    print(knitr::kable(Info, row.names = FALSE, align = c("c"), format = "rst", padding = 1,
        output = TRUE))
    if (heatmap) {
        # check package:
        if (!is.logical(heatmap)) {
            stop("heatmap has to be logical!", call. = FALSE)
        } else {
            if (!requireNamespace("ggplot2", quietly = TRUE) | !requireNamespace("reshape2",
                quietly = TRUE)) {
                stop("ggplot2 and reshape2 are needed for this
                     function to work if heatmap==TRUE.
                     Please install it.",
                  call. = FALSE)
            }
        }
        ChomNames = PET.counts$Chrom
        PET.counts = diag(PET.counts$`Intra-chrom.`)
        colnames(PET.counts) = rownames(PET.counts) = ChomNames
        plot.data = reshape2::melt(as.matrix(PET.counts))
        plot.data$value = plot.data$value/max(plot.data$value)
        RES = ggplot2::ggplot(data = plot.data, ggplot2::aes(x = Var1, y = Var2,
            fill = value)) + ggplot2::geom_tile(color = "red") + ggplot2::ggtitle("Heat-Map plot for Intra-chromosomal Counts") +
            ggplot2::xlab("Chromosome") + ggplot2::ylab("Chromosome") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        return(RES)
    }
}
#' @rdname summary
#' @export
#' @method summary PInter
#' @examples
#'
#' #load Inter-chromosomal data: (class=PInter)
#' load(system.file('extdata', 'MACPET_pinterData.rda', package = 'MACPET'))
#' class(MACPET_pinterData)
#' summary(MACPET_pinterData)
#' requireNamespace('ggplot2')
#' requireNamespace('reshape2')
#' summary(MACPET_pinterData,heatmap=TRUE)#sample data, not good heatmap plot.
summary.PInter = function(object, heatmap = FALSE, ...) {
    # for R check:
    Var1 = NULL
    Var2 = NULL
    value = NULL
    if (!is.logical(heatmap)) {
        stop("heatmap has to be logical!", call. = FALSE)
    }
    #-----------
    # prepare data:
    N = length(object)
    genome = unique(GenomeInfoDb::genome(GenomeInfoDb::seqinfo(object)))
    PET.counts = S4Vectors::metadata(object)$InteractionCounts
    PET.countsPrint = rowSums(PET.counts)
    PET.countsPrint = data.frame(Chrom = names(PET.countsPrint), Counts = PET.countsPrint)
    colnames(PET.countsPrint) = c("Chrom", "Inter-chrom.")
    PET.countsPrint$Chrom = as.character(PET.countsPrint$Chrom)
    Info = data.frame(Total.PETs = N, genome = genome, class = class(object))
    colnames(Info) = c("Tot. Inter-chrom. PETs", "Genome", "class")
    cat("|--Inter-chrom. PETs--|\n")
    cat("|-------Summary-------|")
    print(knitr::kable(PET.countsPrint, row.names = FALSE, align = c("l", "c"), format = "markdown",
        padding = 1, output = TRUE))
    print(knitr::kable(Info, row.names = FALSE, align = c("c"), format = "rst", padding = 1,
        output = TRUE))
    if (heatmap) {
        # check package:
        if (!requireNamespace("ggplot2", quietly = TRUE) | !requireNamespace("reshape2",
            quietly = TRUE)) {
            stop("ggplot2 and reshape2 are needed for this
                 function to work if heatmap==TRUE.
                 Please install it.",
                call. = FALSE)
        }
        plot.data = reshape2::melt(as.matrix(PET.counts))
        plot.data$value = plot.data$value/max(plot.data$value)
        RES = ggplot2::ggplot(data = plot.data, ggplot2::aes(x = Var1, y = Var2,
            fill = value)) + ggplot2::geom_tile(color = "red") + ggplot2::ggtitle("Heat-Map plot for Inter-chromosomal Counts") +
            ggplot2::xlab("Chromosome") + ggplot2::ylab("Chromosome") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        return(RES)
    }
}
#' @rdname summary
#' @export
#' @method summary GenomeMap
#' @examples
#'
#' #load Interaction data: (class=GenomeMap)
#' load(system.file('extdata', 'MACPET_GenomeMapData.rda', package = 'MACPET'))
#' class(MACPET_GenomeMapData)
#' summary(MACPET_GenomeMapData)
summary.GenomeMap = function(object, threshold = NULL, ...) {
    # for R check: prepare data:
    N = length(object) # total interactions
    hg = unique(GenomeInfoDb::genome(GenomeInfoDb::seqinfo(object)))
    # take object:
    object = GetSignInteractions(object = object, threshold = threshold, ReturnedAs = "GInteractions")
    object = as.data.frame(object)
    InterInteractions = length(which(object$seqnames1 != object$seqnames2))
    IntraInteractions = length(which(object$seqnames1 == object$seqnames2))
    # print:
    cat("|", rep("-", 42), "Interactions Summary", rep("-", 42), "|", sep = "")
    Info = data.frame(TotalPeaks = N, TotalIntra = IntraInteractions, TotalInter = InterInteractions,
        hg = hg, class = "GenomeMap")
    colnames(Info) = c("Tot. Peaks used", "Tot. Intra-chrom. interactions", "Tot. Inter-chrom. interactions",
        "Genome", "Class")
    print(knitr::kable(Info, row.names = FALSE, align = c("c"), format = "markdown",
        padding = 1, output = TRUE))
}
