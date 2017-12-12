#' @importFrom GenomicRanges GRanges seqinfo seqnames
#' @importFrom InteractionSet anchors GInteractions pairdist findOverlaps
#' @importFrom GenomeInfoDb seqlengths seqinfo seqlevelsInUse seqnames seqlevels genome
#' @importFrom plyr . dlply
#' @importFrom intervals Intervals interval_included
#' @importFrom S4Vectors metadata mcols queryHits
#' @importFrom GenomicAlignments readGAlignmentPairs first last
#' @importFrom Rsamtools asBam indexBam testPairedEndBam BamFile scanBamFlag ScanBamParam scanBamHeader
#' @importFrom gtools mixedsort
#' @importFrom BiocParallel bplapply
#' @importFrom methods as
#' @importFrom IRanges IRanges
#' @importClassesFrom GenomicAlignments GAlignmentPairs
#' @importClassesFrom Rsamtools ScanBamParam
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom IRanges IRanges
############################################## Main function for Stage 2:
#-------------
#-------------
Stage_2_Main_fun = function(SA_prefix, S2_image, SA_LogFile.dir, S2_AnalysisDir,
    S2_PairedData) {
    # Take time:
    Analysis.time.start = Sys.time()
    #--------------------------------------------
    #---------------Find Inter PETs:
    #--------------------------------------------
    Ninter = FindInter_fun(S2_PairedData = S2_PairedData,
                           SA_LogFile.dir = SA_LogFile.dir,
                           S2_AnalysisDir = S2_AnalysisDir,
                           SA_prefix = SA_prefix)
    #--------------------------------------------
    #---------------Classify Self/Intra:
    #--------------------------------------------
    # remove inter:
    S2_PairedData = subset(S2_PairedData, !is.na(S2_PairedData$Dist))
    # classify:
    SelfIndicator = SepSelfIntra_fun(S2_PairedData = S2_PairedData,
                                     SA_LogFile.dir = SA_LogFile.dir,
                                     SA_prefix = SA_prefix, S2_image = S2_image,
                                     S2_AnalysisDir = S2_AnalysisDir)
    #--------------------------------------------
    #---------------Find Intra PETs:
    #--------------------------------------------
    Nintra = FindIntra_fun(S2_PairedData = S2_PairedData,
                           SelfIndicator = SelfIndicator,
                           SA_LogFile.dir = SA_LogFile.dir,
                           S2_AnalysisDir = S2_AnalysisDir,
                           SA_prefix = SA_prefix)
    #--------------------------------------------
    #---------------Find Self PETs:
    #--------------------------------------------
    Nself = FindSelf_fun(S2_PairedData = S2_PairedData,
                         SelfIndicator = SelfIndicator,
                         SA_LogFile.dir = SA_LogFile.dir,
                         S2_AnalysisDir = S2_AnalysisDir,
                         SA_prefix = SA_prefix)
    #--------------------------------------------
    #---------------plot and print:
    #--------------------------------------------
    LogFile = list()
    LogFile[[1]] = paste("===========>PET statistics<==========\n")
    LogFile[[2]] = paste("Total Self-ligated PETs:", Nself, "\n")
    LogFile[[3]] = paste("Total Intra-chromosomal PETs:", Nintra, "\n")
    LogFile[[4]] = paste("Total Inter-chromosomal PETs:", Ninter, "\n")
    LogFile[[5]] = "=====================================\n"
    for (lf in seq_len(5)) cat(LogFile[[lf]])
    for (lf in seq_len(5)) write(LogFile[[lf]], file = SA_LogFile.dir,
                                 append = TRUE)
    if (S2_image) {
        Get_image_S2_P3_fun(S2_AnalysisDir = S2_AnalysisDir,
                            SA_prefix = SA_prefix, Nself = Nself,
                            Nintra = Nintra, Ninter = Ninter)
    }
    #-----------------------------------------------------------------------#
    LogFile = paste("Stage 2 is done!\n")
    cat(LogFile)
    write(LogFile, file = SA_LogFile.dir, append = TRUE)
    LogFile = paste("Analysis results for stage 2 are in:\n",
                    S2_AnalysisDir, "\n")
    cat(LogFile)
    write(LogFile, file = SA_LogFile.dir, append = TRUE)
    # take time:
    Analysis.time.end = Sys.time()
    Total.Time = Analysis.time.end - Analysis.time.start
    LogFile = paste("Total stage 2 time: ", Total.Time, " ",
                    units(Total.Time), "\n",
                    sep = "")
    cat(LogFile)
    write(LogFile, file = SA_LogFile.dir, append = TRUE)
}
# done
#-------------
#-------------
############## Loading data functions:
#-------------
#-------------
# function for loading BAM file is stage 2 is run only, loaded at InputChecks
LoadBAM_FromInputChecks_fun = function(SA_AnalysisDir, S2_AnalysisDir, SA_prefix,
    S2_PairedEndBAMpath, Format, SA_LogFile.dir, S2_BlackList, S2_image) {
    #------------
    #--Find the file and maybe convert sam to bam:
    #------------
    if (Format == "sam") {
        # print:
        cat("SAM format detected..\n")
        cat("Converting SAM to BAM format and creating BAM index..\n")
        # make folder to save the output bam, should be done at stage 1:
        S1_AnalysisDir = file.path(SA_AnalysisDir, "S1_results")
        if (!dir.exists(S1_AnalysisDir))
            dir.create(S1_AnalysisDir)
        # convert:
        PairedEndBAMpath = file.path(S1_AnalysisDir, paste(SA_prefix, "_Paired_end",
            sep = ""))
        suppressWarnings(Rsamtools::asBam(file = S2_PairedEndBAMpath,
                                          destination = PairedEndBAMpath,
                                          overwrite = FALSE,
                                          indexDestination = TRUE))
        # update path:
        S2_PairedEndBAMpath = paste(PairedEndBAMpath, ".bam", sep = "")
    } else if (Format == "bam") {
        # check if index exists:
        if (!file.exists(paste(S2_PairedEndBAMpath, ".bai", sep = ""))) {
            cat("Creating BAM index...\n")
            invisible(Rsamtools::indexBam(file = S2_PairedEndBAMpath))
        }
    }
    #------------
    #---Test if paired-end data:
    #------------
    IsPairedEnd = Rsamtools::testPairedEndBam(file = S2_PairedEndBAMpath,
                                              index = S2_PairedEndBAMpath)
    if (!IsPairedEnd) {
        stop("S2_PairedEndBAMpath bam file is not paired-end file!", call. = FALSE)
    } else {
        cat("S2_PairedEndBAMpath bam file is paired-end file.\n")
    }
    #------------
    #---Load the data:
    #------------
    S2_PairedData = LoadBAM_FromMACPETUlt_fun(S2_PairedEndBAMpath = S2_PairedEndBAMpath,
        SA_LogFile.dir = SA_LogFile.dir, S2_BlackList = S2_BlackList, S2_image = S2_image,
        S2_AnalysisDir = S2_AnalysisDir, SA_prefix = SA_prefix)
    return(S2_PairedData)
}
# done
#-------------
#-------------
# function for loading BAM file is stage 1:2 are in sequence, used in MACPETUlt
# function
LoadBAM_FromMACPETUlt_fun = function(S2_PairedEndBAMpath, SA_LogFile.dir, S2_BlackList,
    S2_image, S2_AnalysisDir, SA_prefix) {
    # then I know that the BAM file is paired correctly and it is bam format
    #-------------
    # create directory
    #-------------
    if (!dir.exists(S2_AnalysisDir))
        dir.create(S2_AnalysisDir)
    #-------------
    # check the header and return black list:
    #-------------
    S2_BL_genome = Check_BAM_Header_fun(S2_PairedEndBAMpath = S2_PairedEndBAMpath,
        S2_BlackList = S2_BlackList, SA_LogFile.dir = SA_LogFile.dir)
    #-------------
    # Load data:
    #-------------
    cat("Loading PET data...\n")
    TotPETs = Rsamtools::countBam(file = S2_PairedEndBAMpath)
    TotPETs = TotPETs$records/2  #get pets
    #-------------
    # BAM instance:
    #-------------
    bamfile = Rsamtools::BamFile(file = S2_PairedEndBAMpath,
                                 index = S2_PairedEndBAMpath, asMates = TRUE)
    FlagsParam = Rsamtools::scanBamFlag(isPaired = TRUE, isUnmappedQuery = FALSE,
        hasUnmappedMate = FALSE, isSecondaryAlignment = FALSE,
        isNotPassingQualityControls = FALSE, isDuplicate = FALSE)
    ReadParam = Rsamtools::ScanBamParam(flag = FlagsParam)
    S2_PairedData = GenomicAlignments::readGAlignmentPairs(file = bamfile, use.names = FALSE,
        with.which_label = FALSE, strandMode = 1, param = ReadParam)
    # check PCR and before black list:
    TotPETB_BL = length(S2_PairedData)
    TotPCR = TotPETs - TotPETB_BL
    # call Ginteractions convert and remove black list too.:
    S2_PairedData = GInteractionsCovnert_fun(S2_PairedData = S2_PairedData,
                                             S2_BlackList = S2_BL_genome$S2_BlackList,
                                             PselfConvert = FALSE)
    TotPETfinal = length(S2_PairedData)
    TotBL = TotPETB_BL - TotPETfinal
    # set the genome:
    GenomeInfoDb::genome(GenomeInfoDb::seqinfo(S2_PairedData)) = S2_BL_genome$S1_genome
    #------------
    # print
    #------------
    NPETfinal100 = TotPETfinal/TotPETs * 100
    NPCR100 = TotPCR/TotPETs * 100
    NBL100 = TotBL/TotPETs * 100
    LogFile = list()
    LogFile[[1]] = paste("===========>PET statistics<==========\n")
    LogFile[[2]] = paste("Total PETs in data:", TotPETs, "(", 100, "%)", "\n")
    LogFile[[3]] = paste("Total PCR replicates:", TotPCR, "(", NPCR100, "%)", "\n")
    LogFile[[4]] = paste("Total Black-listed PETs:", TotBL, "(", NBL100, "%)", "\n")
    LogFile[[5]] = paste("Total valid PETs left:", TotPETfinal, "(", NPETfinal100,
        "%)", "\n")
    for (lf in seq_len(5)) cat(LogFile[[lf]])
    for (lf in seq_len(5)) write(LogFile[[lf]], file = SA_LogFile.dir, append = TRUE)
    #-------------
    # plot:
    #-------------
    if (S2_image) {
        Get_image_S2_P1_fun(S2_AnalysisDir = S2_AnalysisDir, SA_prefix = SA_prefix,
            NPETfinal100 = NPETfinal100, NPCR100 = NPCR100, NBL100 = NBL100)
    }
    #-------------
    # return:
    #-------------
    return(S2_PairedData)
}
# done
#-------------
#-------------
# function for testing the header of the bam file:
Check_BAM_Header_fun = function(S2_PairedEndBAMpath, S2_BlackList, SA_LogFile.dir) {
    cat("Checking the bam file header for the genome....")
    HeaderBAM = Rsamtools::scanBamHeader(file = S2_PairedEndBAMpath, what = "text")
    HeaderBAM = HeaderBAM[[1]]$text
    HeaderBAM = HeaderBAM[which(names(HeaderBAM) %in% "@SQ")]
    if (length(HeaderBAM) == 0) {
        stop("The bam file is missing the header section!", call. = FALSE)
    }
    HeaderBAM = do.call(rbind, HeaderBAM)
    # check if any header is missing:
    SNpos = which(grepl("SN:", HeaderBAM[1, ]))
    LNpos = which(grepl("LN:", HeaderBAM[1, ]))
    ASpos = which(grepl("AS:", HeaderBAM[1, ]))
    if (length(SNpos) == 0) {
        stop("Bam file header is missing the SN entry!\n", call. = FALSE)
    }
    if (length(LNpos) == 0) {
        stop("Bam file header is missing the LN entry!\n", call. = FALSE)
    }
    if (length(ASpos) == 0 && isTRUE(S2_BlackList)) {
        LogFile = "The bam file is missing the 'AS' genome header. No black-listed regions will be removed from the data\n"
        warning(LogFile)
        write(paste("WARNING:", LogFile), file = SA_LogFile.dir, append = TRUE)
        S2_BlackList = NULL
        S1_genome = NA
    } else if (length(ASpos) != 0) {
        # get the genome:
        S1_genome = HeaderBAM[1, ASpos]
        S1_genome = strsplit(S1_genome, "AS:")
        S1_genome = unlist(S1_genome)[2]
        if (!S1_genome %in% names(sysdata) && isTRUE(S2_BlackList)) {
            LogFile = paste("The bam file genome: ", S1_genome, " is not one of the following: ",
                paste(names(sysdata), collapse = "/"), ". No black listed regions will be removed!\n")
            warning(LogFile)
            write(paste("WARNING:", LogFile), file = SA_LogFile.dir, append = TRUE)
            S2_BlackList = NULL
        } else if (S1_genome %in% names(sysdata) && isTRUE(S2_BlackList)) {
            S2_BlackList = sysdata[[S1_genome]]
        } else if (isFALSE(S2_BlackList)) {
            S2_BlackList = NULL
        }
    }
    cat("OK\n")
    return(list(S1_genome = S1_genome, S2_BlackList = S2_BlackList))
}
# Done
#-------------
#-------------
# create GI object and remove black list
GInteractionsCovnert_fun = function(S2_PairedData, S2_BlackList, PselfConvert) {
    #------------
    # break returned in anchors:
    #------------
    if (!PselfConvert) {
        # then the function called from here
        Anchor1 = GenomicAlignments::first(x = S2_PairedData, real.strand = FALSE)
        Anchor1 = methods::as(Anchor1, "GRanges")
        # to tags
        Anchor2 = GenomicAlignments::last(x = S2_PairedData, real.strand = FALSE)
        Anchor2 = methods::as(Anchor2, "GRanges")
    } else {
        # then the function called by ConvertToPSelf: from tags
        Anchor1 = InteractionSet::anchors(S2_PairedData, type = "first")
        # to tags
        Anchor2 = InteractionSet::anchors(S2_PairedData, type = "second")
    }
    #------------
    # make object:
    #------------
    S2_PairedData = InteractionSet::GInteractions(anchor1 = Anchor1, anchor2 = Anchor2)
    #--------------------------------------------
    #-------Remove black listed regions:
    #--------------------------------------------
    S2_PairedData = BlackListCorrection_fun(S2_PairedData = S2_PairedData, S2_BlackList = S2_BlackList)
    #------------
    # Find spans/Dist:
    #------------
    S2_PairedData$Dist = InteractionSet::pairdist(S2_PairedData, type = "span")
    return(S2_PairedData)
}
# done
#-------------
#-------------
# function for removing black listed PETs
BlackListCorrection_fun = function(S2_PairedData, S2_BlackList) {
    if (!is.null(S2_BlackList)) {
        # then it is a GRanges object
        BlackListed = InteractionSet::findOverlaps(query = S2_PairedData, subject = S2_BlackList,
            maxgap = -1L, minoverlap = 0L, type = c("any"), select = c("all"), ignore.strand = TRUE,
            use.region = "both")
        BlackListed = S4Vectors::queryHits(BlackListed)
        #---------------
        #-----reduce:
        #---------------
        if (length(BlackListed) != 0)
            S2_PairedData = S2_PairedData[-BlackListed]
        if (length(S2_PairedData) == 0) {
            stop("The data contained only black-listed regions!", call. = FALSE)
        }
    }
    return(S2_PairedData)
}
# done
#-------------
#-------------
# function for plotting for stage 2 part 1: PET totals etc.
Get_image_S2_P1_fun = function(S2_AnalysisDir, SA_prefix, NPETs, NPETfinal100, NPCR100,
    NBL100) {
    # Rcheck:
    Value = NULL
    Kind = NULL
    # image dir:
    S2_P1_image_dir = file.path(S2_AnalysisDir, paste(SA_prefix, "_stage_2_p1_image.jpg",
        sep = ""))
    #-------------
    # create data:
    #-------------
    S2_imagedata_1 = data.frame(Kind = c(paste("Final PETs (", round(NPETfinal100,
        digits = 1), "%)", sep = ""), paste("PCR PETs (", round(NPCR100, digits = 1),
        "%)", sep = ""), paste("Black-listed PETs (", round(NBL100, digits = 1),
        "%)", sep = "")), Value = c(round(NPETfinal100, digits = 1), round(NPCR100,
        digits = 1), round(NBL100, digits = 1)))
    #-------------
    # plot the split:
    #-------------
    S2_image_p1 = ggplot2::ggplot(S2_imagedata_1, ggplot2::aes(x = "", y = Value,
        fill = factor(Kind))) + ggplot2::geom_bar(width = 1, stat = "identity") +
        ggplot2::coord_polar(theta = "y") + ggplot2::theme(axis.title = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 20, color = "black"), legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 17), axis.text = ggplot2::element_blank(),
        legend.position = "bottom", legend.direction = "vertical", axis.ticks = ggplot2::element_blank()) +
        ggplot2::ggtitle("Pie chart for PET removal") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::scale_fill_brewer(palette = "Dark2")
    # save:
    ggplot2::ggsave(plot = S2_image_p1, file = S2_P1_image_dir, scale = 2)
}
# done
#-------------
#-------------
############################################## PET classification functions:
#-------------
#-------------
FindInter_fun = function(S2_PairedData, SA_LogFile.dir, S2_AnalysisDir, SA_prefix) {
    cat("=====================================\n")
    cat("Separating Inter-chromosomal data...")
    #------------
    # keep inter pets:
    #------------
    pinterData = subset(S2_PairedData, is.na(S2_PairedData$Dist))
    Ninter = length(pinterData)
    #------------
    if (Ninter != 0) {
        #------------
        # remove distance
        #------------
        pinterData$Dist = NULL
        #------------
        # reduce the DataSeqinfo based on the chromosomes inside the data:
        #------------
        Anchor1 = InteractionSet::anchors(pinterData, type = "first")
        Anchor2 = InteractionSet::anchors(pinterData, type = "second")
        Anchor1LevelsUsed = GenomeInfoDb::seqlevelsInUse(Anchor1)
        Anchor2LevelsUsed = GenomeInfoDb::seqlevelsInUse(Anchor2)
        LevelsUsed = unique(c(Anchor1LevelsUsed, Anchor2LevelsUsed))
        LevelsUsed = gtools::mixedsort(LevelsUsed)
        #------------
        # update Anchors:
        #------------
        GenomeInfoDb::seqlevels(Anchor1) = LevelsUsed
        GenomeInfoDb::seqlevels(Anchor2) = LevelsUsed
        pinterData = InteractionSet::GInteractions(Anchor1, Anchor2)
        #------------
        # find counts:
        #------------
        Inter.Counts = data.frame(seqnames1 = GenomeInfoDb::seqnames(Anchor1), seqnames2 = GenomeInfoDb::seqnames(Anchor2),
            stringsAsFactors = FALSE)
        Inter.Counts = table(Inter.Counts)
        Inter.Counts = as.data.frame.matrix(Inter.Counts, stringsAsFactors = FALSE)
        # order:
        Order.col = match(gtools::mixedsort(colnames(Inter.Counts)), colnames(Inter.Counts))
        Inter.Counts = Inter.Counts[, Order.col]
        Order.row = match(gtools::mixedsort(rownames(Inter.Counts)), rownames(Inter.Counts))
        Inter.Counts = Inter.Counts[Order.row, ]
        S4Vectors::metadata(pinterData) = list(InteractionCounts = Inter.Counts)
        #------------
        # save:
        #------------
        class(pinterData) = "PInter"
        NamepinterData = paste(SA_prefix, "_pinterData", sep = "")
        assign(NamepinterData, pinterData)  #assign value.
        save(list = NamepinterData, file = file.path(S2_AnalysisDir, NamepinterData))
    } else {
        warning("Inter-chromosomal data is empty!")
        write("WARNING: Inter-chromosomal data is empty!\n", file = SA_LogFile.dir,
            append = TRUE)
    }
    cat("Done\n")
    write("=====================================\n", file = SA_LogFile.dir, append = TRUE)
    return(Ninter)
}
# done
#-------------
#-------------
# Function For SelfIntra classification using the elbow method:
SepSelfIntra_fun = function(S2_PairedData, SA_LogFile.dir, SA_prefix, S2_image, S2_AnalysisDir) {
    # global variables for Rcheck:
    Freq = NULL
    Size = NULL
    PointDist = NULL
    logSize = NULL
    #--------------------------
    cat("Finding Self-Intra cut-off...")
    #------------
    # Find Elbow Point:
    #------------
    MaxSpan = max(S2_PairedData$Dist)
    MinSpan = min(S2_PairedData$Dist)
    SeqSpan = seq(from = 0, to = MaxSpan, by = 100)
    SeqSpan[which(SeqSpan == max(SeqSpan))] = MaxSpan
    IntSpan = cbind(SeqSpan[-length(SeqSpan)], SeqSpan[-1])
    IntSpan = intervals::Intervals(IntSpan, closed = rep(TRUE, 2))
    DataInt = intervals::Intervals(cbind(S2_PairedData$Dist, S2_PairedData$Dist),
        closed = rep(TRUE, 2))
    IntIncl = intervals::interval_included(IntSpan, DataInt)
    SpanDF = data.frame(Size = SeqSpan[-1], Freq = lengths(IntIncl))
    SpanDF = subset(SpanDF, Freq != 0)
    SpanDF$logSize = log(SpanDF$Size)
    #---take the points for finding the line:
    point1 = subset(SpanDF, Freq == max(Freq))  #max frequency point
    point2 = subset(SpanDF, Size == max(Size))  #last point
    #----find slope and intercept:
    slope = (point2$Freq - point1$Freq)/(point2$logSize - point1$logSize)  #slope
    intercept = point1$Freq - slope * point1$logSize
    #-----find distance to line
    SpanDF$PointDist = abs(-slope * SpanDF$logSize + 1 * SpanDF$Freq - intercept)/sqrt(slope^2 +
        1)
    # keep only dist to those after the pean-top
    SpanDFSub = subset(SpanDF, Size >= point1$Size)
    ElbowPoint = subset(SpanDFSub, PointDist == max(PointDist))
    # ElbowPoint
    ElbowPoint = subset(ElbowPoint, Size == min(Size))
    #------------
    # find Self_indicator:
    #------------
    SelfIndicator = which(S2_PairedData$Dist <= ElbowPoint$Size)
    MAXcut = max(S2_PairedData$Dist[SelfIndicator])
    MINcut = min(S2_PairedData$Dist[SelfIndicator])
    SelfBorder = data.frame(MAX = MAXcut, MIN = MINcut)
    #------------
    # print::
    #------------
    cat("Done\n")
    LogFile = paste("Self-ligated cut-off at:", SelfBorder$MAX, "bp\n")
    cat(LogFile)
    write(LogFile, file = SA_LogFile.dir, append = TRUE)
    #------------
    # image:
    #------------
    if (S2_image) {
        Get_image_S2_P2_fun(S2_AnalysisDir = S2_AnalysisDir, SA_prefix = SA_prefix,
            SpanDF = SpanDF, ElbowPoint = ElbowPoint, SelfBorder = SelfBorder)
    }
    return(SelfIndicator)
}
# done
#-------------
#-------------
# function for plotting for stage 2 part 2: self0ligated cut-off
Get_image_S2_P2_fun = function(S2_AnalysisDir, SA_prefix, SpanDF, ElbowPoint, SelfBorder) {
    # Rcheck:
    Size = NULL
    Freq = NULL
    # image dir:
    S2_P2_image_dir = file.path(S2_AnalysisDir, paste(SA_prefix, "_stage_2_p2_image.jpg",
        sep = ""))
    # Plot the elbow:
    S2_image_p2 = ggplot2::ggplot(SpanDF, ggplot2::aes(x = Size, y = Freq)) + ggplot2::geom_rect(ggplot2::aes(xmin = SpanDF$Size -
        49, xmax = SpanDF$Size + 49, ymin = 0, ymax = SpanDF$Freq), size = 0.3, fill = "grey69",
        color = "black") + ggplot2::geom_line(color = "blue", size = 0.6) + ggplot2::geom_vline(xintercept = ElbowPoint$Size,
        linetype = "dashed", color = "red", size = 0.6) + ggplot2::scale_x_log10(labels = function(co) round(log10(co)), expand = c(0,0)) +
        ggplot2::ggtitle(paste("Elbow-point cut-off for Self/Intra Pets at: ", SelfBorder$MAX,
            " bp")) + ggplot2::xlab("Sorted log-sizes of PET") + ggplot2::ylab("Frequency") +
        ggplot2::theme_bw() + ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
        panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank()) + ggplot2::theme(axis.text = ggplot2::element_text(size = 15,
        color = "black"), axis.title = ggplot2::element_text(size = 18, color = "black"),
        plot.title = ggplot2::element_text(size = 18, color = "black")) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+
        ggplot2::scale_y_continuous(expand = c(0,0))
    # save:
    ggplot2::ggsave(plot = S2_image_p2, file = S2_P2_image_dir, scale = 2)
}
# done
#-------------
#-------------
# Function for Intra PETs
FindIntra_fun = function(S2_PairedData, SelfIndicator, SA_LogFile.dir, S2_AnalysisDir,
    SA_prefix) {
    cat("Separating Intra-chromosomal data...")
    #------------
    # keep intra:
    #------------
    pintraData = S2_PairedData[-SelfIndicator, ]
    Nintra = length(pintraData)
    if (Nintra != 0) {
        #------------
        # keep to save again
        #------------
        Mcols = S4Vectors::mcols(pintraData)
        #------------
        # reduce the DataSeqinfo based on the chromosomes inside the data:
        #------------
        Anchor1 = InteractionSet::anchors(pintraData, type = "first")
        Anchor2 = InteractionSet::anchors(pintraData, type = "second")
        Anchor1LevelsUsed = GenomeInfoDb::seqlevelsInUse(Anchor1)
        Anchor2LevelsUsed = GenomeInfoDb::seqlevelsInUse(Anchor2)
        LevelsUsed = unique(c(Anchor1LevelsUsed, Anchor2LevelsUsed))
        LevelsUsed = gtools::mixedsort(LevelsUsed)
        #------------
        # update Anchors:
        #------------
        GenomeInfoDb::seqlevels(Anchor1) = LevelsUsed
        GenomeInfoDb::seqlevels(Anchor2) = LevelsUsed
        pintraData = InteractionSet::GInteractions(Anchor1, Anchor2)
        S4Vectors::mcols(pintraData) = Mcols  #return
        #------------
        # find counts
        #------------
        Intra.counts = table(GenomeInfoDb::seqnames(Anchor1))
        ChrNames = names(Intra.counts)
        CountValues = as.numeric(Intra.counts)
        Intra.counts = data.frame(Chrom = ChrNames, Counts = CountValues, stringsAsFactors = FALSE)
        # sort:
        Intra.counts = Intra.counts[match(gtools::mixedsort(ChrNames), ChrNames),
            ]
        S4Vectors::metadata(pintraData) = list(InteractionCounts = Intra.counts)
        #------------
        # remove Dist variable:
        #------------
        pintraData$Dist = NULL
        #------------
        # change class and save:
        #------------
        class(pintraData) = "PIntra"
        NamepintraData = paste(SA_prefix, "_pintraData", sep = "")
        assign(NamepintraData, pintraData)  #assign value.
        save(list = NamepintraData, file = file.path(S2_AnalysisDir, NamepintraData))
    } else {
        warning("Intra-chromosomal data is empty!")
        write("WARNING: Intra-chromosomal data is empty!\n", file = SA_LogFile.dir,
            append = TRUE)
    }
    cat("Done\n")
    return(Nintra)
}
#-------------
#-------------
# Function for Self PETs
FindSelf_fun = function(S2_PairedData, SelfIndicator, SA_LogFile.dir, S2_AnalysisDir,
    SA_prefix) {
    cat("Separating Self-ligated data...")
    #------------
    # keep self:
    #------------
    pselfData = S2_PairedData[SelfIndicator, ]
    Nself = length(pselfData)
    if (Nself != 0) {
        Mcols = S4Vectors::mcols(pselfData)
        #------------
        # reduce the DataSeqinfo based on the chromosomes inside the data:
        #------------
        Anchor1 = InteractionSet::anchors(pselfData, type = "first")
        Anchor2 = InteractionSet::anchors(pselfData, type = "second")
        Anchor1LevelsUsed = GenomeInfoDb::seqlevelsInUse(Anchor1)
        Anchor2LevelsUsed = GenomeInfoDb::seqlevelsInUse(Anchor2)
        LevelsUsed = unique(c(Anchor1LevelsUsed, Anchor2LevelsUsed))
        LevelsUsed = gtools::mixedsort(LevelsUsed)
        # update Anchors:
        GenomeInfoDb::seqlevels(Anchor1) = LevelsUsed
        GenomeInfoDb::seqlevels(Anchor2) = LevelsUsed
        pselfData = InteractionSet::GInteractions(Anchor1, Anchor2)
        S4Vectors::mcols(pselfData) = Mcols  #return mcols infor
        #------------
        # find self-ligated length:
        #------------
        SLmean = round(mean(pselfData$Dist))
        #------------
        # Take info:
        #------------
        Self_info = table(GenomicRanges::seqnames(InteractionSet::anchors(pselfData,
            type = "first")))
        Self_info = data.frame(Chrom = as.character(names(Self_info)), PET.counts = as.numeric(Self_info))
        S4Vectors::metadata(pselfData) = list(Self_info = Self_info)
        MaxSize = max(pselfData$Dist)
        MinSize = min(pselfData$Dist)
        pselfData$Dist = NULL  #remove distance
        S4Vectors::metadata(pselfData)$SLmean = SLmean
        S4Vectors::metadata(pselfData)$MaxSize = MaxSize
        S4Vectors::metadata(pselfData)$MinSize = MinSize
        #------------
        # change class and save:
        #------------
        class(pselfData) = "PSelf"
        NamepselfData = paste(SA_prefix, "_pselfData", sep = "")
        assign(NamepselfData, pselfData)  #assign value.
        save(list = NamepselfData, file = file.path(S2_AnalysisDir, NamepselfData))
    } else {
        warning("Self-ligated data is empty!")
        if (!is.null(SA_LogFile.dir)) {
            write("WARNING: Self-ligated data is empty!\n", file = SA_LogFile.dir,
                append = TRUE)
        }
    }
    cat("Done\n")
    LogFile = paste("Self-ligated mean length: ", SLmean, "\n")
    cat(LogFile)
    if (!is.null(SA_LogFile.dir)) {
        write(LogFile, file = SA_LogFile.dir, append = TRUE)
    }
    return(Nself)
}
# done
#-------------
#-------------
# function for plotting for stage 2 part 3: PET classification
Get_image_S2_P3_fun = function(S2_AnalysisDir, SA_prefix, Nself, Nintra, Ninter) {
    # Rcheck:
    Value = NULL
    Kind = NULL
    # find percentage:
    Ntot = Nself + Nintra + Ninter
    Nself100 = Nself/Ntot * 100
    Nintra100 = Nintra/Ntot * 100
    Ninter100 = Ninter/Ntot * 100
    # image dir:
    S2_P3_image_dir = file.path(S2_AnalysisDir, paste(SA_prefix, "_stage_2_p3_image.jpg",
        sep = ""))
    #-------------
    # create data:
    #-------------
    S2_imagedata_3 = data.frame(Kind = c(paste("Self-ligated PETs (", round(Nself100,
        digits = 1), "%)", sep = ""), paste("Intra-chromosomal PETs (", round(Nintra100,
        digits = 1), "%)", sep = ""), paste("Inter-chromosomal PETs (", round(Ninter100,
        digits = 1), "%)", sep = "")), Value = c(round(Nself100, digits = 1), round(Nintra100,
        digits = 1), round(Ninter100, digits = 1)))
    #-------------
    # plot the split:
    #-------------
    S2_image_p3 = ggplot2::ggplot(S2_imagedata_3, ggplot2::aes(x = "", y = Value,
        fill = factor(Kind))) + ggplot2::geom_bar(width = 1, stat = "identity") +
        ggplot2::coord_polar(theta = "y") + ggplot2::theme(axis.title = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 18, color = "black"), legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 17), axis.text = ggplot2::element_blank(),
        legend.position = "bottom", legend.direction = "vertical", axis.ticks = ggplot2::element_blank()) +
        ggplot2::ggtitle("Pie chart for PET classification") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::scale_fill_brewer(palette = "Dark2")
    # save:
    ggplot2::ggsave(plot = S2_image_p3, file = S2_P3_image_dir, scale = 2)
}
# done
#-------------
#-------------
