#' @importFrom GenomicRanges GRanges seqinfo seqnames start end strand
#' @importFrom InteractionSet anchors GInteractions pairdist findOverlaps
#' @importFrom GenomeInfoDb seqlengths seqinfo seqlevelsInUse seqnames seqlevels
#' @importFrom plyr . dlply
#' @importFrom intervals Intervals interval_overlap interval_included
#' @importFrom S4Vectors metadata unique mcols queryHits
#' @importFrom GenomicAlignments readGAlignmentPairs first last
#' @importFrom Rsamtools asBam indexBam testPairedEndBam BamFile scanBamFlag ScanBamParam
#' @importFrom gtools mixedsort
#' @importFrom BiocParallel bplapply
#' @importFrom methods as
#' @importFrom IRanges IRanges
####################### 
#-------------
#-------------
# Main function for Stage 0:
PETClassification_fun = function(DataDir, DataFile, AnalysisDir, fileSelf, fileIntra, 
    fileInter, PopImage, BlackList, LogFile.dir, Format, GenInfo, ChromLengths) {
    #--------------------------------------------
    #---------------Take and reorder Input:
    #--------------------------------------------
    # Take time:
    Analysis.time.start = Sys.time()
    # check inputs:
    InArg = list(DataDir = DataDir, DataFile = DataFile, AnalysisDir = AnalysisDir, 
        fileSelf = fileSelf, fileIntra = fileIntra, fileInter = fileInter, PopImage = PopImage, 
        BlackList = BlackList, LogFile.dir = LogFile.dir, Format = Format, GenInfo = GenInfo, 
        ChromLengths = ChromLengths)
    #--------------------------------------------
    #---------------Load input file:
    #--------------------------------------------
    if (!is.na(InArg$Format[length(InArg$Format)])) {
        # then load data
        ChIAPET = LoadData_fun(ArgLoadData = InArg)
    } else {
        # else the object is loaded already
        ChIAPET = InArg$DataFile
    }
    #--------------------------------------------
    #---------------Create GInteractions object:
    #--------------------------------------------
    PairedData = GInteractionsCovnert_fun(ChIAPET = ChIAPET, ArgPairedData = InArg)
    #--------------------------------------------
    #---------------Find Inter PETs:
    #--------------------------------------------
    pinterData = FindInter_fun(x = PairedData, ArgFindInter = InArg)
    LogFile = paste("Total ", length(pinterData), "Inter-chromosomal PETs found.\n")
    cat(LogFile)
    write(LogFile, file = InArg$LogFile.dir, append = TRUE)  #write in log file.
    if (length(pinterData) != 0) {
        cat("Saving Inter-chromosomal data...\n")
        assign(InArg$fileInter, pinterData)  #assign value.
        save(list = InArg$fileInter, file = file.path(InArg$AnalysisDir, InArg$fileInter))  #save
    } else {
        # dont save any
        LogFile = "Empty Inter-PETs are not saved.\n"
        cat(LogFile)
        write(LogFile, file = InArg$LogFile.dir, append = TRUE)  #write in log file.
    }
    #--------------------------------------------
    #---------------Classify Self/Intra:
    #--------------------------------------------
    # remove inter:
    SIPairedData = subset(PairedData, !is.na(PairedData$Dist))
    # classify:
    SelfIndicator = FindSelfIntra_fun(x = SIPairedData, ArgFindSelfIntra = InArg)
    #--------------------------------------------
    #---------------Find Intra PETs:
    #--------------------------------------------
    pintraData = FindIntra_fun(x = SIPairedData, SelfIndicator = SelfIndicator, ArgFindIntra = InArg)
    LogFile = paste("Total ", length(pintraData), "Intra-chromosomal PETs found.\n")
    cat(LogFile)
    write(LogFile, file = InArg$LogFile.dir, append = TRUE)  #write in log file.
    if (length(pintraData) != 0) {
        cat("Saving Intra-chromosomal data...\n")
        assign(InArg$fileIntra, pintraData)  #assign value.
        save(list = InArg$fileIntra, file = file.path(InArg$AnalysisDir, InArg$fileIntra))  #save
    } else {
        # dont save any
        LogFile = "Empty Intra-PETs are not saved.\n"
        cat(LogFile)
        write(LogFile, file = InArg$LogFile.dir, append = TRUE)  #write in log file.
    }
    #--------------------------------------------
    #---------------Find Self PETs:
    #--------------------------------------------
    pselfData = FindSelf_fun(x = SIPairedData, SelfIndicator = SelfIndicator, ArgFindSelf = InArg, 
        PSelf.Convert = FALSE)
    LogFile = paste("Total ", length(pselfData), "Self-ligated PETs found.\n")
    cat(LogFile)
    write(LogFile, file = InArg$LogFile.dir, append = TRUE)  #write in log file.
    cat("Saving Self-ligated data...\n")
    assign(InArg$fileSelf, pselfData)  #assign value.
    save(list = InArg$fileSelf, file = file.path(InArg$AnalysisDir, InArg$fileSelf))  #save
    #-----------------------------------------------------------------------#
    # take time:
    Analysis.time.end = Sys.time()
    Total.Time = Analysis.time.end - Analysis.time.start
    LogFile = paste("Total splitting time: ", Total.Time, " ", units(Total.Time), 
        "\n", sep = "")
    cat(LogFile)
    write(LogFile, file = InArg$LogFile.dir, append = TRUE)  #write in log file.
    return("The split is done!\n")
}
#-------------
#-------------
####################### Loading data:
#-------------
#-------------
# function for loading BAM/SAM
LoadData_fun = function(ArgLoadData) {
    #------------
    #--Find the file and maybe convert sam to bam:
    #------------
    if (ArgLoadData$Format[length(ArgLoadData$Format)] == "sam") {
        cat("Converting SAM to BAM format and creating BAM index..\n")
        file = file.path(ArgLoadData$DataDir, ArgLoadData$DataFile)
        destination = file.path(ArgLoadData$DataDir, ArgLoadData$Format[1])
        Rsamtools::asBam(file = file, destination = destination, overwrite = FALSE, 
            indexDestination = TRUE)
        ArgLoadData$DataFile = paste(ArgLoadData$Format[1], ".bam", sep = "")
    } else if (ArgLoadData$Format[length(ArgLoadData$Format)] == "bam") {
        cat("Creating BAM index...\n")
        Rsamtools::indexBam(file.path(ArgLoadData$DataDir, ArgLoadData$DataFile))
    }
    #------------
    #---Test if paired-end data:
    #------------
    file = file.path(ArgLoadData$DataDir, ArgLoadData$DataFile)
    PairedEnd = Rsamtools::testPairedEndBam(file = file, index = file)
    if (!PairedEnd) {
        stop("BAM file is not paired-end file!", call. = FALSE)
    } else {
        cat("BAM file is paired-end file.\n")
    }
    #------------
    #---Load the data, by default it will only load PETs with both ends mapped
    #------------
    cat("Loading PET data...\n")
    bamfile = Rsamtools::BamFile(file = file, index = file, asMates = TRUE)
    FlagsParam = Rsamtools::scanBamFlag(isPaired = TRUE, isUnmappedQuery = FALSE, 
        hasUnmappedMate = FALSE, isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE, 
        isDuplicate = FALSE)
    ReadParam = Rsamtools::ScanBamParam(flag = FlagsParam)
    ChIAPET = GenomicAlignments::readGAlignmentPairs(file = bamfile, index = file, 
        use.names = FALSE, with.which_label = FALSE, strandMode = 1, param = ReadParam)
    #------------
    #--Remove index:
    #------------
    cat("Removing indexbam...\n")
    remove_bai = list.files(ArgLoadData$DataDir, pattern = ".bai")
    unlink(x = file.path(ArgLoadData$DataDir, remove_bai), recursive = TRUE)
    LogFile = paste("Number of PETs in data:", length(ChIAPET), "\n")
    cat(LogFile)
    write(LogFile, file = ArgLoadData$LogFile.dir, append = TRUE)  #write in log file.
    return(ChIAPET)
}
#-------------
#-------------
####################### Converting to GInteractions
#-------------
#-------------
# trim anchors, dist and create GInteractions object:
GInteractionsCovnert_fun = function(ChIAPET, ArgPairedData) {
    #------------
    # break returned in anchors:
    #------------
    if (!is.null(ArgPairedData$LogFile.dir)) {
        # then the function called by PETClassification from tags
        Anchor1 = GenomicAlignments::first(x = ChIAPET, real.strand = FALSE)
        Anchor1 = methods::as(Anchor1, "GRanges")
        # to tags
        Anchor2 = GenomicAlignments::last(x = ChIAPET, real.strand = FALSE)
        Anchor2 = methods::as(Anchor2, "GRanges")
    } else {
        # then the function called by ConvertToPSelf: from tags
        Anchor1 = InteractionSet::anchors(ChIAPET, type = "first")
        # to tags
        Anchor2 = InteractionSet::anchors(ChIAPET, type = "second")
    }
    #--------------------------------------------
    #-------Remove black listed regions:
    #--------------------------------------------
    if (!is.null(ArgPairedData$BlackList)) {
        BlackListRM = BlackListCorrection_fun(Anchor1 = Anchor1, Anchor2 = Anchor2, 
            BlackList = ArgPairedData$BlackList, LogFile.dir = ArgPairedData$LogFile.dir)
    } else {
        BlackListRM = c()
    }
    #------------
    # remove false mapped PETs, outside of genomic regions, and unspecified strand
    # anchors.
    #------------
    TrimmedRM = TrimAnchors_fn(Anchor1 = Anchor1, Anchor2 = Anchor2, ArgPairedData = ArgPairedData)
    # total removals:
    TotalRM = unique(c(BlackListRM, TrimmedRM))
    # returns error if all are emptied
    if (length(TotalRM) != 0) {
        Anchor1 = Anchor1[-TotalRM]
        Anchor2 = Anchor2[-TotalRM]
    }
    if (length(Anchor1) == 0) {
        stop("All PETs ar removed from the data, exiting.", call. = FALSE)
    }
    #------------
    # make object:
    #------------
    cat("Converting to GInteractions class...\n")
    PairedData = InteractionSet::GInteractions(anchor1 = Anchor1, anchor2 = Anchor2)
    #------------
    # reduce PCR:
    #------------
    Nbefore = length(PairedData)  #before PCR removal
    cat("Reducing noise from PCR amplification procedures...")
    PairedData = S4Vectors::unique(PairedData)
    Nafter = length(PairedData)
    LogFile1 = paste("Total PCR replicates removed: ", Nbefore - Nafter, "\n")
    cat(LogFile1)
    LogFile2 = paste("Total PETs left: ", Nafter, "\n")
    cat(LogFile2)
    if (!is.null(ArgPairedData$LogFile.dir)) {
        write(LogFile1, file = ArgPairedData$LogFile.dir, append = TRUE)  #write in log file.
        write(LogFile2, file = ArgPairedData$LogFile.dir, append = TRUE)  #write in log file.
    }
    #------------
    # add DataSeqinfo in the data:
    #------------
    ChromInData = GenomeInfoDb::seqnames(GenomicRanges::seqinfo(PairedData))
    Matchsl = match(ChromInData, ArgPairedData$ChromLengths$Chrom)  #match names
    DataSeqinfo = GenomeInfoDb::Seqinfo(seqnames = ChromInData, seqlengths = ArgPairedData$ChromLengths$size[Matchsl], 
        genome = ArgPairedData$GenInfo$provider_version)
    GenomicRanges::seqinfo(PairedData) = DataSeqinfo
    #------------
    # Find spans/Dist:
    #------------
    PairedData$Dist = InteractionSet::pairdist(PairedData, type = "span")
    return(PairedData)
}  #done
#-------------
#-------------
# function for removing black listed PETs
BlackListCorrection_fun = function(Anchor1, Anchor2, BlackList, LogFile.dir) {
    cat("Removing black-listed PETs...")
    #---------------
    # make Granges for the black list:
    #---------------
    BlackList = GenomicRanges::GRanges(seqnames = BlackList$Chrom, ranges = IRanges::IRanges(start = BlackList$Region.Start, 
        end = BlackList$Region.End))
    #---------------
    #-----Mark First achor:
    #---------------
    cat("(checking first Anchor...)")
    Anchor1bl = InteractionSet::findOverlaps(Anchor1, BlackList)
    Anchor1bl = S4Vectors::queryHits(Anchor1bl)
    #---------------
    #-----Mark second achor:
    #---------------
    cat("(checking second Anchor...)\n")
    Anchor2bl = InteractionSet::findOverlaps(Anchor2, BlackList)
    Anchor2bl = S4Vectors::queryHits(Anchor2bl)
    #---------------
    #-----Merge and reduce:
    #---------------
    BlackListRM = unique(c(Anchor1bl, Anchor2bl))
    LogFile1 = paste("Total black-listed PETs removed:", length(BlackListRM), "\n")
    cat(LogFile1)
    LogFile2 = paste("Total PETs left:", length(Anchor1) - length(BlackListRM), "\n")
    cat(LogFile2)
    if (!is.null(LogFile.dir)) {
        write(LogFile1, file = LogFile.dir, append = TRUE)  #write in log file.
        write(LogFile2, file = LogFile.dir, append = TRUE)  #write in log file.
    }
    return(BlackListRM)
}  #done
#-------------
#-------------
# remove false mapped PETs, outside of genomic regions:
TrimAnchors_fn = function(Anchor1, Anchor2, ArgPairedData) {
    cat("Checking if any PETs have to be removed...")
    #--------------------------------------------
    #----------trim on coordinates outside borders
    #--------------------------------------------
    #------------
    # trim Anchor1:
    #------------
    StartAnchor1 = which(GenomicRanges::start(Anchor1) < 1)
    SeqPlace1 = match(as.character(GenomicRanges::seqnames(Anchor1)), ArgPairedData$ChromLengths$Chrom)
    EndAnchor1 = which(GenomicRanges::end(Anchor1) > ArgPairedData$ChromLengths$size[SeqPlace1])
    StrandAnchor1 = which(as.character(GenomicRanges::strand(Anchor1)) == "*")
    #------------
    # trim Anchor2:
    #------------
    StartAnchor2 = which(GenomicRanges::start(Anchor2) < 1)
    SeqPlace2 = match(as.character(GenomicRanges::seqnames(Anchor2)), ArgPairedData$ChromLengths$Chrom)
    EndAnchor2 = which(GenomicRanges::end(Anchor2) > ArgPairedData$ChromLengths$size[SeqPlace2])
    StrandAnchor2 = which(as.character(GenomicRanges::strand(Anchor2)) == "*")
    #------------
    # remove PETs:
    #------------
    RemAnchors = unique(c(StartAnchor1, EndAnchor1, StrandAnchor1, StartAnchor2, 
        EndAnchor2, StrandAnchor2))
    if (length(RemAnchors) != 0) {
        LogFile1 = paste(length(RemAnchors), " PETs were removed.\n")
        cat(LogFile1)
    } else {
        LogFile1 = "No PETs needed to be removed.\n"
        cat(LogFile1)
    }
    if (!is.null(ArgPairedData$LogFile.dir)) {
        write(LogFile1, file = ArgPairedData$LogFile.dir, append = TRUE)  #write in log file.
    }
    return(RemAnchors)
}  #done
#-------------
#-------------
####################### Function for Inter PETs
#-------------
#-------------
FindInter_fun = function(x, ArgFindInter) {
    cat("Separating Inter-chromosomal data...\n")
    #------------
    # keep inter pets:
    #------------
    interpets = subset(x, is.na(x$Dist))
    #------------
    if (length(interpets) != 0) {
        #------------
        # remove distance
        #------------
        interpets$Dist = NULL
        #------------
        # make them unique:
        #------------
        interpets = S4Vectors::unique(interpets)
        #------------
        # reduce the DataSeqinfo based on the chromosomes inside the data:
        #------------
        Anchor1 = InteractionSet::anchors(interpets, type = "first")
        Anchor2 = InteractionSet::anchors(interpets, type = "second")
        Anchor1LevelsUsed = GenomeInfoDb::seqlevelsInUse(Anchor1)
        Anchor2LevelsUsed = GenomeInfoDb::seqlevelsInUse(Anchor2)
        LevelsUsed = unique(c(Anchor1LevelsUsed, Anchor2LevelsUsed))
        LevelsUsed = gtools::mixedsort(LevelsUsed)
        #------------
        # update Anchors:
        #------------
        GenomeInfoDb::seqlevels(Anchor1) = LevelsUsed
        GenomeInfoDb::seqlevels(Anchor2) = LevelsUsed
        interpets = InteractionSet::GInteractions(Anchor1, Anchor2)
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
        #------------
        # save:
        #------------
        S4Vectors::metadata(interpets) = list(InteractionCounts = Inter.Counts)
        S4Vectors::metadata(interpets)$GenInfo = ArgFindInter$GenInfo
        class(interpets) = "PInter"
    } else {
        interpets = NULL  #empty
    }
    return(interpets)
}
#-------------
#-------------
####################### 
#-------------
#-------------
# Function For SelfIntra classification using the elbow method:
FindSelfIntra_fun = function(x, ArgFindSelfIntra) {
    # global variables for Rcheck:
    Freq = NULL
    Size = NULL
    PointDist = NULL
    logSize = NULL
    #--------------------------
    cat("Separating Self and Intra PETs...\n")
    #------------
    # Find Elbow Point:
    #------------
    MaxSpan = max(x$Dist)
    MinSpan = min(x$Dist)
    SeqSpan = seq(from = 0, to = MaxSpan, by = 100)
    SeqSpan[which(SeqSpan == max(SeqSpan))] = MaxSpan
    IntSpan = cbind(SeqSpan[-length(SeqSpan)], SeqSpan[-1])
    IntSpan = intervals::Intervals(IntSpan, closed = rep(TRUE, 2))
    DataInt = intervals::Intervals(cbind(x$Dist, x$Dist), closed = rep(TRUE, 2))
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
    SelfIndicator = which(x$Dist <= ElbowPoint$Size)
    MAXcut = max(x$Dist[SelfIndicator])
    MINcut = min(x$Dist[SelfIndicator])
    SelfBorder = data.frame(MAX = MAXcut, MIN = MINcut)
    LogFile1 = paste("Self-ligated cut-offs at: ", SelfBorder$MIN, "bp and ", SelfBorder$MAX, 
        "bp\n")
    cat(LogFile1)
    write(LogFile1, file = ArgFindSelfIntra$LogFile.dir, append = TRUE)  #write in log file.
    #------------
    # image:
    #------------
    if (ArgFindSelfIntra$PopImage) {
        # Plot the elbow:
        Self_ligated_borders_plot = ggplot2::ggplot(SpanDF, ggplot2::aes(x = Size, 
            y = Freq)) + # add histogram:
        ggplot2::geom_rect(ggplot2::aes(xmin = SpanDF$Size - 49, xmax = SpanDF$Size + 
            49, ymin = 0, ymax = SpanDF$Freq), size = 0.3, fill = "grey69", color = "black") + 
            # add line:
        ggplot2::geom_line(color = "blue", size = 0.6) + # add cut-off line:
        ggplot2::geom_vline(xintercept = ElbowPoint$Size, linetype = "dashed", color = "red", 
            size = 0.6) + # change scale:
        ggplot2::scale_x_log10(labels = function(co) round(log10(co))) + # titles
        ggplot2::ggtitle(paste("Elbow-point cut-off for Self/Intra Pets at: ", SelfBorder$MAX, 
            " bp")) + ggplot2::xlab("Sorted log-sizes of PET") + ggplot2::ylab("Frequency") + 
            ggplot2::theme_bw() + ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"), 
            panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), 
            panel.background = ggplot2::element_blank()) + ggplot2::theme(axis.text = ggplot2::element_text(size = 15, 
            color = "black"), axis.title = ggplot2::element_text(size = 18, color = "black"), 
            plot.title = ggplot2::element_text(size = 20, color = "black")) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        ggplot2::ggsave(plot = Self_ligated_borders_plot, file = file.path(ArgFindSelfIntra$AnalysisDir, 
            "Self_ligated_borders_plot.jpg"), scale = 2)
    }
    return(SelfIndicator)
}
#-------------
#-------------
####################### Functions for Intra PETs:
#-------------
#-------------
# Function for Intra PETs
FindIntra_fun = function(x, SelfIndicator, ArgFindIntra) {
    cat("Separating Intra-chromosomal data...\n")
    #------------
    # keep intra:
    #------------
    intrapets = x[-SelfIndicator, ]
    if (length(intrapets) != 0) {
        #------------
        # remove repeats:
        #------------
        intrapets = S4Vectors::unique(intrapets)
        #------------
        # keep to save
        #------------
        Mcols = S4Vectors::mcols(intrapets)
        #------------
        # reduce the DataSeqinfo based on the chromosomes inside the data:
        #------------
        Anchor1 = InteractionSet::anchors(intrapets, type = "first")
        Anchor2 = InteractionSet::anchors(intrapets, type = "second")
        Anchor1LevelsUsed = GenomeInfoDb::seqlevelsInUse(Anchor1)
        Anchor2LevelsUsed = GenomeInfoDb::seqlevelsInUse(Anchor2)
        LevelsUsed = unique(c(Anchor1LevelsUsed, Anchor2LevelsUsed))
        LevelsUsed = gtools::mixedsort(LevelsUsed)
        #------------
        # update Anchors:
        #------------
        GenomeInfoDb::seqlevels(Anchor1) = LevelsUsed
        GenomeInfoDb::seqlevels(Anchor2) = LevelsUsed
        intrapets = InteractionSet::GInteractions(Anchor1, Anchor2)
        S4Vectors::mcols(intrapets) = Mcols  #return
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
        S4Vectors::metadata(intrapets) = list(InteractionCounts = Intra.counts)
        S4Vectors::metadata(intrapets)$GenInfo = ArgFindIntra$GenInfo
        #------------
        # remove Dist variable:
        #------------
        intrapets$Dist = NULL
        #------------
        # change class
        #------------
        class(intrapets) = "PIntra"
    } else {
        intrapets = NULL  #empty.
    }
    return(intrapets)
}  #done
#-------------
#-------------
####################### Functions for Self PETs:
#-------------
#-------------
# Function for Self PETs
FindSelf_fun = function(x, SelfIndicator, ArgFindSelf, PSelf.Convert) {
    cat("Separating Self-ligated data...\n")
    #------------
    # keep self:
    #------------
    selfpets = x[SelfIndicator, ]
    Mcols = S4Vectors::mcols(selfpets)
    #------------
    # reduce the DataSeqinfo based on the chromosomes inside the data:
    #------------
    Anchor1 = InteractionSet::anchors(selfpets, type = "first")
    Anchor2 = InteractionSet::anchors(selfpets, type = "second")
    Anchor1LevelsUsed = GenomeInfoDb::seqlevelsInUse(Anchor1)
    Anchor2LevelsUsed = GenomeInfoDb::seqlevelsInUse(Anchor2)
    LevelsUsed = unique(c(Anchor1LevelsUsed, Anchor2LevelsUsed))
    LevelsUsed = gtools::mixedsort(LevelsUsed)
    # update Anchors:
    GenomeInfoDb::seqlevels(Anchor1) = LevelsUsed
    GenomeInfoDb::seqlevels(Anchor2) = LevelsUsed
    selfpets = InteractionSet::GInteractions(Anchor1, Anchor2)
    S4Vectors::mcols(selfpets) = Mcols  #return mcols infor
    #------------
    # find self-ligated length:
    #------------
    SLmean = round(mean(selfpets$Dist))
    LogFile1 = paste("Self-ligated mean length: ", SLmean, "\n")
    cat(LogFile1)
    if (!PSelf.Convert) {
        write(LogFile1, file = ArgFindSelf$LogFile.dir, append = TRUE)  #write in log file.
    }
    #------------
    # Take info:
    #------------
    Self_info = table(GenomicRanges::seqnames(InteractionSet::anchors(selfpets, type = "first")))
    Self_info = data.frame(Chrom = as.character(names(Self_info)), PET.counts = as.numeric(Self_info))
    S4Vectors::metadata(selfpets) = list(Self_info = Self_info)
    MaxSize = max(selfpets$Dist)
    MinSize = min(selfpets$Dist)
    selfpets$Dist = NULL  #remove distance
    S4Vectors::metadata(selfpets)$SLmean = SLmean
    S4Vectors::metadata(selfpets)$GenInfo = ArgFindSelf$GenInfo
    S4Vectors::metadata(selfpets)$MaxSize = MaxSize
    S4Vectors::metadata(selfpets)$MinSize = MinSize
    # update class:
    class(selfpets) = "PSelf"
    return(selfpets)
}  #done
#-------------
#-------------
####################### 
