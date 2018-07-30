#' @importFrom S4Vectors metadata DataFrame
#' @importFrom GenomeInfoDb seqinfo seqnames seqlevelsInUse seqlevels
#' @importFrom IRanges IRanges
#' @importFrom InteractionSet GInteractions anchors
#' @importFrom GenomicRanges GRanges start end seqnames findOverlaps reduce
#' @importFrom plyr dlply . alply ddply
#' @importFrom stats p.adjust smooth.spline dpois ppois quantile
#' @importFrom futile.logger flog.info
#' @importFrom methods is
#' @importFrom BiocParallel bpstart bpstop bpisup bpparam bpworkers bplapply
#' @import BH
#' @importFrom bigmemory filebacked.big.matrix describe attach.big.matrix as.matrix
#' @importClassesFrom bigmemory big.matrix
############################################# Main function for stage 4:
#--------------------------------------------
# Main function for interaction calls: GENERAL NOTE: Currently the inter are not
# supported.  If not intra is found either the algorithm will return error of
# empty data.
#--------------------------------------------
Stage_4_Main_fun = function(SA_prefix, S4_AnalysisDir, S4_FitObject, S4_IntraObject = NULL,
    S4_InterObject = NULL, S4_FDR_peak, S4_method, S4_image, S4_minPETs, S4_PeakExt) {
    #---------------
    # NOTE: inter not included at the time so set them to NULL:
    #---------------
    S4_InterObject = NULL
    # global variables for Rcheck: Take time:
    Analysis.time.start = Sys.time()
    # create folder to save data:
    if (!dir.exists(S4_AnalysisDir))
        dir.create(S4_AnalysisDir)
    #----------------------------
    # Write in file:
    #----------------------------
    futile.logger::flog.info(paste("Minimum number of allowed interaction PETs is set to:",
        S4_minPETs), name = "SA_LogFile", capture = FALSE)
    futile.logger::flog.info(paste("FDR cut-off of peaks to be used in the analysis: ",
        S4_FDR_peak), name = "SA_LogFile", capture = FALSE)
    #--------------------------------------------
    #------Prepare the data for the interactions:
    #--------------------------------------------
    cat("Preparing interactions data...\n")
    InteractionData = Prepare_InteractionsData_fun(SA_prefix = SA_prefix, S4_AnalysisDir = S4_AnalysisDir,
        S4_FitObject = S4_FitObject, S4_IntraObject = S4_IntraObject, S4_InterObject = S4_InterObject,
        S4_FDR_peak = S4_FDR_peak, S4_image = S4_image, S4_minPETs = S4_minPETs,
        S4_PeakExt = S4_PeakExt)
    #--------------------------------------------
    #------Run interaction analysis:
    #--------------------------------------------
    cat("|=================== Running interactions analysis ===================|\n")
    InteractionResults = Run_InteractionAnalysis_fun(InteractionData = InteractionData,
        S4_method = S4_method)
    futile.logger::flog.info("=====================================", name = "SA_LogFile",
        capture = FALSE)
    futile.logger::flog.info(paste("Total interactions processed: ", InteractionResults$TotIntAdded),
        name = "SA_LogFile", capture = FALSE)
    futile.logger::flog.info(paste("Total bi-products removed: ", InteractionResults$TotBiRem),
        name = "SA_LogFile", capture = FALSE)
    #--------------------------------------------
    #----Create Genome Map object (suppressWarnings for trimming)
    #--------------------------------------------
    suppressWarnings(Create_GenomeMapObject(InteractionResults = InteractionResults,
        SA_prefix = SA_prefix, S4_AnalysisDir = S4_AnalysisDir))
    futile.logger::flog.info("The Genome map is successfully built!", name = "SA_LogFile",
        capture = FALSE)
    #--------------------------------------------
    #------print and write the end:
    #--------------------------------------------
    futile.logger::flog.info("=====================================", name = "SA_LogFile",
        capture = FALSE)
    futile.logger::flog.info("Stage 4 is done!", name = "SA_LogFile", capture = FALSE)
    futile.logger::flog.info(paste("Analysis results for stage 4 are in:\n", S4_AnalysisDir),
        name = "SA_LogFile", capture = FALSE)
    # save time:
    Analysis.time.end = Sys.time()
    Total.Time = Analysis.time.end - Analysis.time.start
    LogFile = paste("Total stage 4 time:", Total.Time, " ", units(Total.Time))
    futile.logger::flog.info(LogFile, name = "SA_LogFile", capture = FALSE)
}
# done
#----------------------------
# Main function for preparing the data
#----------------------------
Prepare_InteractionsData_fun = function(SA_prefix, S4_AnalysisDir, S4_FitObject,
    S4_IntraObject, S4_InterObject, S4_FDR_peak, S4_image, S4_minPETs, S4_PeakExt) {
    #----------------------------
    # Prepare the peaks data:
    # suppress trimming warnings which I dont need
    #----------------------------
    PeaksSplit = suppressWarnings(Get_PeaksSplit_fun(S4_FitObject = S4_FitObject, S4_FDR_peak = S4_FDR_peak,
        S4_PeakExt = S4_PeakExt))
    #----------------------------
    # Get the PETs split data
    #----------------------------
    PETsSplit = Get_PETsSplit_fun(S4_IntraObject = S4_IntraObject, S4_InterObject = S4_InterObject)
    #----------------------------
    # Find overlaps of PETs with the peaks
    #----------------------------
    ConnectionPETsData = Create_ConnectionPETsData_fun(PETsData = PETsSplit$PETsData,
        PeaksGR = PeaksSplit$PeaksGR)
    #----------------------------
    # Create the structures used in the analysis. BigMatrices, Networks etc.
    #----------------------------
    InteractionData = Create_InteractionStructures_fun(SA_prefix = SA_prefix, S4_AnalysisDir = S4_AnalysisDir,
        ConnectionPETsData = ConnectionPETsData, PeaksSplit = PeaksSplit, PETsSplit = PETsSplit,
        S4_image = S4_image, S4_minPETs = S4_minPETs)
    return(InteractionData)
}
# done
#----------------------------
# Function for preparing the Peaks
#----------------------------
Get_PeaksSplit_fun = function(S4_FitObject, S4_FDR_peak, S4_PeakExt) {
    # Rcheck:
    FDR=NULL
    queryHits=NULL
    #------
    # take Peak data and subset by significance:
    #------
    PeakData = S4Vectors::metadata(S4_FitObject)$Peaks.Info
    NGlobalPeaks = nrow(PeakData)  #for print, and image
    PeakData = subset(PeakData, FDR < S4_FDR_peak)  #subset the data
    SeqInfo = GenomeInfoDb::seqinfo(S4_FitObject)  #to save again
    # keep what you need for the analysis:
    PeakData = PeakData[, c("Chrom", "Peak.Summit", "CIQ.Up.start", "CIQ.Down.end", "FDR")]
    NFDRPeaks = nrow(PeakData)  #the total number of Peaks left in the data after FDR
    #------
    # write in file and print:
    #------
    if (NFDRPeaks == 0)
        stop("No peaks to use for interaction analysis on FDR ", S4_FDR_peak, "! Use a higher threshold!",
            call. = FALSE)
    futile.logger::flog.info(paste("Total peaks passing the FDR cut-off:", NFDRPeaks,
        "(", NFDRPeaks/NGlobalPeaks * 100, "% of the total peaks )"), name = "SA_LogFile",
        capture = FALSE)
    #------
    # Extend the peaks from their center:
    #------
    futile.logger::flog.info(paste("Extending peak intervals by", S4_PeakExt, "bp on either side."),
        name = "SA_LogFile", capture = FALSE)
    LeftInterval = PeakData$CIQ.Up.start - S4_PeakExt
    RightInterval = PeakData$CIQ.Down.end + S4_PeakExt
    #------
    # Merge overlapping and create the final PeaksGR object
    #------
    futile.logger::flog.info("Merging overlapping peaks...Done", name = "SA_LogFile",
        capture = FALSE)
    NM_PeaksIR = IRanges::IRanges(start = LeftInterval, end = RightInterval)
    NM_PeaksGR = GenomicRanges::GRanges(seqnames = PeakData$Chrom, ranges = NM_PeaksIR,
        seqinfo = SeqInfo) #not merged yet
    PeaksGR = GenomicRanges::reduce(x = NM_PeaksGR, min.gapwidth = 0, ignore.strand = TRUE) #merged
    NPeaksMerged = length(PeaksGR) #total peaks after merging
    # Find overlaps from PeaksGR to NM_PeaksGR to decide the summit based on FDR
    MergedNotMergedOv = GenomicRanges::findOverlaps(query = PeaksGR, subject = NM_PeaksGR,
        ignore.strand = TRUE)
    MergedNotMergedOv = as.data.frame(MergedNotMergedOv)
    MergedNotMergedOv = MergedNotMergedOv[with(MergedNotMergedOv, order(queryHits)),]
    # find the summit:
    PeakSummit = Get_NewPeakSummit_fun_Rcpp(MergedNotMergedOv$queryHits,
        MergedNotMergedOv$subjectHits, PeakData$Peak.Summit, PeakData$FDR,
        nrow(MergedNotMergedOv), NPeaksMerged)
    # save:
    PeaksGR$PeakSummit = PeakSummit
    PeaksGR$LID = seq_len(NPeaksMerged)
    futile.logger::flog.info(paste("Total peaks after merging: ", NPeaksMerged),
        name = "SA_LogFile", capture = FALSE)
    return(list(PeaksGR = PeaksGR, NFDRPeaks = NFDRPeaks, NGlobalPeaks = NGlobalPeaks,
        NPeaksMerged = NPeaksMerged))
}
# done
#----------------------------
#--Function for the intra/inter ligated pets
#----------------------------
Get_PETsSplit_fun = function(S4_IntraObject, S4_InterObject) {
    #----------------------------
    # Merge the Intra-inter into one data: this will be used to check overlaps: But
    # first take cases for what exists:
    #----------------------------
    if (!is.null(S4_IntraObject) & !is.null(S4_InterObject)) {
        cat("Merging Intra- and Inter-chromosomal PETs data...")
        Anchor_1_intra = InteractionSet::anchors(S4_IntraObject, type = "first")
        Anchor_2_intra = InteractionSet::anchors(S4_IntraObject, type = "second")
        Anchor_1_inter = InteractionSet::anchors(S4_InterObject, type = "first")
        Anchor_2_inter = InteractionSet::anchors(S4_InterObject, type = "second")
        Anchor_1_both = c(Anchor_1_intra, Anchor_1_inter)
        Anchor_2_both = c(Anchor_2_intra, Anchor_2_inter)
        PETsData = InteractionSet::GInteractions(anchor1 = Anchor_1_both, anchor2 = Anchor_2_both)
    } else if (!is.null(S4_IntraObject)) {
        cat("Including only Intra-ligated PETs in the analysis (Inter-ligated are empty)...")
        Anchor_1_intra = InteractionSet::anchors(S4_IntraObject, type = "first")
        Anchor_2_intra = InteractionSet::anchors(S4_IntraObject, type = "second")
        PETsData = InteractionSet::GInteractions(anchor1 = Anchor_1_intra, anchor2 = Anchor_2_intra)
    } else {
        cat("Including only Inter-ligated PETs in the analysis (Intra-ligated are empty)...")
        Anchor_1_inter = InteractionSet::anchors(S4_InterObject, type = "first")
        Anchor_2_inter = InteractionSet::anchors(S4_InterObject, type = "second")
        PETsData = InteractionSet::GInteractions(anchor1 = Anchor_1_inter, anchor2 = Anchor_2_inter)
    }
    NGlobalInterPETs = length(PETsData)  #the total connection PETs in data
    if (NGlobalInterPETs == 0)
        stop("No Intra/Inter-chromosomal PETs to use for interaction analysis!",
            call. = FALSE)
    cat("Done\n")
    return(list(PETsData = PETsData, NGlobalInterPETs = NGlobalInterPETs))
}
# done
#----------------------------
# function creating the connection data: This function uses the Anchor 1 and 2 of
# each interaction PET in PETsData And it seperately finds overlaps with each
# Peak in PeaksGR . The resulting data will then be used to classify each hit/Tag
# in interactions PETs by the Get_PETsInfoMat_fun_Rcpp c++ funtion
#----------------------------
Create_ConnectionPETsData_fun = function(PETsData, PeaksGR) {
    # Rcheck:
    query = NULL
    # -------
    cat("Intersecting the PETs with the peaks...")
    #----------------------------
    # Get single overlap for Anchor 1 in PETsData with PeaksGR
    #----------------------------
    Anchor1 = InteractionSet::anchors(x = PETsData, type = "first")
    Anchor1ovlp = GenomicRanges::findOverlaps(query = Anchor1, subject = PeaksGR,
        ignore.strand = TRUE, type = "any", select = "all")
    Anchor1ovlp = as(Anchor1ovlp, "DataFrame")
    if (nrow(Anchor1ovlp) == 0)
        stop("None of the intra/inter-ligation PETs overlaps with any Peak!", call. = FALSE)
    names(Anchor1ovlp) = c("query", "subject")
    Anchor1ovlp$Type = 1  #means anchor 1.
    # Find the tag to be used in the classification process
    Anchor1ovlp$Tag = (GenomicRanges::start(Anchor1[Anchor1ovlp$query]) + GenomicRanges::end(Anchor1[Anchor1ovlp$query]))/2
    #----------------------------
    # Get single overlap for Anchor 2 with PeaksGR
    #----------------------------
    Anchor2 = InteractionSet::anchors(x = PETsData, type = "second")
    Anchor2ovlp = GenomicRanges::findOverlaps(query = Anchor2, subject = PeaksGR,
        ignore.strand = TRUE, type = "any", select = "all")
    Anchor2ovlp = as(Anchor2ovlp, "DataFrame")
    if (nrow(Anchor2ovlp) == 0)
        stop("None of the intra/inter-ligation PETs overlaps with any Peak!", call. = FALSE)
    names(Anchor2ovlp) = c("query", "subject")
    Anchor2ovlp$Type = 2  #means anchor 2
    # Find the tag to be used in the classification process
    Anchor2ovlp$Tag = (GenomicRanges::start(Anchor2[Anchor2ovlp$query]) + GenomicRanges::end(Anchor2[Anchor2ovlp$query]))/2
    #----------------------------
    # Merge the Anchors and sort by query name (Tag): Then it is easier to subset by
    # query in c++
    #----------------------------
    ConPETsData = rbind(Anchor1ovlp, Anchor2ovlp)
    NIntTagsloop = nrow(ConPETsData)  #the total number of tags to classify in c++
    OrderQuery = order(ConPETsData$query)
    ConPETsData = ConPETsData[OrderQuery, ]
    #----------------------------
    # Give information that you need to run the classification in c++:
    #----------------------------
    ConPETsData$LID = PeaksGR$LID[ConPETsData$subject]
    ConPETsData$PeakSummit = PeaksGR$PeakSummit[ConPETsData$subject]
    cat("Done\n")
    return(list(ConPETsData = ConPETsData, NIntTagsloop = NIntTagsloop))
}
# done
#----------------------------
# Main function for creating the interaction structures used in the analysis
#----------------------------
Create_InteractionStructures_fun = function(SA_prefix, S4_AnalysisDir, ConnectionPETsData,
    PeaksSplit, PETsSplit, S4_image, S4_minPETs) {
    # Rcheck:
    V1 = NULL
    V2 = NULL
    #----------------------------
    # Create the PETsInfoMat:
    #----------------------------
    cat("Counting PETs in the peaks...")
    PETsInfoMat = Get_PETsInfoMat_fun_Rcpp(ConnectionPETsData$ConPETsData$query,
        ConnectionPETsData$ConPETsData$Type, ConnectionPETsData$ConPETsData$Tag,
        ConnectionPETsData$ConPETsData$LID, ConnectionPETsData$ConPETsData$PeakSummit,
        PETsSplit$NGlobalInterPETs, ConnectionPETsData$NIntTagsloop)
    if (is.null(PETsInfoMat))
        stop("None of the intra/inter-ligation PETs overlaps with two Peaks!", call. = FALSE)
    #----------------------------
    # Turn PETsInfoMat into a data.frame and split by V1 (PBS-i),and V2 (PBS-j) to
    # get statistics. At the end of the following the InteractionInfo Will have
    # col-1, PBS-i, col2-PBSj and col3 total PETs. The PBS ids are in R index
    #----------------------------
    PETsInfoMat = as.data.frame(PETsInfoMat)
    InteractionInfo = plyr::dlply(.data = PETsInfoMat, plyr::.(V1, V2), .fun = Get_PETsStats_fun,
        S4_minPETs = S4_minPETs)
    InteractionInfo = do.call(rbind, InteractionInfo)
    NInteractions = nrow(InteractionInfo)  # the potential interactions
    if (NInteractions == 0)
        stop("No candidate interactions found at threshold S4_minPETs = ", S4_minPETs,
            ". Try reducing S4_minPETs.", call. = FALSE)
    rownames(InteractionInfo) = NULL
    NPETsInvolved = sum(InteractionInfo[, c(3)])  # total number of PETs involved
    #----------------------------
    # Find the chromosome combinations and return the data.
    #----------------------------
    ChromData = as.character(GenomicRanges::seqnames(PeaksSplit$PeaksGR))
    ChromCombData = Get_ChromCombData_fun(InteractionInfo = InteractionInfo, ChromData = ChromData)
    InteractionInfo = ChromCombData$InteractionInfo  #used for building the InteractionsMatrix
    UniqueChromComb = ChromCombData$UniqueChromComb  #Info about the Chromosome combinations, not used after next line
    NPeaksInvolved = sum(UniqueChromComb$TotPeaks)  #Peaks involved in general
    NetworkBuildData = ChromCombData$NetworkBuildData  #used for building the genomic networks
    cat("Done\n")
    #----------------------------
    # Summarize
    #----------------------------
    futile.logger::flog.info(paste("Total", NInteractions, "candidate interactions will be processed"),
        name = "SA_LogFile", capture = FALSE)
    futile.logger::flog.info(paste("Total", NPeaksInvolved, "peaks are involved in potential interactions",
        "(", NPeaksInvolved/PeaksSplit$NPeaksMerged * 100, "% of the total FDR peaks )"),
        name = "SA_LogFile", capture = FALSE)
    futile.logger::flog.info(paste("Total", NPETsInvolved, "PETs are involved in potential interactions",
        "(", NPETsInvolved/PETsSplit$NGlobalInterPETs * 100, "% of the total interaction PETs )"),
        name = "SA_LogFile", capture = FALSE)
    #----------------------------
    # create the InteractionInfMat:
    #----------------------------
    # rows=NInteractions, each row corresponds to one interaction col 1/2: the left
    # LID on the PeaksData (R indeces) col 3/4 are the ids on the network element (R
    # indeces for node names, changing) those are used for bi-products and for
    # finding paths when peaks are merged.  col 5/6 are the node ids for the Adj
    # matrices, not changing, col 7/8 are node IDS for the NiNjMat, not changing. The
    # NiNjMat is global in case you also include inter afterwards.  col 9/10 p-value
    # and FDR, col 11, the nij of the interaction col 12, the QCell ID in the
    # expected matrix for the DVijij(NA on inter) col 13:
    # order of the interaction col 14 is the Chromosome combination ID of the
    # interaction col 15 is the intra indicator(intra=1, inter=0)
    cat("Summarizing interaction information...")
    InteractionInfMat = matrix(0, nrow = NInteractions, ncol = 15)
    colnames(InteractionInfMat) = c("PBS_i", "PBS_j", "AdjBiNode_i", "AdjBiNode_j",
        "AdjNode_i", "AdjNode_j", "NiNjNode_i", "NiNjNode_j", "pvalue", "FDR", "nij",
        "QCell", "Order", "Chrom12ID", "IntraID")
    #----------------------------
    # fill up InteractionInfMat, dont need to return it, it is by reference Also
    # count the observed PETs in each combination.
    #----------------------------
    Indeces_Vij = Initiate_InteractionInfMat_fun_Rcpp(InteractionInfMat, InteractionInfo,
        NPeaksInvolved, NInteractions)
    cat("Done\n")
    #----------------------------
    # plot:
    #----------------------------
    if (S4_image) {
        Get_image_S4_P1_fun(S4_AnalysisDir = S4_AnalysisDir, SA_prefix = SA_prefix,
            NGlobalInterPETs = PETsSplit$NGlobalInterPETs, NPETsInvolved = NPETsInvolved,
            NFDRPeaks = PeaksSplit$NFDRPeaks, NGlobalPeaks = PeaksSplit$NGlobalPeaks,
            NPeaksInvolved = NPeaksInvolved)
    }
    #----------------------------
    # Create the Newtorks as well as the Adjucency matrices
    #----------------------------
    NetworksData = Create_Networks_fun(NetworkBuildData = NetworkBuildData, PeakSummit = PeaksSplit$PeaksGR$PeakSummit,
        SA_prefix = SA_prefix, S4_AnalysisDir = S4_AnalysisDir)
    #----------------------------
    # return:
    #----------------------------
    return(list(PeaksGR = PeaksSplit$PeaksGR, InteractionInfMat = InteractionInfMat,
        NInteractions = NInteractions, AllInteIndeces = Indeces_Vij$AllInteIndeces,
        NiNjMat = Indeces_Vij$NiNjMat, NetworksData = NetworksData))
}
# done
#----------------------------
# Function for getting the sampled statistics for each PET group
#----------------------------
Get_PETsStats_fun = function(x, S4_minPETs) {
    #----------------------------
    # check if Total PETs are above the S4_minPETs:
    #----------------------------
    N_x = nrow(x)  #Pets in the interaction
    if (N_x < S4_minPETs)
        return(NULL)
    #----------------------------
    # Else return the total PETs for each interaction:
    #----------------------------
    InteractionInfo_x = c(x$V1[1], x$V2[1], N_x)
    return(InteractionInfo_x)
}
# done
#----------------------------
# Function to give new node names to peaks and separate them by chromosome
# combination Important: This functions works for intra. It will not crash if
# Intra are not empty.  If intra are empty, the algorithm will return an error at
# a previous point.  If you include inter at any time, this function has to be
# updated
#----------------------------
Get_ChromCombData_fun = function(InteractionInfo, ChromData) {
    #----------------------------
    # First Create the Chromosome Combination data:
    #----------------------------
    ChromComb = data.frame(Chrom1 = ChromData[InteractionInfo[, c(1)]], Chrom2 = ChromData[InteractionInfo[,
        c(2)]])
    ChromComb$Chrom12 = paste(ChromComb$Chrom1, "-", ChromComb$Chrom2, sep = "")  #merged comb
    # Get unique Comb, make data frame and give unique IDS:
    UniqueChromComb = table(ChromComb$Chrom12)
    UniqueChromComb = as.data.frame(UniqueChromComb, stringsAsFactors = FALSE)
    names(UniqueChromComb) = c("Chrom12", "Totals")  #note the Totals are total interactions, NOT Peaks
    # Check which intra/inter:
    Chrom12 = do.call(rbind, strsplit(x = UniqueChromComb$Chrom12, split = "-"))
    IntraComb = which(Chrom12[, c(1)] == Chrom12[, c(2)])
    UniqueChromComb$Intra = 0  #indicating no intra combination
    if (length(IntraComb) != 0)
        UniqueChromComb$Intra[IntraComb] = 1
    # sort intra/inter combination before giving IDS: This helps to create a list
    # with no empty entries afterwards As inter=0 will go at the end
    UniqueChromComb = UniqueChromComb[with(UniqueChromComb, order(Intra, decreasing = TRUE)),
        ]
    UniqueChromComb$Chrom12ID = seq_len(nrow(UniqueChromComb))  #give unique IDS
    #----------------------------
    # Assign the combination ID to the InteractionInfo
    #----------------------------
    MatchChrom12 = match(ChromComb$Chrom12, UniqueChromComb$Chrom12)
    Intra = UniqueChromComb$Intra[MatchChrom12]
    Chrom12ID = UniqueChromComb$Chrom12ID[MatchChrom12]
    #----------------------------
    # Make the InteractionInfo to data.frame. Assign the Chrom12ID and Intra Split
    # the data frame to also given Node names for the big.matrices
    #----------------------------
    InteractionInfo = as.data.frame(InteractionInfo)
    InteractionInfo$Chrom12ID = Chrom12ID
    InteractionInfo$Intra = Intra
    colnames(InteractionInfo) = c("PBS_i", "PBS_j", "nij", "Chrom12ID", "Intra")
    # Assign the new node names, this will also assign node names to inter, but it is
    # ok:
    InteractionInfo = plyr::ddply(InteractionInfo, plyr::.(Chrom12ID), function(x) {
        # take the PBS:
        PBS_x = sort(unique(c(x$PBS_i, x$PBS_j)))
        # Match the PBS_i/j to the PBS_x and this will give the new node IDS:
        x$AdjNode_i = match(x$PBS_i, PBS_x)
        x$AdjNode_j = match(x$PBS_j, PBS_x)
        # return:
        return(x)
    })
    # Also take the Total Peaks involved in each combination for creating the
    # matrices:
    TotPeaksComb = plyr::ddply(InteractionInfo, plyr::.(Chrom12ID), function(x) {
        # take the PBS:
        Tot_x = length(unique(c(x$PBS_i, x$PBS_j)))
        # return:
        return(data.frame(Chrom12ID = x$Chrom12ID[1], Tot_x = Tot_x))
    })
    UniqueChromComb$TotPeaks = TotPeaksComb$Tot_x[match(UniqueChromComb$Chrom12ID,
        TotPeaksComb$Chrom12ID)]
    #----------------------------
    # Then, create the NiNjNodes for the NiNj vector which is global:
    #----------------------------
    PBSunique = sort(unique(c(InteractionInfo$PBS_i, InteractionInfo$PBS_j)))
    InteractionInfo$NiNjNode_i = match(InteractionInfo$PBS_i, PBSunique)
    InteractionInfo$NiNjNode_j = match(InteractionInfo$PBS_j, PBSunique)
    #----------------------------
    # Also load the PBS and AdjNodes indeces to create List used for building the
    # networks This will only return intra-ligated networks The Chrom12ID are sorted,
    # so the list will have in Chrom12ID_comb position the Chrom12ID_comb of the
    # specific combination. So access is used by Chrom12ID_comb
    #----------------------------
    NetworkBuildData = plyr::dlply(InteractionInfo, plyr::.(Chrom12ID), function(x) {
        # If inter return NULL:
        if (x$Intra[1] == 0)
            return(NULL)
        # else take info
        PBS_ij = unique(c(x$PBS_i, x$PBS_j))
        AdjNode_ij = unique(c(x$AdjNode_i, x$AdjNode_j))
        NiNjNode_ij = unique(c(x$NiNjNode_i, x$NiNjNode_j))
        # take order to sort all the correct way:
        Order_PBS_ij = order(PBS_ij, decreasing = FALSE)
        PBS_comb = PBS_ij[Order_PBS_ij]
        AdjNode_comb = AdjNode_ij[Order_PBS_ij]
        NiNjNode_comb = NiNjNode_ij[Order_PBS_ij]
        Chrom12ID_comb = x$Chrom12ID[1]
        return(list(PBS_comb = PBS_comb, AdjNode_comb = AdjNode_comb, NiNjNode_comb = NiNjNode_comb,
            Chrom12ID_comb = Chrom12ID_comb))
    })
    # remove the NULL if any:
    RmNULL = unlist(lapply(NetworkBuildData, length))
    NetworkBuildData = NetworkBuildData[which(RmNULL != 0)]
    #----------------------------
    # Transform the InteractionInfo to matrix and return:
    #----------------------------
    InteractionInfo = as.matrix(InteractionInfo)
    return(list(InteractionInfo = InteractionInfo, UniqueChromComb = UniqueChromComb,
        NetworkBuildData = NetworkBuildData))
}
# done
#----------------------------
# function for plotting for stage 4 part 1: Used peaks/involved Peaks and used
# PETs
#----------------------------
Get_image_S4_P1_fun = function(S4_AnalysisDir, SA_prefix, NGlobalInterPETs, NPETsInvolved,
    NFDRPeaks, NGlobalPeaks, NPeaksInvolved) {
    # Rcheck:
    Value = NULL
    Kind = NULL
    # image dir:
    S4_P1_image_dir = file.path(S4_AnalysisDir, paste(SA_prefix, "_stage_4_p1_image.jpg",
        sep = ""))
    #-------------
    # create data:
    #-------------
    PeaksAboveFDR = (NGlobalPeaks - NFDRPeaks)/NGlobalPeaks * 100
    PeaksBelowFDRnotInv = (NFDRPeaks - NPeaksInvolved)/NGlobalPeaks * 100
    PeaksInvolved = NPeaksInvolved/NGlobalPeaks * 100
    PETSNotInvolved = (NGlobalInterPETs - NPETsInvolved)/NGlobalInterPETs * 100
    PETSInvolved = NPETsInvolved/NGlobalInterPETs * 100
    S4_imagedata_1 = data.frame(Kind = c(paste("Peaks not passing FDR (", round(PeaksAboveFDR,
        digits = 1), "%)", sep = ""), paste("Peaks passing FDR, but not involved in interactions (",
        round(PeaksBelowFDRnotInv, digits = 1), "%)", sep = ""), paste("Peaks passing FDR, and involved in interactions (",
        round(PeaksInvolved, digits = 1), "%)", sep = ""), paste("Interaction PETs not involved in interactions (",
        round(PETSNotInvolved, digits = 1), "%)", sep = ""), paste("Interaction PETs involved in interactions (",
        round(PETSInvolved, digits = 1), "%)", sep = "")), Value = c(PeaksAboveFDR,
        PeaksBelowFDRnotInv, PeaksInvolved, PETSNotInvolved, PETSInvolved), facet = c("Peaks Summary",
        "Peaks Summary", "Peaks Summary", "PETs summary", "PETs summary"))
    #-------------
    # plot the split:
    #-------------
    S4_image_p1 = ggplot2::ggplot(S4_imagedata_1, ggplot2::aes(x = "", y = Value,
        fill = factor(Kind))) + ggplot2::geom_bar(width = 1, stat = "identity") +
        ggplot2::facet_wrap(~facet) + ggplot2::coord_polar(theta = "y") + ggplot2::theme(axis.title = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 20, color = "black"), legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 17), axis.text = ggplot2::element_blank(),
        legend.position = "bottom", legend.direction = "vertical", axis.ticks = ggplot2::element_blank()) +
        ggplot2::ggtitle("Pie chart for the summary of the potential interactions.") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_fill_brewer(palette = "Dark2")
    # save:
    ggplot2::ggsave(plot = S4_image_p1, file = S4_P1_image_dir, scale = 2)
}
# done
#----------------------------
# Function for creating the Network structures and the Adj matrices
#----------------------------
Create_Networks_fun = function(NetworkBuildData, PeakSummit, SA_prefix, S4_AnalysisDir) {
    #----------------------------
    # Initiate the structures you will use:
    #----------------------------
    TotNetworks = length(NetworkBuildData)
    NetworkList = vector(mode = "list", length = TotNetworks)  #list of networks
    UpdateIndecesList = vector(mode = "list", length = TotNetworks)  #list of updating indeces for merging nodes
    MergedNodesIndexList = vector(mode = "list", length = TotNetworks)  #list of merging nodes IDs
    BigInfoMatDescList = vector(mode = "list", length = TotNetworks)  #list of description of the AdjMatrices
    NadjList = vector(mode = "list", length = TotNetworks)  #list of the length of the AdjMatrices
    NPeaksInvolvedList = vector(mode = "list", length = TotNetworks)  #list of the total Peaks in each combination
    InteractionPairsList = vector(mode = "list", length = TotNetworks)  #list of the total interaction pairs involved(reducing), note this does NOT include inter
    NiNjIndecesList = vector(mode = "list", length = TotNetworks)  #indeces in R, used to access the NiNjMat
    #----------------------------
    # create directory for the big.matrix:
    #----------------------------
    S4_MatrixDir = file.path(S4_AnalysisDir, "GenomeMapMatrices")
    if (dir.exists(S4_MatrixDir))
        unlink(S4_MatrixDir, recursive = TRUE)
    dir.create(S4_MatrixDir)
    #----------------------------
    # Loop in NetworkBuildData elements and create what you need:
    #----------------------------
    for (Net in seq_len(TotNetworks)) {
        #----------------------------
        # Take network elements:
        #----------------------------
        PBS_Net = NetworkBuildData[[Net]]$PBS_comb
        AdjNode_Net = NetworkBuildData[[Net]]$AdjNode_comb
        Chrom12ID_Net = NetworkBuildData[[Net]]$Chrom12ID_comb
        NPeaksInvolved_Net = length(AdjNode_Net)  #total peaks for network
        NiNjNode_Net = NetworkBuildData[[Net]]$NiNjNode_comb
        #----------------------------
        # Initiate Network structure:
        #----------------------------
        Network_Net = Initiate_GenomeMap_fun_Rcpp(NPeaksInvolved_Net, AdjNode_Net,
            PBS_Net, PeakSummit, Chrom12ID_Net)
        #----------------------------
        # Create AdjMatrix for the chromosome
        #----------------------------
        Nadj_Net = NPeaksInvolved_Net * (NPeaksInvolved_Net - 1)/2
        bkFile_BigInfoMat_Net = paste(SA_prefix, "_BigInfoMat_", Chrom12ID_Net, ".bk",
            sep = "")
        descFile_BigInfoMat_Net = paste(SA_prefix, "_BigInfoMat_", Chrom12ID_Net,
            ".desc", sep = "")
        BigInfoMat_Net = bigmemory::filebacked.big.matrix(nrow = Nadj_Net, ncol = 1,
            type = "double", backingfile = bkFile_BigInfoMat_Net, backingpath = S4_MatrixDir,
            descriptorfile = descFile_BigInfoMat_Net, dimnames = c(NULL, NULL))
        # get description:
        BigInfoMatDesc_Net = bigmemory::describe(BigInfoMat_Net)
        #----------------------------
        # create the indeces to run in parallel for AdjMat:
        #----------------------------
        PeakSeq_Net = seq_len(NPeaksInvolved_Net - 1)  #dont need the last one
        if (BiocParallel::bpisup(BiocParallel::bpparam())) {
            # then parallel backhead, so split in tasks get number of registered cores:
            Cores_Net = max(1, BiocParallel::bpworkers(BiocParallel::bpparam()))
            # get correct dimentions for the vector by extending with NAs
            length(PeakSeq_Net) = suppressWarnings(prod(dim(matrix(PeakSeq_Net, ncol = Cores_Net,
                byrow = TRUE))))
            # create the matrix:
            UpdateIndeces_Net = matrix(PeakSeq_Net, byrow = TRUE, ncol = Cores_Net)
            UpdateIndeces_Net = plyr::alply(.data = UpdateIndeces_Net, .margins = 2,
                .fun = function(x) x[!is.na(x)])
        } else {
            # no parallel backhead sto do it linearly
            UpdateIndeces_Net = list(PeakSeq_Net)
        }
        #----------------------------
        # Create also the MergeNodesIndex to keep track of the merging nodes
        #----------------------------
        MergedNodesIndex_Net = as.list(seq_len(NPeaksInvolved_Net))
        #----------------------------
        # Save element:
        #----------------------------
        NetworkList[[Net]] = Network_Net
        UpdateIndecesList[[Net]] = UpdateIndeces_Net
        MergedNodesIndexList[[Net]] = MergedNodesIndex_Net
        BigInfoMatDescList[[Net]] = BigInfoMatDesc_Net
        NadjList[[Net]] = Nadj_Net
        NPeaksInvolvedList[[Net]] = NPeaksInvolved_Net
        InteractionPairsList[[Net]] = Nadj_Net
        NiNjIndecesList[[Net]] = NiNjNode_Net
    }
    cat("|================ Network Initialization is finished =================|\n")
    # return:
    return(list(NetworkList = NetworkList, UpdateIndecesList = UpdateIndecesList,
        MergedNodesIndexList = MergedNodesIndexList, BigInfoMatDescList = BigInfoMatDescList,
        NadjList = NadjList, NPeaksInvolvedList = NPeaksInvolvedList, S4_MatrixDir = S4_MatrixDir,
        TotNetworks = TotNetworks, InteractionPairsList = InteractionPairsList, NiNjIndecesList = NiNjIndecesList))
}
# done
#----------------------------
# Main Interaction analysis function:
#----------------------------
Run_InteractionAnalysis_fun = function(InteractionData, S4_method) {
    #----------------------------
    # Initiate by breaking the InteractionData input:
    #----------------------------
    PeaksGR = InteractionData$PeaksGR  # peaks data for return
    InteractionInfMat = InteractionData$InteractionInfMat  #the matrix of information
    NInteractions = InteractionData$NInteractions  #total interactions to test
    AllInteIndeces = InteractionData$AllInteIndeces  #interaction ids
    NiNjMat = InteractionData$NiNjMat  #matrix with the ni/nj
    # network related data:
    NetworksData = InteractionData$NetworksData
    S4_MatrixDir = NetworksData$S4_MatrixDir
    NetworksData$S4_MatrixDir = NULL
    NetUpdateIndicator = seq_len(NetworksData$TotNetworks)  # Network Update indicators(for not updating all distances all the time)
    # Quantiles:
    QuantileProbs = seq(from = 0, to = 1, length = 1001)[-c(1)]  #the size of the bins(plus one will be the last bin with the inter)
    SavedQuantiles = lapply(as.list(seq_len(NetworksData$TotNetworks)), function(x) return(list(BinsVij = NA,
        ChunkSizeVij = NA)))  #save for not computing again
    # counters
    TotIntAdded = 0  #the total interactions added in the model
    TotBiRem = 0  #total bi-products removed from the model
    OrdersCount = 1  #Counter for the interactions added in the model(multiple interactions counter)
    #----------------------------
    # run the loop to add interactions:
    #----------------------------
    while (TotIntAdded + TotBiRem < NInteractions) {
        #----------------------------
        # Compute the quantiles and the expected PETs under H0:
        #----------------------------
        QuantRes = Get_ExpectedPETs_fun(InteractionInfMat = InteractionInfMat, AllInteIndeces = AllInteIndeces,
            NiNjMat = NiNjMat, QuantileProbs = QuantileProbs, NetUpdateIndicator = NetUpdateIndicator,
            NetworksData = NetworksData, SavedQuantiles = SavedQuantiles)
        # break output:
        BinMatVij = QuantRes$BinMatVij
        NetworksData = QuantRes$NetworksData
        InteractionInfMat = QuantRes$InteractionInfMat
        SavedQuantiles = QuantRes$SavedQuantiles
        #----------------------------
        # Generate p-values:
        #----------------------------
        pValues_round = BiocParallel::bplapply(X = as.list(AllInteIndeces), FUN = Assess_Interaction_fun_Rcpp,
            InteractionInfMat = InteractionInfMat, Poiss_fun = Poiss_fun, BinMatVij = BinMatVij)
        pValues_round = do.call(rbind, pValues_round)
        #----------------------------
        # get latest interaction information
        #----------------------------
        LaIn = Get_LatestInteraction_fun(pValues_round = pValues_round, S4_method = S4_method,
            TotPairs = sum(unlist(NetworksData$InteractionPairsList)), InteractionInfMat = InteractionInfMat)
        LastInteractions = LaIn$LastInteractions
        InteractionInfMat = LaIn$InteractionInfMat
        TR_Si = length(LastInteractions)  #total round significant interactions
        #----------------------------
        # Update the network for the newley added interactions
        #----------------------------
        NetworkUpdates = Update_Network_fun(AllInteIndeces = AllInteIndeces, TotIntAdded = TotIntAdded,
            TotBiRem = TotBiRem, NInteractions = NInteractions, LastInteractions = LastInteractions,
            TR_Si = TR_Si, InteractionInfMat = InteractionInfMat, OrdersCount = OrdersCount,
            NetworksData = NetworksData)
        NetworksData = NetworkUpdates$NetworksData
        TotIntAdded = NetworkUpdates$TotIntAdded
        TotBiRem = NetworkUpdates$TotBiRem
        AllInteIndeces = NetworkUpdates$AllInteIndeces
        InteractionInfMat = NetworkUpdates$InteractionInfMat
        NetUpdateIndicator = NetworkUpdates$NetUpdateIndicator
        # check if break:
        if (TotIntAdded + TotBiRem == NInteractions) {
            cat("\n")
            break
        }
        #----------------------------
        # increase counter:
        #----------------------------
        OrdersCount = OrdersCount + 1
    }
    #----------------------------
    # print and write:
    #----------------------------
    futile.logger::flog.info("Interaction analysis completed!", name = "SA_LogFile",
        capture = FALSE)
    #--------------------------------------------
    #------Delete matrices
    #--------------------------------------------
    unlink(x = S4_MatrixDir, recursive = TRUE, force = TRUE)
    #----------------------------
    # return:
    #----------------------------
    return(list(PeaksGR = PeaksGR, NInteractions = NInteractions, InteractionInfMat = InteractionInfMat,
        TotIntAdded = TotIntAdded, TotBiRem = TotBiRem))
}
# done
#----------------------------
# Function for finding the expected number of PETs if random
#----------------------------
Get_ExpectedPETs_fun = function(InteractionInfMat, AllInteIndeces, NiNjMat, QuantileProbs,
    NetUpdateIndicator, NetworksData, SavedQuantiles) {
    #----------------------------
    # Compute the SP for the need changed chromosomes:
    #----------------------------
    for (Net in NetUpdateIndicator) {
        # Gather network data:
        Network_Net = NetworksData$NetworkList[[Net]]
        UpdateIndeces_Net = NetworksData$UpdateIndecesList[[Net]]
        MergedNodesIndex_Net = NetworksData$MergedNodesIndexList[[Net]]
        BigInfoMatDesc_Net = NetworksData$BigInfoMatDescList[[Net]]
        Nadj_Net = NetworksData$NadjList[[Net]]
        NPeaksInvolved_Net = NetworksData$NPeaksInvolvedList[[Net]]
        NiNjIndeces_Net = NetworksData$NiNjIndecesList[[Net]]
        # update SP:
        InteractionPairs_Net = BiocParallel::bplapply(X = UpdateIndeces_Net, FUN = BigMat_SP_fun,
            NPeaksInvolved = NPeaksInvolved_Net, Nadj = Nadj_Net, BigInfoMatDesc = BigInfoMatDesc_Net,
            Network = Network_Net, MergedNodesIndex = MergedNodesIndex_Net, NiNjIndeces = NiNjIndeces_Net,
            NiNjMat = NiNjMat)
        InteractionPairs_Net = Reduce("+", InteractionPairs_Net)
        # Update InteractionPairs:
        NetworksData$InteractionPairsList[[Net]] = InteractionPairs_Net
    }
    # ---------------------------- Find quantiles for Vij=ni*nj/ Dij
    # ----------------------------
    BigMatQuantiles = BiocParallel::bplapply(X = as.list(NetUpdateIndicator), FUN = Get_QuantileChunks_fun,
        QuantileProbs = QuantileProbs, BigInfoMatDescList = NetworksData$BigInfoMatDescList)
    # save in the SavedQuantiles
    SavedQuantiles[NetUpdateIndicator] = BigMatQuantiles
    # split quantiles:
    BinsVij = Reduce("+", lapply(SavedQuantiles, "[[", 1))
    ChunkSizeVij = Reduce("+", lapply(SavedQuantiles, "[[", 2))  #the chunk total
    # finalize quantiles:
    BinsVij = unique(BinsVij/ChunkSizeVij)
    BinsVijSize = length(BinsVij)
    # ---------------------------- Load it from the data itself:
    # ----------------------------
    ObsVij = Get_Vij_Data_fun(InteractionInfMat = InteractionInfMat, NetworksData = NetworksData,
        AllInteIndeces = AllInteIndeces)
    #----------------------------
    # Count PETs and assign the QCells
    #----------------------------
    QCellPETCountsVij = numeric(BinsVijSize)
    Get_QCellPETCounts_fun_Rcpp(BinsVij, BinsVijSize, ObsVij, InteractionInfMat,
        AllInteIndeces, QCellPETCountsVij)
    #----------------------------
    # Find Cell counts, now you will loop all the networks:
    #----------------------------
    QCellCombCountsVij = numeric(BinsVijSize)
    for (Net in seq_len(NetworksData$TotNetworks)) {
        # Gather network data:
        UpdateIndeces_Net = NetworksData$UpdateIndecesList[[Net]]
        BigInfoMatDesc_Net = NetworksData$BigInfoMatDescList[[Net]]
        Nadj_Net = NetworksData$NadjList[[Net]]
        NPeaksInvolved_Net = NetworksData$NPeaksInvolvedList[[Net]]
        # Count in cells
        QCellCombCountsVij_Net = BiocParallel::bplapply(X = UpdateIndeces_Net, FUN = Get_QCellCombCounts_fun,
            BinsVij = BinsVij, BinsVijSize = BinsVijSize, BigInfoMatDesc = BigInfoMatDesc_Net,
            NPeaksInvolved = NPeaksInvolved_Net, Nadj = Nadj_Net)
        # Sum
        QCellCombCountsVij_Net = Reduce("+", QCellCombCountsVij_Net)
        QCellCombCountsVij = QCellCombCountsVij + QCellCombCountsVij_Net
    }
    #----------------------------
    # Create Bin Matrix for Vij
    #----------------------------
    BinMatVij = cbind(seq_len(BinsVijSize), BinsVij, QCellPETCountsVij/QCellCombCountsVij)
    colnames(BinMatVij) = NULL
    #----------------------------
    # Create and smooth bins
    #----------------------------
    BinMatVij = SmoothSplines_fun(BinMat = BinMatVij, Which = 3)
    # return:
    return(list(BinMatVij = BinMatVij, NetworksData = NetworksData,
        InteractionInfMat = InteractionInfMat, SavedQuantiles = SavedQuantiles))
}
# done
#----------------------------
# Function for filling the Big.matrix with the shortest paths.  IF the distance
# is zero, it becomes NA, as the nodes are merged
#----------------------------
BigMat_SP_fun = function(Indeces_x, NPeaksInvolved, Nadj, BigInfoMatDesc, Network,
    MergedNodesIndex, NiNjIndeces, NiNjMat) {
    #----------------------------
    # attach the matrices:
    #----------------------------
    BigInfoMatDescInst = bigmemory::attach.big.matrix(BigInfoMatDesc)
    InteractionPairs = numeric(1)  #the intra combinations
    #----------------------------
    # loop through the indeces:
    #----------------------------
    for (ind in Indeces_x) {
        #----------------------------
        # Take all the k from the MergedNodesIndex:
        #----------------------------
        k_vect = MergedNodesIndex[[ind]]
        if (is.na(k_vect[1]))
            next  #skip element if NA, merged, will be updated another time
        k_vect = sort(x = k_vect, decreasing = FALSE, na.last = TRUE)  #sort, need the bigger k in the beginning
        #----------------------------
        # Find SP for the first k index:
        #----------------------------
        GlobalNodesDist = Dijkstra_GSP_fun_Rcpp(k_vect[1], Network, NPeaksInvolved)
        #----------------------------
        # save the entries in the AdjucencyMatDescInst for all the k indeces
        #----------------------------
        for (k in k_vect) {
            #----------------------------
            # Get start and end of the saving indeces
            #----------------------------
            StartInd = Get_VectPosIndex_fun_Rcpp(NPeaksInvolved, Nadj, k - 1, k)  #c++
            EndInd = Get_VectPosIndex_fun_Rcpp(NPeaksInvolved, Nadj, k - 1, NPeaksInvolved -
                1)  #c++
            #----------------------------
            # Fill in the SP values:
            #----------------------------
            Save_BigMat_fun_Rcpp(BigInfoMatDescInst@address, GlobalNodesDist, k,
                StartInd, EndInd, InteractionPairs, NiNjIndeces, NiNjMat)
        }
    }
    return(InteractionPairs)
}
# done
#----------------------------
# Function for the Quantiles for Vij
#----------------------------
Get_QuantileChunks_fun = function(Net, QuantileProbs, BigInfoMatDescList) {
    #----------------------------
    # Attach and load matrix:
    #----------------------------
    BigInfoMatDescInst = bigmemory::attach.big.matrix(BigInfoMatDescList[[Net]])
    BigInfoMatDescInst = bigmemory::as.matrix(BigInfoMatDescInst)
    ChunkVij = which(!is.na(BigInfoMatDescInst))
    ChunkSizeVij = length(ChunkVij)
    #----------------------------
    # Get quantiles:
    #----------------------------
    if (ChunkSizeVij != 0) {
        # For Distance:
        BinsVij = stats::quantile(x = BigInfoMatDescInst[ChunkVij], probs = QuantileProbs,
            na.rm = TRUE, names = FALSE, type = 7) * ChunkSizeVij
    } else {
        BinsVij = rep(0, length(QuantileProbs))
    }
    return(list(BinsVij = BinsVij, ChunkSizeVij = ChunkSizeVij))
}
# done
#----------------------------
# Function for returning the observed Vij in the order of the interactions
#----------------------------
Get_Vij_Data_fun = function(InteractionInfMat, NetworksData, AllInteIndeces) {
    #----------------------------
    # First subset the data and keep what you need:
    #----------------------------
    InteractionInfMatSub = InteractionInfMat[AllInteIndeces, c("AdjNode_i", "AdjNode_j",
        "NiNjNode_i", "NiNjNode_j", "Chrom12ID", "IntraID")]
    if (!methods::is(InteractionInfMatSub, "matrix")){
        InteractionInfMatSub = t(as.matrix(InteractionInfMatSub))
    }
    InteractionInfMatSub = InteractionInfMatSub[which(InteractionInfMatSub[, c("IntraID")] ==
        1), ]  #keep intra in that part
    if (!methods::is(InteractionInfMatSub, "matrix")){
        InteractionInfMatSub = t(as.matrix(InteractionInfMatSub))
    }
    #----------------------------
    # Create ObsDVij vector and loop to update
    #----------------------------
    ObsVij = numeric(length(AllInteIndeces))
    Chrom12IDUnique = unique(InteractionInfMatSub[, c("Chrom12ID")])
    for (Net in Chrom12IDUnique) {
        # attach:
        BigInfoMatDescInst = bigmemory::attach.big.matrix(NetworksData$BigInfoMatDescList[[Net]])
        # Get AdjNodes/NiNjNode:
        Pos_Net = which(InteractionInfMatSub[, c("Chrom12ID")] == Net)
        AdjNode_i_Net = InteractionInfMatSub[Pos_Net, c("AdjNode_i")]
        AdjNode_j_Net = InteractionInfMatSub[Pos_Net, c("AdjNode_j")]
        # Get Nadj and NPeaks:
        Nadj_Net = NetworksData$NadjList[[Net]]
        NPeaksInvolved_Net = NetworksData$NPeaksInvolvedList[[Net]]
        # Get vect position on the matrix in R indeces!
        AdjNode_ij_Net = Get_VectPosIndex_Vectorized_fun_Rcpp(NPeaksInvolved_Net,
            Nadj_Net, AdjNode_i_Net, AdjNode_j_Net)
        # load the Vij and save them:
        ObsVij[Pos_Net] = BigInfoMatDescInst[AdjNode_ij_Net, c(1)]
    }
    return(ObsVij)
}
# done
#----------------------------
# Function for finding the QCellIDs and counting in the cells
#----------------------------
Get_QCellCombCounts_fun = function(Indeces_x, BinsVij, BinsVijSize, BigInfoMatDesc,
    NPeaksInvolved, Nadj) {
    #----------------------------
    # attach the matrices:
    #----------------------------
    BigInfoMatDescInst = bigmemory::attach.big.matrix(BigInfoMatDesc)
    #----------------------------
    # Initialize:
    #----------------------------
    QCellCombCountsVij_Net = numeric(BinsVijSize)
    #----------------------------
    # loop through the indeces:
    #----------------------------
    for (ind in Indeces_x) {
        #----------------------------
        # Find start and end indeces on the matrices
        #----------------------------
        StartInd = Get_VectPosIndex_fun_Rcpp(NPeaksInvolved, Nadj, ind - 1, ind)  #c++
        EndInd = Get_VectPosIndex_fun_Rcpp(NPeaksInvolved, Nadj, ind - 1, NPeaksInvolved -
            1)  #c++
        #----------------------------
        # load the Dij and take order:
        #----------------------------
        VkhOrder = order(x = BigInfoMatDescInst[seq(from = StartInd + 1, to = EndInd +
            1), c(1)], na.last = TRUE, decreasing = FALSE)
        #----------------------------
        # Fill inn the values:
        #----------------------------
        Get_QCellCombCounts_fun_Rcpp(BinsVij, BinsVijSize, BigInfoMatDescInst@address,
            VkhOrder, QCellCombCountsVij_Net, StartInd, EndInd)
    }
    return(QCellCombCountsVij_Net)
}
# done
#----------------------------
# Splines
#----------------------------
SmoothSplines_fun = function(BinMat, Which) {
    #----------------------------
    # replace inf with 0:
    #----------------------------
    InfBins = which(is.infinite(BinMat[, Which]))
    if (length(InfBins) != 0)
        BinMat[InfBins, Which] = 0
    #----------------------------
    # Remove Na bins, those are with no pairs in them, if you set them to zero they
    # will affect the nearby values So dont consider them just.
    #----------------------------
    NoNABins = which(!is.na(BinMat[, Which]))  #the non NA
    #----------------------------
    # smooth with splines:
    #----------------------------
    if (length(NoNABins) > 4) {
        # keep cells which are used because the lambda ij should not become zero:
        UsedQuantiles = which(BinMat[NoNABins, Which] != 0)  #cells that are used
        # splines
        Splines = stats::smooth.spline(x = log10(BinMat[NoNABins, c(2)]), y = BinMat[NoNABins,
            Which], all.knots = TRUE, spar = 0.75)
        # fix zero:
        if (any(Splines$y < 0))
            Splines$y[which(Splines$y < 0)] = 0
        BinMat[NoNABins, Which] = Splines$y
        #----------------------------
        # Splines might set a Cell which lij!=0, to lij==0, and this will be a problem
        # with the p-values So fix those by setting them to .Machine$eps
        #----------------------------
        BinMat[NoNABins[UsedQuantiles], Which] = pmax(.Machine$double.eps, BinMat[NoNABins[UsedQuantiles],
            Which])
    }
    return(BinMat)
}
# done
#----------------------------
# Poisson p-value
#----------------------------
Poiss_fun = function(nij, lij) {
    pval = stats::dpois(x = nij, lambda = lij, log = FALSE) + stats::ppois(q = nij,
        lambda = lij, lower.tail = FALSE, log.p = FALSE)
    return(pval)
}
# done
#----------------------------
# Function for finding which interaction IDs will be added based on their FDR
#----------------------------
Get_LatestInteraction_fun = function(pValues_round, S4_method, TotPairs, InteractionInfMat) {
    #----------------------------
    # Adjust with FDR, note that the p-values of the non-interacting are also needed,
    # which are pval=1:
    #----------------------------
    FDR_round = stats::p.adjust(p = pValues_round[, c(2)], method = S4_method, n = TotPairs)
    #----------------------------
    # Update the InteractionsInfMat too:
    #----------------------------
    InteractionInfMat[pValues_round[,c(1)], "pvalue"] = pValues_round[, c(2)]
    InteractionInfMat[pValues_round[,c(1)], "FDR"] = FDR_round
    #----------------------------
    # Find the most significant FDR and the interactions who have it
    #----------------------------
    Min_round_FDR = min(FDR_round)
    WhichMinFDR_round = which(FDR_round == Min_round_FDR)
    LastInteractions = pValues_round[WhichMinFDR_round, c(1)]  #R index
    return(list(LastInteractions = LastInteractions, InteractionInfMat = InteractionInfMat))
}
# Done
#----------------------------
# Function for updating the network
#----------------------------
Update_Network_fun = function(AllInteIndeces, TotIntAdded, TotBiRem, NInteractions,
    LastInteractions, TR_Si, InteractionInfMat, OrdersCount, NetworksData) {
    #----------------------------
    # First remove all the interactions to be added from the interaction from
    # looping, and initiate which chrom comb will be updated next:
    #----------------------------
    AllInteIndeces = AllInteIndeces[-which(AllInteIndeces %in% LastInteractions)]
    NetUpdateIndicator = c()
    #----------------------------
    # loop on the interactions to add:
    #----------------------------
    for (i in seq_len(TR_Si)) {
        #----------------------------
        # print:
        #----------------------------
        TotIntAdded = TotIntAdded + 1
        TotIntAddedPrintPers = TotIntAdded/NInteractions * 100
        cat("|---- Total interactions processed: ", TotIntAdded, " (", TotIntAddedPrintPers,
            "%) ----|\r")
        #----------------------------
        # Take interaction information ID
        #----------------------------
        ID_i = LastInteractions[i]  #row ID of InteractionInfMat
        #----------------------------
        # update InteractionInfMat:
        #----------------------------
        InteractionInfMat[ID_i, c("Order")] = OrdersCount
        Chrom12ID_i = InteractionInfMat[ID_i, c("Chrom12ID")]  #for updating paths again
        NetUpdateIndicator = c(NetUpdateIndicator, Chrom12ID_i)
        #----------------------------
        # take the nodes, you will merge the biggest to the smallest in the network:
        #----------------------------
        Node_kh = sort(InteractionInfMat[ID_i, c("AdjBiNode_i", "AdjBiNode_j")],
            decreasing = FALSE)
        # Check if nodes are merged already(not this is NOT a bi-product)
        if (Node_kh[1] == Node_kh[2])
            next  #everything is done so go to next
        #----------------------------
        # Move all nodes to the MergedNodesIndex:
        #----------------------------
        k = Node_kh[1]
        names(k) = NULL
        h = Node_kh[2]
        names(h) = NULL
        NetworksData$MergedNodesIndexList[[Chrom12ID_i]][[k]] = c(NetworksData$MergedNodesIndexList[[Chrom12ID_i]][[k]],
            NetworksData$MergedNodesIndexList[[Chrom12ID_i]][[h]])
        NetworksData$MergedNodesIndexList[[Chrom12ID_i]][[h]] = NA
        #----------------------------
        # Update the Nodes of the rest of the interactions to be added: Those are the
        # interactions which are significant this round but not added yet, not the rest
        # of the interactions left.
        #----------------------------
        Update_ToBeAddedInter_fun_Rcpp(InteractionInfMat, k, h, i, TR_Si, LastInteractions,
            Chrom12ID_i)
        #----------------------------
        # update the rest of the interactions, not added yet.  update the index and check
        # if bi-products (same index) The bi-products are removed from the model, their
        # p-value/FDR will be NA and the order will be NA too Their index will be removed
        # from the AllInteIndeces
        #----------------------------
        BiPorductsInfo = Check_BiProd_fun_Rcpp(InteractionInfMat, k, h, AllInteIndeces,
            TotBiRem, Chrom12ID_i, OrdersCount)
        # Check bi-product rejected:
        if (BiPorductsInfo$TotBiRem != TotBiRem) {
            # then at least one bi-product. remove from the rest of the interactions
            AllInteIndeces = AllInteIndeces[-which(AllInteIndeces %in% BiPorductsInfo$BiProductIDSreject)]
            TotBiRem = BiPorductsInfo$TotBiRem
        }
        # check bi-products accepted:
        if (BiPorductsInfo$TotBiAcc != 0) {
            # then at least one bi-product. remove from the rest of the interactions
            AllInteIndeces = AllInteIndeces[-which(AllInteIndeces %in% BiPorductsInfo$BiProductIDSaccepted)]
            TotIntAdded = TotIntAdded + BiPorductsInfo$TotBiAcc
        }
        #----------------------------
        # Update the network:
        #----------------------------
        Network_i = NetworksData$NetworkList[[Chrom12ID_i]]
        Network_i = Update_Network_kh_fun(Network = Network_i, k = k, h = h)
        NetworksData$NetworkList[[Chrom12ID_i]] = Network_i
    }
    # finalize NetUpdateIndicator:
    NetUpdateIndicator = sort(unique(NetUpdateIndicator))
    #----------------------------
    # return
    #----------------------------
    return(list(NetworksData = NetworksData, TotIntAdded = TotIntAdded, TotBiRem = TotBiRem,
        AllInteIndeces = AllInteIndeces, InteractionInfMat = InteractionInfMat,
        NetUpdateIndicator = NetUpdateIndicator))
}
# done
#----------------------------
# Function to update the Network by merging h to k node
#----------------------------
Update_Network_kh_fun = function(Network, k, h) {
    # --------------------- NOTE: First move all h edges to k (in commons keep min
    # dist) In h set k edge only, weight 0. (dont delete element since then you need
    # to reduce all indeces which takes time) (just leave it empty pointing to k
    # only) Dont remove the h edge from k! Because you need the distance The
    # InteractionsMat will not have h anymore since it is replaced with k for all.
    # All edges connected to h should also be traversed and the h name replaced to k
    # (and the dist to min again.)  -------------------- Take the k and h networks:
    # --------------------- k:
    Network_k = Network[[k]]
    edges_k = Network_k$edges
    weights_k = Network_k$weights
    # h:
    Network_h = Network[[h]]
    edges_h = Network_h$edges
    weights_h = Network_h$weights
    # -------------------- Remove h edge from k if there ---------------------
    if (h %in% edges_k) {
        Rmh = which(edges_k == h)
        edges_k = edges_k[-c(Rmh)]
        weights_k = weights_k[-c(Rmh)]
    }
    # -------------------- Remove k edge from h if there ---------------------
    if (k %in% edges_h) {
        Rmk = which(edges_h == k)
        edges_h = edges_h[-c(Rmk)]
        weights_h = weights_h[-c(Rmk)]
    }
    # -------------------- update the common elements to the min
    # ---------------------
    if (any(edges_h %in% edges_k)) {
        # take common edges
        edges_hk = intersect(edges_h, edges_k)
        # take weights of h:
        Common_h = which(edges_h %in% edges_hk)
        weights_hk_h = weights_h[Common_h]
        # take weights of k:
        Common_k = which(edges_k %in% edges_hk)
        weights_hk_k = weights_k[Common_k]
        # take pairwise min:
        weights_min = pmin(weights_hk_k, weights_hk_h)
        # update:
        weights_k[Common_k] = weights_min
    }
    # -------------------- Add exclusive elements of h into k ---------------------
    if (any(!edges_h %in% edges_k)) {
        Exclusive_h = which(!edges_h %in% edges_k)
        edges_k = c(edges_k, edges_h[Exclusive_h])
        weights_k = c(weights_k, weights_h[Exclusive_h])
    }
    # -------------------- Add h to the k elements and sort them:
    # ---------------------
    edges_k = c(edges_k, h)
    weights_k = c(weights_k, 0)
    Orderk = order(weights_k, decreasing = FALSE)
    edges_k = edges_k[Orderk]
    weights_k = weights_k[Orderk]
    # -------------------- Go into the elements of h and update h edge to k This will
    # reduce the recurrence ---------------------
    for (el in edges_h) {
        # note k is removed from edges_h Also all the el are in k, so set k in el too
        # with same distance take network:
        Network_el = Network[[el]]
        edges_el = Network_el$edges
        weights_el = Network_el$weights
        # remove the k and h edges completely:
        Rmkh = which(edges_el %in% c(k, h))
        edges_el = edges_el[-c(Rmkh)]
        weights_el = weights_el[-c(Rmkh)]
        # add the k edge, weights_k[Addk] has the min already:
        Addk = which(edges_k %in% c(el))
        edges_el = c(edges_el, k)
        weights_el = c(weights_el, weights_k[Addk])
        # sort:
        Orderel = order(weights_el, decreasing = FALSE)
        edges_el = edges_el[Orderel]
        weights_el = weights_el[Orderel]
        # save to network:
        Network_el$edges = edges_el
        Network_el$weights = weights_el
        Network[[el]] = Network_el
    }
    # -------------------- update network for k and h: ---------------------
    Network_k$edges = edges_k
    Network_k$weights = weights_k
    Network[[k]] = Network_k
    Network_h$edges = k
    Network_h$weights = 0
    Network[[h]] = Network_h
    return(Network)
}
# done
#----------------------------
# Function for creating the genome map object
#----------------------------
Create_GenomeMapObject = function(InteractionResults, SA_prefix, S4_AnalysisDir) {
    cat("Creating the GenomeMap Object...")
    # R check:
    Order = NULL
    #----------------------------
    # Create the significant information matrix nrow=Tot_interactions, col1=PBS from,
    # col2=PBS to (rows in the PeaksData) col3=pvalue, col4=FDR, col5=order, col6=Tot
    # interaction PETs between the two PBS
    #----------------------------
    InteractionInfo = Get_InteractionInfo_fun_Rcpp(InteractionResults$InteractionInfMat,
        InteractionResults$NInteractions)
    InteractionInfo = as.data.frame(InteractionInfo)
    colnames(InteractionInfo) = c("PBSF", "PBST", "pvalue", "FDR", "Order", "TotalInterPETs")
    # subset the order being >0, remove not added interactions/bi products:
    InteractionInfo = subset(InteractionInfo, !(is.na(Order) | Order == 0))
    # sort by order:
    GetOrder = order(InteractionInfo$Order)
    InteractionInfo = InteractionInfo[GetOrder, ]
    #----------------------------
    # create the GRanges object from the peaks:
    #----------------------------
    # take the PeaksGR:
    PeaksGR = InteractionResults$PeaksGR
    # reduce unused levels:
    LevelsUsed = GenomeInfoDb::seqlevelsInUse(PeaksGR)
    GenomeInfoDb::seqlevels(PeaksGR) = LevelsUsed
    #----------------------------
    # create the GInteractions object from the PeaksGR:
    #----------------------------
    AnchorSummits = round(PeaksGR$PeakSummit)
    PeaksGR$PeakSummit = NULL
    PeaksGR$LID = NULL
    # Take Anchor 1:
    Anchor1 = PeaksGR[InteractionInfo$PBSF]
    Anchor1Summit = AnchorSummits[InteractionInfo$PBSF]
    # Take Anchor 2:
    Anchor2 = PeaksGR[InteractionInfo$PBST]
    Anchor2Summit = AnchorSummits[InteractionInfo$PBST]
    #----------------------------
    # create the GenomeMapData
    #----------------------------
    GenomeMapData = InteractionSet::GInteractions(anchor1 = Anchor1, anchor2 = Anchor2)
    # Add the summits:
    GenomeMapData$Anchor1Summit = Anchor1Summit
    GenomeMapData$Anchor2Summit = Anchor2Summit
    # make InteractionInfo DataFrame:
    InteractionInfo = S4Vectors::DataFrame(InteractionInfo[, c("pvalue", "FDR", "Order",
        "TotalInterPETs")])
    #----------------------------
    # Add Metadata:
    #----------------------------
    S4Vectors::metadata(GenomeMapData)$InteractionInfo = InteractionInfo
    #----------------------------
    # Convert To class and save:
    #----------------------------
    class(GenomeMapData) = "GenomeMap"
    NameGenomeMapData = paste(SA_prefix, "_GenomeMapData", sep = "")
    assign(NameGenomeMapData, GenomeMapData)  #assign value.
    save(list = NameGenomeMapData, file = file.path(S4_AnalysisDir, NameGenomeMapData))
    cat("Done\n")
}
# done
