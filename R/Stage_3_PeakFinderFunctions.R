# Helping functions for PeakFinder.R
#' @importFrom plyr dlply . ldply llply ddply
#' @importFrom intervals Intervals clusters interval_included size
#' @importFrom stats ppois p.adjust
#' @importFrom S4Vectors metadata
#' @importFrom GenomeInfoDb seqlengths seqinfo
#' @importFrom BiocParallel bplapply
#' @importClassesFrom S4Vectors Annotated
#' @importClassesFrom intervals Intervals
############################################## Main function for stage 3:
#-------------
#-------------
Stage_3_Main_fun = function(SA_prefix, S3_method, S3_AnalysisDir, S3_Selfobject, 
    S3_image) {
    # global variables for Rcheck:
    #--------------------------------------------
    #---------------Take and reorder Input:
    #--------------------------------------------
    # Take time:
    Analysis.time.start = Sys.time()
    # create folder to save data:
    if (!dir.exists(S3_AnalysisDir)) 
        dir.create(S3_AnalysisDir)
    #--------------------------------------------
    #---------------Keep data you need only:
    #--------------------------------------------
    FitDataInfo = TakeFitAnalysisData_fun(S3_Selfobject = S3_Selfobject)
    PETsData = FitDataInfo$PETsData  #all the PETs data, used in inference too.
    ChromInf = FitDataInfo$ChromInf  #used in inference and fitting.
    #--------------------------------------------
    #---------------Segment to regions:
    #--------------------------------------------
    SegmRes = FindRegions_fun(PETsData = PETsData, ChromInf = ChromInf)
    SegPETS = SegmRes$SegPETS  #segmented PETs data, used in fitting.
    RegionCounts = SegmRes$RegionCounts  #take the counts
    #------------
    #-----save the region counts on main data:
    #------------
    Ordermatch = match(RegionCounts$Chrom, S4Vectors::metadata(S3_Selfobject)$Self_info$Chrom)
    S4Vectors::metadata(S3_Selfobject)$Self_info$Region.counts[Ordermatch] = RegionCounts$Counts
    #--------------------------------------------
    #---------------Fit:
    #--------------------------------------------
    GlobalPeakCallRes = FitCallGlobal_fun(SegPETS = SegPETS, ChromInf = ChromInf)
    #--------------------------------------------
    #---------------update the results
    #--------------------------------------------
    #------classification information
    S4Vectors::metadata(S3_Selfobject)$Classification.Info = GlobalPeakCallRes$Classification.Info
    #------Peak counts
    S4Vectors::metadata(S3_Selfobject)$Self_info$Peak.counts = 0
    Match = match(GlobalPeakCallRes$Peak.counts$Chrom, S4Vectors::metadata(S3_Selfobject)$Self_info$Chrom)
    S4Vectors::metadata(S3_Selfobject)$Self_info$Peak.counts[Match] = GlobalPeakCallRes$Peak.counts$V1
    #--------------------------------------------
    #---------------Run inference:
    #--------------------------------------------
    Peaks.Info = GlobalPeakCallRes$Peaks.Info  #peaks found by fit
    # run inference
    InfRes = SignificanceCall_fun(Peaks.Info = Peaks.Info, PETsData = PETsData, ChromInf = ChromInf, 
        S3_method = S3_method)
    #---add Peak information for peaks found in data with their FDR etc
    S4Vectors::metadata(S3_Selfobject)$Peaks.Info = InfRes
    #------------------------
    #-----update class and save:
    #------------------------
    class(S3_Selfobject) = "PSFit"
    NamepsfitData = paste(SA_prefix, "_psfitData", sep = "")
    assign(NamepsfitData, S3_Selfobject)  #assign value.
    save(list = NamepsfitData, file = file.path(S3_AnalysisDir, NamepsfitData))
    #------------------------
    #-----plot images
    #------------------------
    if (S3_image) {
        #-------------
        # create data:
        #-------------
        PeakInfo = S4Vectors::metadata(S3_Selfobject)
        PeakInfo = PeakInfo$Peaks.Info
        PeakInfo$FDR0.05 = "FDR >= 0.05"
        PeakInfo$FDR0.05[which(PeakInfo$FDR < 0.05)] = "FDR < 0.05"
        PeakInfo$ID = 1:nrow(PeakInfo)
        # image 1:
        Get_image_S3_P1_fun(S3_AnalysisDir = S3_AnalysisDir, SA_prefix = SA_prefix, 
            PeakInfo = PeakInfo)
        # image 2:
        Get_image_S3_P2_fun(S3_AnalysisDir = S3_AnalysisDir, SA_prefix = SA_prefix, 
            PeakInfo = PeakInfo)
        # image 3:
        Get_image_S3_P3_fun(S3_AnalysisDir = S3_AnalysisDir, SA_prefix = SA_prefix, 
            PeakInfo = PeakInfo)
        # image 4:
        Get_image_S3_P4_fun(S3_AnalysisDir = S3_AnalysisDir, SA_prefix = SA_prefix, 
            PeakInfo = PeakInfo)
    }
    #------------------------
    # save log:
    #------------------------
    futile.logger::flog.info("=====================================", name = "SA_LogFile", 
        capture = FALSE)
    futile.logger::flog.info("Stage 3 is done!", name = "SA_LogFile", capture = FALSE)
    futile.logger::flog.info(paste("Analysis results for stage 3 are in:\n", S3_AnalysisDir), 
        name = "SA_LogFile", capture = FALSE)
    # save time:
    Analysis.time.end = Sys.time()
    Total.Time = Analysis.time.end - Analysis.time.start
    LogFile = paste("Total stage 3 time:", Total.Time, " ", units(Total.Time))
    futile.logger::flog.info(LogFile, name = "SA_LogFile", capture = FALSE)
}
# done
#-------------
#-------------
############## Functions for inputs and data
#-------------
#-------------
# function for converting the data in the analysis format:
TakeFitAnalysisData_fun = function(S3_Selfobject) {
    cat("Converting data for analysis...")
    #------------
    # take Chromosome information(some used in inference):
    #------------
    ChromInf = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(S3_Selfobject))
    ChromInf = data.frame(Chrom = names(ChromInf), size = as.numeric(ChromInf), stringsAsFactors = FALSE)
    Self_info = S4Vectors::metadata(S3_Selfobject)$Self_info
    MatchChrom = match(Self_info$Chrom, ChromInf$Chrom)
    ChromInf$PET.counts[MatchChrom] = Self_info$PET.counts
    #------------
    # convert to data frame and Swap Anchors:
    #------------
    PETsData = as.data.frame(S3_Selfobject)
    # use Tag mids:
    PETsData$UTag = (PETsData$start1 + PETsData$end1)/2
    PETsData$DTag = (PETsData$start2 + PETsData$end2)/2
    Swap = which(PETsData$UTag > PETsData$DTag)
    if (length(Swap) != 0) {
        UTagswap = PETsData$UTag[Swap]
        DTagswap = PETsData$DTag[Swap]
        PETsData$UTag[Swap] = DTagswap
        PETsData$DTag[Swap] = UTagswap
    }
    #------------
    # add mainindex
    #------------
    PETsData$MainIndex = 1:nrow(PETsData)
    #------------
    # keep left part of downstream and right of upstream
    #------------
    PETsData$Chrom = as.character(PETsData$seqnames1)
    PETsData = PETsData[, c("Chrom", "UTag", "DTag", "MainIndex")]
    #------------
    # return:
    #------------
    cat("Done\n")
    return(list(PETsData = PETsData, ChromInf = ChromInf))
}
# done
#-------------
#-------------
############## Functions for region segmentation:
#-------------
#-------------
# main function for region segmentation:
FindRegions_fun = function(PETsData, ChromInf) {
    # rcheck:
    Chrom = Region = NULL
    #------------
    #------------
    # split by chromosome to list to pass to bplapply
    bppass = plyr::dlply(PETsData, plyr::.(Chrom), function(x) x)
    #------------
    #------------
    #--Segment to regions:
    #------------
    cat("Segmenting into regions...")
    SegmRes = BiocParallel::bplapply(X = bppass, FUN = function(y, ChromInf) {
        ChromCur = as.character(unique(y$Chrom))
        ChromInfCur = subset(ChromInf, Chrom == ChromCur)
        #------------
        # Give the region ids:
        #------------
        y = RegionIds_fun(PETinf = y, ChromCur = ChromCur, ChromInfCur = ChromInfCur)
        return(y)
    }, ChromInf = ChromInf)
    cat("Done\n")
    #------------
    #---break and rbind:
    #------------
    # take the region info:
    PETinf = lapply(SegmRes, "[[", 1)
    PETinf = plyr::ldply(PETinf, function(y) y)
    PETinf = subset(PETinf, !is.na(Region))  #keep not NA
    if (nrow(PETinf) != 0) {
        # update data with the regions:
        PETsData$Region = NA
        PETsData$Region[PETinf$MainIndex] = PETinf$Region
        PETsData = subset(PETsData, !is.na(Region))
    } else {
        stop("No regions found in data. The data is probably too sparse.", call. = FALSE)
    }
    #------------
    #--take total regions info:
    #------------
    RegionCounts = lapply(SegmRes, "[[", 2)
    RegionCounts = plyr::ldply(RegionCounts, function(y) y)
    RegionCounts = RegionCounts[, c("Chrom", "Counts")]
    futile.logger::flog.info(paste("Total Regions found:", sum(RegionCounts$Counts)), 
        name = "SA_LogFile", capture = FALSE)
    return(list(SegPETS = PETsData, RegionCounts = RegionCounts))
}
# done
#-------------
#-------------
# Fuction for creating region ids by segmenting into regions:
RegionIds_fun = function(PETinf, ChromCur, ChromInfCur) {
    # global variables for Rcheck:
    #--------------------------
    #------------
    # initiate
    #------------
    PETinf$Region = NA
    TotRegions = data.frame(Chrom = ChromCur, Counts = 0)
    #------------
    # take Pet Pmid to create global intervals:
    #------------
    PETInt = intervals::Intervals(PETinf[, c("UTag", "DTag")], closed = c(TRUE, TRUE))
    #------------
    # Take clusters:
    #------------
    PETclst = intervals::clusters(PETInt, which = TRUE, w = 0)
    TotRegions$Counts = length(PETclst)  #total regions
    if (TotRegions$Counts != 0) 
        {
            #------------
            # region ids:
            #------------
            Region = rep(1:length(PETclst), lengths(PETclst))
            #------------
            # pets in the regions:
            #------------
            TagsIncluded = unlist(PETclst)
            #------------
            # add them according to the rownames of the PETinf
            #------------
            PETinf$Region[TagsIncluded] = Region
        }  #else no region found
    PETinf = PETinf[, c("Chrom", "MainIndex", "Region")]
    return(list(PETinf = PETinf, TotRegions = TotRegions))
}
# done
#-------------
#-------------
############## Functions for fitting the main model:
#-------------
#-------------
# main function for fitting: breaking each region
FitCallGlobal_fun = function(SegPETS, ChromInf) {
    #---R-check
    Region = Chrom = NULL
    #----
    bppass = plyr::dlply(SegPETS, plyr::.(Chrom, Region), function(y) y)
    #---------------------------------------
    #--------call the function in parallel for fitting regions:
    #---------------------------------------
    cat("Running peak calling process...")
    # the following runs in parallel:
    GlobalFitRes = BiocParallel::bplapply(X = bppass, FUN = FitCallLocal_fun_Rcpp, 
        ChromInf = ChromInf)
    cat("Done\n")
    #---------------------------
    #--------save Classification:
    #---------------------------
    Peak.Id.Inf = plyr::llply(GlobalFitRes, function(x) x[[1]])
    Peak.Id.Inf = do.call(rbind, Peak.Id.Inf)
    if (is.null(Peak.Id.Inf)) {
        stop("No Peaks found by the algorithm!", call. = FALSE)
    }
    #---------------------------
    #--match classes with SegPETS
    #---------------------------
    SegPETS$Peak.ID = NA
    Match = match(Peak.Id.Inf[, 2], SegPETS$MainIndex)
    SegPETS$Peak.ID[Match] = Peak.Id.Inf[, 1]
    Classification.Info = SegPETS[, c("MainIndex", "Region", "Peak.ID")]
    #---------------------------
    #--------save Peaks:
    #---------------------------
    Peaks.Info = plyr::ldply(GlobalFitRes, function(x) x[[2]])
    if (is.null(Peaks.Info)) {
        stop("No Peaks found by the algorithm!", call. = FALSE)
    }
    Peaks.Info = Peaks.Info[, which(colnames(Peaks.Info) != ".id")]
    futile.logger::flog.info(paste("Total", nrow(Peaks.Info), "candidate peaks found in data"), 
        name = "SA_LogFile", capture = FALSE)
    #---------------------------
    #--------Find counts:
    #---------------------------
    Peak.counts = plyr::ddply(Peaks.Info, plyr::.(Chrom), nrow)
    return(list(Peaks.Info = Peaks.Info, Classification.Info = Classification.Info, 
        Peak.counts = Peak.counts))
}
# done
#-------------
#-------------
############## Functions for the Inference:
#-------------
#-------------
# Main significance call function:
SignificanceCall_fun = function(Peaks.Info, PETsData, ChromInf, S3_method) {
    # global variables for Rcheck:
    Chrom = NULL
    #--------------------------
    #--------------------------------------------
    #---------------Take information for the analysis:
    #--------------------------------------------
    Peaks.Info$Chrom = as.character(Peaks.Info$Chrom)
    PETsData$MainIndex = NULL  #remove index not needed
    #--------------------------------------------
    #---------------run inference:
    #--------------------------------------------
    cat("Splitting data by chromosome for inference...\n")
    # split by chromosome all the three inputs: then merge them by chromosome
    bppass = Get_SignificanceCall_Data_fun(Peaks.Info = Peaks.Info, PETsData = PETsData, 
        ChromInf = ChromInf)
    #------------
    # Run inference(call c++):
    #------------
    futile.logger::flog.info("Computing p-values...", name = "SA_LogFile", capture = FALSE)
    InfRes = BiocParallel::bplapply(X = bppass, FUN = PoissonLocalInference_fun, 
        windows = c(10, 15))
    #------------
    # merge and return:
    #------------
    InfRes = do.call(rbind.data.frame, InfRes)
    rownames(InfRes) = NULL
    #------
    # FDR at 0.05:
    #------
    futile.logger::flog.info("FDR adjusting p-values...", name = "SA_LogFile", capture = FALSE)
    InfRes$FDRUp = stats::p.adjust(p = InfRes$p.valueUp, method = S3_method)
    InfRes$FDRDown = stats::p.adjust(p = InfRes$p.valueDown, method = S3_method)
    InfRes$FDR = stats::p.adjust(p = InfRes$p.value, method = S3_method)
    #------
    #---return:
    #------
    return(InfRes)
}
# done
#-------------
#-------------
# Function for breaking the inputs for significance calling:
Get_SignificanceCall_Data_fun = function(Peaks.Info, PETsData, ChromInf) {
    # Rcheck:
    Chrom = NULL
    #-------------
    # break Peaks.Info by chrom:
    #-------------
    Peaks.Info_list = plyr::dlply(Peaks.Info, plyr::.(Chrom), function(x) x)
    #-------------
    # break PETsData by chrom:
    #-------------
    PETsData$MainIndex = NULL  #dont need index
    PETsData_list = plyr::dlply(PETsData, plyr::.(Chrom), function(x) x)
    #-------------
    # break ChromInf by chrom:
    #-------------
    ChromInf_list = plyr::dlply(ChromInf, plyr::.(Chrom), function(x) x)
    #-------------
    # create the bppass input:
    #-------------
    bppass = list()
    for (chri in seq_along(Peaks.Info_list)) {
        # take Peaks_Info_x:
        Peaks_Info_x = Peaks.Info_list[[chri]]
        Chr_x = Peaks_Info_x$Chrom[1]
        # take PETsData_x:
        PETsData_x = PETsData_list[[which(names(PETsData_list) == Chr_x)]]
        PETsData_x$Chrom = NULL  #dont need that anymore
        PETsData_x = as.matrix(PETsData_x)
        # take ChromInf_x:
        ChromInf_x = ChromInf_list[[which(names(ChromInf_list) == Chr_x)]]
        ChromInf_x$Chrom = NULL  #dont need that
        bppass[[chri]] = list(Peaks_Info_x = Peaks_Info_x, PETsData_x = PETsData_x, 
            ChromInf_x = ChromInf_x)
    }
    return(bppass)
}
# done
#-------------
#-------------
# Function for Inference:
PoissonLocalInference_fun = function(bppass_x, windows) {
    #-------------
    # break data:
    #-------------
    Peaks_Info_x = bppass_x$Peaks_Info_x
    PETsData_x = bppass_x$PETsData_x
    # take UTag
    UTag_x = intervals::Intervals(PETsData_x[, c(1, 1)], closed = c(TRUE, TRUE))
    # take DTag
    DTag_x = intervals::Intervals(PETsData_x[, c(2, 2)], closed = c(TRUE, TRUE))
    # take chromosome info:
    ChromSize = bppass_x$ChromInf_x$size
    ChromPETs = bppass_x$ChromInf_x$PET.counts
    #-------------
    # ------ Up window 1:
    #-------------
    UpW1 = cbind(Peaks_Info_x$CIQ.Up.start - Peaks_Info_x$CIQ.Up.size * windows[1]/2, 
        Peaks_Info_x$CIQ.Up.start + Peaks_Info_x$CIQ.Up.size * windows[1]/2)
    UpW1 = intervals::Intervals(UpW1, closed = c(TRUE, TRUE))
    # observed pets in intervals
    UpW1_obs = lengths(intervals::interval_included(UpW1, UTag_x))
    lambdaUpW1 = UpW1_obs * Peaks_Info_x$CIQ.Up.size/(intervals::size(UpW1) + 1)
    #-------------
    # ------ Up window 2:
    #-------------
    UpW2 = cbind(Peaks_Info_x$CIQ.Up.start - Peaks_Info_x$CIQ.Up.size * windows[2]/2, 
        Peaks_Info_x$CIQ.Up.start + Peaks_Info_x$CIQ.Up.size * windows[2]/2)
    UpW2 = intervals::Intervals(UpW2, closed = c(TRUE, TRUE))
    # observed pets in intervals
    UpW2_obs = lengths(intervals::interval_included(UpW2, UTag_x))
    lambdaUpW2 = UpW2_obs * Peaks_Info_x$CIQ.Up.size/(intervals::size(UpW2) + 1)
    #-------------
    # ------ Down window 1:
    #-------------
    DownW1 = cbind(Peaks_Info_x$CIQ.Down.start - Peaks_Info_x$CIQ.Down.size * windows[1]/2, 
        Peaks_Info_x$CIQ.Down.start + Peaks_Info_x$CIQ.Down.size * windows[1]/2)
    DownW1 = intervals::Intervals(DownW1, closed = c(TRUE, TRUE))
    # observed pets in intervals
    DownW1_obs = lengths(intervals::interval_included(DownW1, DTag_x))
    lambdaDownW1 = DownW1_obs * Peaks_Info_x$CIQ.Down.size/(intervals::size(DownW1) + 
        1)
    #-------------
    # ------ Down window 2:
    #-------------
    DownW2 = cbind(Peaks_Info_x$CIQ.Down.start - Peaks_Info_x$CIQ.Down.size * windows[2]/2, 
        Peaks_Info_x$CIQ.Down.start + Peaks_Info_x$CIQ.Down.size * windows[2]/2)
    DownW2 = intervals::Intervals(DownW2, closed = c(TRUE, TRUE))
    # observed pets in intervals
    DownW2_obs = lengths(intervals::interval_included(DownW2, DTag_x))
    lambdaDownW2 = DownW2_obs * Peaks_Info_x$CIQ.Down.size/(intervals::size(DownW2) + 
        1)
    #-------------
    # ------ genome backgrounds:
    #-------------
    lambdaBGC_Up = ChromPETs * Peaks_Info_x$CIQ.Up.size/ChromSize
    lambdaBGC_Down = ChromPETs * Peaks_Info_x$CIQ.Down.size/ChromSize
    #-------------
    # ------ find maximum:
    #-------------
    OptimalLambda = plyr::llply(seq_len(nrow(Peaks_Info_x)), function(i, lambdaUpW1, 
        lambdaUpW2, lambdaBGC_Up, lambdaDownW1, lambdaDownW2, lambdaBGC_Down) {
        UpMax = max(c(2, lambdaUpW1[i], lambdaUpW2[i], lambdaBGC_Up[i]))
        DownMax = max(c(2, lambdaDownW1[i], lambdaDownW2[i], lambdaBGC_Down[i]))
        return(c(UpMax, DownMax))
    }, lambdaUpW1 = lambdaUpW1, lambdaUpW2 = lambdaUpW2, lambdaBGC_Up = lambdaBGC_Up, 
        lambdaDownW1 = lambdaDownW1, lambdaDownW2 = lambdaDownW2, lambdaBGC_Down = lambdaBGC_Down)
    OptimalLambda = do.call(rbind, OptimalLambda)
    #-------------
    # ------ Find folds and p.values:
    #-------------
    # Upstream
    FoldEnrichUp = Peaks_Info_x$Pets/OptimalLambda[, 1]
    pvalueUp = stats::ppois(q = Peaks_Info_x$Pets, lambda = OptimalLambda[, 1], lower.tail = FALSE)
    # Downstream
    FoldEnrichDown = Peaks_Info_x$Pets/OptimalLambda[, 2]
    pvalueDown = stats::ppois(q = Peaks_Info_x$Pets, lambda = OptimalLambda[, 2], 
        lower.tail = FALSE)
    # merged p-value:
    pvalue = pvalueUp * pvalueDown
    #-------------
    # ------Return
    #-------------
    Peaks_Info_x$lambdaUp = OptimalLambda[, 1]
    Peaks_Info_x$FoldEnrichUp = FoldEnrichUp
    Peaks_Info_x$p.valueUp = pvalueUp
    Peaks_Info_x$lambdaDown = OptimalLambda[, 2]
    Peaks_Info_x$FoldEnrichDown = FoldEnrichDown
    Peaks_Info_x$p.valueDown = pvalueDown
    Peaks_Info_x$p.value = pvalue
    return(Peaks_Info_x)
}
#-------------
#-------------
# done Functions for plotting:
#-------------
#-------------
Get_image_S3_P1_fun = function(S3_AnalysisDir, SA_prefix, PeakInfo) {
    # Rcheck:
    CIQ.Up.size = NULL
    CIQ.Down.size = NULL
    FDR0.05 = NULL
    # image dir:
    S3_P1_image_dir = file.path(S3_AnalysisDir, paste(SA_prefix, "_stage_3_p1_image.jpg", 
        sep = ""))
    # plot
    S3_image_p1 = ggplot2::ggplot(PeakInfo, ggplot2::aes(CIQ.Up.size, CIQ.Down.size, 
        color = factor(FDR0.05))) + ggplot2::geom_point(size = 1) + ggplot2::theme_bw() + 
        ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"), panel.grid.major = ggplot2::element_blank(), 
            panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank()) + 
        ggplot2::theme(axis.text = ggplot2::element_text(size = 15, color = "black"), 
            axis.title = ggplot2::element_text(size = 18, color = "black"), plot.title = ggplot2::element_text(size = 18, 
                color = "black"), legend.text = ggplot2::element_text(size = 18), 
            legend.title = ggplot2::element_text(size = 18)) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
        ggplot2::scale_color_discrete(guide = ggplot2::guide_legend(title = "Binding Site FDR")) + 
        ggplot2::scale_y_continuous(expand = c(0, 0)) + ggplot2::scale_x_continuous(expand = c(0, 
        0)) + ggplot2::xlab("Upstream Peak Size") + ggplot2::ylab("Downstream Peak Size") + 
        ggplot2::geom_abline(intercept = 0, slope = 1) + ggplot2::ggtitle("Upstream/Downstream peak sizes")
    # save:
    ggplot2::ggsave(plot = S3_image_p1, file = S3_P1_image_dir, scale = 2)
}
#-------------
#-------------
Get_image_S3_P2_fun = function(S3_AnalysisDir, SA_prefix, PeakInfo) {
    # Rcheck:
    ID = NULL
    FDR = NULL
    # image dir:
    S3_P2_image_dir = file.path(S3_AnalysisDir, paste(SA_prefix, "_stage_3_p2_image.jpg", 
        sep = ""))
    S3_image_p2 = ggplot2::ggplot(PeakInfo, ggplot2::aes(ID, sort(FDR, decreasing = FALSE))) + 
        ggplot2::geom_line(size = 1) + ggplot2::geom_hline(yintercept = 0.05, color = "red", 
        size = 1) + ggplot2::theme_bw() + ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"), 
        panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), 
        panel.background = ggplot2::element_blank()) + ggplot2::theme(axis.text = ggplot2::element_text(size = 15, 
        color = "black"), axis.title = ggplot2::element_text(size = 18, color = "black"), 
        plot.title = ggplot2::element_text(size = 18, color = "black")) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
        ggplot2::scale_y_continuous(expand = c(0, 0)) + ggplot2::scale_x_continuous(expand = c(0, 
        0)) + ggplot2::xlab("Sorted Binding Site IDs") + ggplot2::ylab("FDR") + ggplot2::ggtitle("FDR values for the binding sites")
    # save:
    ggplot2::ggsave(plot = S3_image_p2, file = S3_P2_image_dir, scale = 2)
}
#-------------
#-------------
Get_image_S3_P3_fun = function(S3_AnalysisDir, SA_prefix, PeakInfo) {
    # Rcheck:
    CIQ.Peak.size = NULL
    FDR0.05 = NULL
    # image dir:
    S3_P3_image_dir = file.path(S3_AnalysisDir, paste(SA_prefix, "_stage_3_p3_image.jpg", 
        sep = ""))
    # plot
    S3_image_p3 = ggplot2::ggplot(PeakInfo, ggplot2::aes(CIQ.Peak.size, fill = factor(FDR0.05))) + 
        ggplot2::geom_density(size = 1, alpha = 0.5) + ggplot2::theme_bw() + ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"), 
        panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), 
        panel.background = ggplot2::element_blank()) + ggplot2::theme(axis.text = ggplot2::element_text(size = 15, 
        color = "black"), axis.title = ggplot2::element_text(size = 18, color = "black"), 
        plot.title = ggplot2::element_text(size = 18, color = "black"), legend.text = ggplot2::element_text(size = 18), 
        legend.title = ggplot2::element_blank()) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
        ggplot2::scale_y_continuous(expand = c(0, 0), labels = function(y) format(y, 
            scientific = FALSE)) + ggplot2::scale_x_continuous(expand = c(0, 0)) + 
        ggplot2::xlab("Binding Site Sizes") + ggplot2::ylab("Density") + ggplot2::geom_abline(intercept = 0, 
        slope = 1) + ggplot2::ggtitle("Densities for the binding site sizes")
    # save:
    ggplot2::ggsave(plot = S3_image_p3, file = S3_P3_image_dir, scale = 2)
}
#-------------
#-------------
Get_image_S3_P4_fun = function(S3_AnalysisDir, SA_prefix, PeakInfo) {
    # Rcheck:
    FDRUp = NULL
    FDRDown = NULL
    FDR0.05 = NULL
    # image dir:
    S3_P4_image_dir = file.path(S3_AnalysisDir, paste(SA_prefix, "_stage_3_p4_image.jpg", 
        sep = ""))
    # plot
    S3_image_p4 = ggplot2::ggplot(PeakInfo, ggplot2::aes(FDRUp, FDRDown, color = factor(FDR0.05))) + 
        ggplot2::geom_point(size = 1) + ggplot2::theme_bw() + ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"), 
        panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), 
        panel.background = ggplot2::element_blank()) + ggplot2::theme(axis.text = ggplot2::element_text(size = 15, 
        color = "black"), axis.title = ggplot2::element_text(size = 18, color = "black"), 
        plot.title = ggplot2::element_text(size = 18, color = "black"), legend.text = ggplot2::element_text(size = 18), 
        legend.title = ggplot2::element_text(size = 18)) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
        ggplot2::scale_color_discrete(guide = ggplot2::guide_legend(title = "Binding Site FDR")) + 
        ggplot2::scale_y_continuous(expand = c(0, 0)) + ggplot2::scale_x_continuous(expand = c(0, 
        0)) + ggplot2::xlab("FDR for Upstream Binding Site Peak") + ggplot2::ylab("FDR for Downstream Binding Site Peak") + 
        ggplot2::geom_abline(intercept = 0, slope = 1) + ggplot2::ggtitle("Upstream/Downstream peak FDR")
    # save:
    ggplot2::ggsave(plot = S3_image_p4, file = S3_P4_image_dir, scale = 2)
}
