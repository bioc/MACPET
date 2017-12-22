#' @importFrom ShortRead FastqStreamer yield writeFastq
#' @importFrom Biostrings DNAStringSet narrow width BStringSet
#' @importClassesFrom ShortRead FastqStreamer ShortReadQ SFastqQuality SRFilterResult
#' @importClassesFrom Biostrings DNAStringSet BStringSet
############################################## Main functions for stage 0
Stage_0_Main_fun = function(SA_prefix, S0_fastq1, S0_fastq2, S0_LinkerA, S0_LinkerB, 
    S0_MinReadLength, S0_MaxReadLength, S0_LinkerOccurence, S0_image, S0_fastqStream, 
    S0_Totfastqreads, S0_AnalysisDir) {
    # Take time:
    Analysis.time.start = Sys.time()
    #----------------
    # create names to save the resulted fastq:
    #----------------
    if (!dir.exists(S0_AnalysisDir)) 
        dir.create(S0_AnalysisDir)
    FastqWriteList = Create_fastq_list_fun(S0_AnalysisDir = S0_AnalysisDir, SA_prefix = SA_prefix)
    # image dir:
    S0_image_dir = file.path(S0_AnalysisDir, paste(SA_prefix, "_stage_0_image.jpg", 
        sep = ""))
    #----------------
    # break and filter fastq files:
    #----------------
    cat("Filtering fastq files...\n")
    Parse_fastqfiles_main_fun(S0_fastq1 = S0_fastq1, S0_fastq2 = S0_fastq2, S0_LinkerA = S0_LinkerA, 
        S0_LinkerB = S0_LinkerB, S0_MinReadLength = S0_MinReadLength, S0_MaxReadLength = S0_MaxReadLength, 
        S0_LinkerOccurence = S0_LinkerOccurence, S0_fastqStream = S0_fastqStream, 
        S0_Totfastqreads = S0_Totfastqreads, FastqWriteList = FastqWriteList, S0_image = S0_image, 
        SA_prefix = SA_prefix, S0_image_dir = S0_image_dir)
    # print:
    futile.logger::flog.info("=====================================", name = "SA_LogFile", 
        capture = FALSE)
    futile.logger::flog.info("Stage 0 is done!", name = "SA_LogFile", capture = FALSE)
    futile.logger::flog.info(paste("Analysis results for stage 0 are in:\n", S0_AnalysisDir), 
        name = "SA_LogFile", capture = FALSE)
    # save time:
    Analysis.time.end = Sys.time()
    Total.Time = Analysis.time.end - Analysis.time.start
    LogFile = paste("Total stage 0 time:", Total.Time, " ", units(Total.Time))
    futile.logger::flog.info(LogFile, name = "SA_LogFile", capture = FALSE)
}
# done
#-------------
#-------------
# function for creating the fastq files for the ouput to be saved:
Create_fastq_list_fun = function(S0_AnalysisDir, SA_prefix) {
    #-------------
    # those with A/A or B/B linkers:(used in subsequent stages)
    #-------------
    fastq1_usable_dir = file.path(S0_AnalysisDir, paste(SA_prefix, "_usable_1.fastq.gz", 
        sep = ""))
    if (file.exists(fastq1_usable_dir)) 
        unlink(x = fastq1_usable_dir, recursive = TRUE, force = TRUE)
    fastq2_usable_dir = file.path(S0_AnalysisDir, paste(SA_prefix, "_usable_2.fastq.gz", 
        sep = ""))
    if (file.exists(fastq2_usable_dir)) 
        unlink(x = fastq2_usable_dir, recursive = TRUE, force = TRUE)
    #-------------
    # those with A/B or B/A or one/two parts with no linkers
    #-------------
    fastq1_chimeric_dir = file.path(S0_AnalysisDir, paste(SA_prefix, "_chimeric_1.fastq.gz", 
        sep = ""))
    if (file.exists(fastq1_chimeric_dir)) 
        unlink(x = fastq1_chimeric_dir, recursive = TRUE, force = TRUE)
    fastq2_chimeric_dir = file.path(S0_AnalysisDir, paste(SA_prefix, "_chimeric_2.fastq.gz", 
        sep = ""))
    if (file.exists(fastq2_chimeric_dir)) 
        unlink(x = fastq2_chimeric_dir, recursive = TRUE, force = TRUE)
    #-------------
    # those with lack of linkers in each read sequence. (always discarded,except if
    # S0_KeepEmptyReads)
    #-------------
    fastq1_ambiguous_dir = file.path(S0_AnalysisDir, paste(SA_prefix, "_ambiguous_1.fastq.gz", 
        sep = ""))
    if (file.exists(fastq1_ambiguous_dir)) 
        unlink(x = fastq1_ambiguous_dir, recursive = TRUE, force = TRUE)
    fastq2_ambiguous_dir = file.path(S0_AnalysisDir, paste(SA_prefix, "_ambiguous_2.fastq.gz", 
        sep = ""))
    if (file.exists(fastq2_ambiguous_dir)) 
        unlink(x = fastq2_ambiguous_dir, recursive = TRUE, force = TRUE)
    #-------------
    # add to list:
    #-------------
    FastqWriteList = list()
    FastqWriteList[[1]] = list(object = NA, file = fastq1_usable_dir)
    FastqWriteList[[2]] = list(object = NA, file = fastq2_usable_dir)
    FastqWriteList[[3]] = list(object = NA, file = fastq1_chimeric_dir)
    FastqWriteList[[4]] = list(object = NA, file = fastq2_chimeric_dir)
    FastqWriteList[[5]] = list(object = NA, file = fastq1_ambiguous_dir)
    FastqWriteList[[6]] = list(object = NA, file = fastq2_ambiguous_dir)
    return(FastqWriteList)
}
# done
#-------------
#-------------
# function for cutting linkers from the yield fastq files
Parse_fastqfiles_main_fun = function(S0_fastq1, S0_fastq2, S0_LinkerA, S0_LinkerB, 
    S0_MinReadLength, S0_MaxReadLength, S0_LinkerOccurence, S0_fastqStream, S0_Totfastqreads, 
    FastqWriteList, S0_image, SA_prefix, S0_image_dir) {
    #----------------
    # Stream the fastq files:
    #----------------
    Streamfastq1 = ShortRead::FastqStreamer(con = S0_fastq1, n = S0_fastqStream)
    Streamfastq2 = ShortRead::FastqStreamer(con = S0_fastq2, n = S0_fastqStream)
    # the first of loop you need to create the files, then change to append
    WrittingMode = "w"
    # count yield/chimeric/usable/ambiguous/N reads:
    TotfastqLinesRead = 0
    Totchimeric = 0
    Totusable = 0
    Totambi = 0
    TotNs = 0
    # the following are for printing progress to the console:
    options(scipen = 999, digits = 2)  #for printing
    cat1 = "||lines scanned:"
    cat2 = "||usable reads:"
    cat3 = "||chimeric reads:"
    cat4 = "||ambiguous reads:"
    cat5 = "||N reads:"
    #----------------
    # loop though the files:
    #----------------
    repeat {
        #----------------
        # read the fastq files:
        #----------------
        # read yield for fastq1
        fastq1yield = ShortRead::yield(Streamfastq1)
        # read yield for fastq2
        fastq2yield = ShortRead::yield(Streamfastq2)
        # break if empty, all is read
        Curfastqyieldsize = length(fastq1yield)
        if (Curfastqyieldsize == 0) 
            break
        TotfastqLinesRead = TotfastqLinesRead + Curfastqyieldsize
        #----------------
        # check if fastq files have the same ids and remove /1,/2.
        #----------------
        SortedReads = Check_If_sorted_fun(fastq1yield = fastq1yield, fastq2yield = fastq2yield)
        if (!SortedReads) {
            stop("S0_fastq1 and S0_fastq2 files are not sorted by ID, or their IDs before /1 and /2 are not identical.", 
                call. = FALSE)
        }
        #----------------
        # Find fastq linker and split (c++11):
        #----------------
        FilteringResults = FilterFastqYield_fun_Rcpp(Curfastqyieldsize, as.character(fastq1yield@sread), 
            Biostrings::width(fastq1yield@sread), as.character(fastq2yield@sread), 
            Biostrings::width(fastq2yield@sread), S0_LinkerA, S0_LinkerB, S0_LinkerOccurence, 
            S0_MinReadLength, S0_MaxReadLength)
        #----------------
        # check the whole Yield is NNs and skip:
        #----------------
        if (FilteringResults$TotNsYield == Curfastqyieldsize) {
            TotNs = TotNs + FilteringResults$TotNsYield
            # print:
            cat(cat1, TotfastqLinesRead, "(", TotfastqLinesRead/S0_Totfastqreads * 
                100, "%)|| ", cat2, Totusable, "(", Totusable/TotfastqLinesRead * 
                100, "%)|| ", cat3, Totchimeric, "(", Totchimeric/TotfastqLinesRead * 
                100, "%)|| ", cat4, Totambi, "(", Totambi/TotfastqLinesRead * 100, 
                "%)|| ", cat5, TotNs, "(", TotNs/TotfastqLinesRead * 100, "%)||\r", 
                sep = "")
            # go to next
            next
        }
        #----------------
        # narrow results
        #----------------
        fastq1yield@sread = Biostrings::narrow(fastq1yield@sread, start = 1, end = FilteringResults$NarrowingPos[, 
            1])
        fastq1yield@quality = Biostrings::narrow(fastq1yield@quality, start = 1, 
            end = FilteringResults$NarrowingPos[, 1])
        fastq2yield@sread = Biostrings::narrow(fastq2yield@sread, start = 1, end = FilteringResults$NarrowingPos[, 
            2])
        fastq2yield@quality = Biostrings::narrow(fastq2yield@quality, start = 1, 
            end = FilteringResults$NarrowingPos[, 2])
        #----------------
        # add counts:
        #----------------
        Totusable = Totusable + FilteringResults$TotusableYield
        Totchimeric = Totchimeric + FilteringResults$TotchimericYield
        Totambi = Totambi + FilteringResults$TotambiYield
        TotNs = TotNs + FilteringResults$TotNsYield
        #----------------
        # split chimeric, usable and ambiguous for the fastq files:
        #----------------
        # take Ids:
        Usable_id = which(FilteringResults$FilterClasses == 1)
        Chim_id = which(FilteringResults$FilterClasses == 2)
        Amb_id = which(FilteringResults$FilterClasses == 0)
        # split:
        FastqWriteList[[1]]$object = fastq1yield[Usable_id]
        FastqWriteList[[2]]$object = fastq2yield[Usable_id]
        FastqWriteList[[3]]$object = fastq1yield[Chim_id]
        FastqWriteList[[4]]$object = fastq2yield[Chim_id]
        FastqWriteList[[5]]$object = fastq1yield[Amb_id]
        FastqWriteList[[6]]$object = fastq2yield[Amb_id]
        #----------------
        # save data: fastq1yield_usable, fastq2yield_usable, fastq1yield_chimeric,
        # fastq2yield_chimeric
        #----------------
        for (i in seq_len(6)) {
            ShortRead::writeFastq(object = FastqWriteList[[i]]$object, file = FastqWriteList[[i]]$file, 
                mode = WrittingMode, compress = TRUE)
        }
        # change mode to append now:
        WrittingMode = "a"
        # print:
        cat(cat1, TotfastqLinesRead, "(", TotfastqLinesRead/S0_Totfastqreads * 100, 
            "%)|| ", cat2, Totusable, "(", Totusable/TotfastqLinesRead * 100, "%)|| ", 
            cat3, Totchimeric, "(", Totchimeric/TotfastqLinesRead * 100, "%)|| ", 
            cat4, Totambi, "(", Totambi/TotfastqLinesRead * 100, "%)|| ", cat5, TotNs, 
            "(", TotNs/TotfastqLinesRead * 100, "%)||\r", sep = "")
    }
    # close connections:
    close(Streamfastq1)
    close(Streamfastq2)
    #----------------
    # write in log file:
    #----------------
    TotfastqLinesRead100 = TotfastqLinesRead/S0_Totfastqreads * 100
    Totusable100 = Totusable/TotfastqLinesRead * 100
    Totchimeric100 = Totchimeric/TotfastqLinesRead * 100
    Totambi100 = Totambi/TotfastqLinesRead * 100
    TotNs100 = TotNs/TotfastqLinesRead * 100
    LogFile = list()
    LogFile[[1]] = paste("Total lines processed:", TotfastqLinesRead, "(", TotfastqLinesRead100, 
        "%)")
    LogFile[[2]] = paste("Total usable reads:", Totusable, "(", Totusable100, "%)")
    LogFile[[3]] = paste("Total chimeric reads:", Totchimeric, "(", Totchimeric100, 
        "%)")
    LogFile[[4]] = paste("Total ambiguous reads:", Totambi, "(", Totambi100, "%)")
    LogFile[[5]] = paste("Total NNs reads:", TotNs, "(", TotNs100, "%)")
    for (lf in seq_len(5)) futile.logger::flog.info(LogFile[[lf]], name = "SA_LogFile", 
        capture = FALSE)
    #----------------
    # plot image:
    #----------------
    if (S0_image) {
        Get_image_S0_fun(Totusable100 = Totusable100, Totchimeric100 = Totchimeric100, 
            Totambi100 = Totambi100, TotNs100 = TotNs100, S0_image_dir = S0_image_dir)
    }
}
# done
#-------------
#-------------
# function for checking if the fastq yields are sorted/have same id
Check_If_sorted_fun = function(fastq1yield, fastq2yield) {
    # get widths of IDs:
    fastq1_idwidth = Biostrings::width(fastq1yield@id)
    fastq2_idwidth = Biostrings::width(fastq2yield@id)
    # remove the /1 and /2 from the fastq files:
    id1 = Biostrings::narrow(fastq1yield@id, start = 1, end = fastq1_idwidth - 1)
    id1 = as.character(id1)
    id2 = Biostrings::narrow(fastq2yield@id, start = 1, end = fastq2_idwidth - 1)
    id2 = as.character(id2)
    # return:
    return(identical(id1, id2))
}
# done
#-------------
#-------------
# function for plotting for stage 0
Get_image_S0_fun = function(Totusable100, Totchimeric100, Totambi100, TotNs100, S0_image_dir) {
    # Rcheck:
    Value = NULL
    Kind = NULL
    #-------------
    # create data:
    #-------------
    S0_imagedata = data.frame(Kind = c(paste("Usable (", round(Totusable100, digits = 1), 
        "%)", sep = ""), paste("Chimeric (", round(Totchimeric100, digits = 1), "%)", 
        sep = ""), paste("Ambiguous (", round(Totambi100, digits = 1), "%)", sep = ""), 
        paste("NNs (", round(TotNs100, digits = 1), "%)", sep = "")), Value = c(round(Totusable100), 
        round(Totchimeric100), round(Totambi100), round(TotNs100)))
    #-------------
    # plot the split:
    #-------------
    S0_image = ggplot2::ggplot(S0_imagedata, ggplot2::aes(x = "", y = Value, fill = factor(Kind))) + 
        ggplot2::geom_bar(width = 1, stat = "identity") + ggplot2::coord_polar(theta = "y") + 
        ggplot2::theme(axis.title = ggplot2::element_blank(), plot.title = ggplot2::element_text(size = 20, 
            color = "black"), legend.title = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = 17), 
            axis.text = ggplot2::element_blank(), legend.position = "bottom", legend.direction = "vertical", 
            axis.ticks = ggplot2::element_blank()) + ggplot2::ggtitle("Pie chart for fastq files") + 
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_fill_brewer(palette = "Dark2")
    # save:
    ggplot2::ggsave(plot = S0_image, file = S0_image_dir, scale = 2)
}
# done
