#' @importFrom Rbowtie bowtie_build bowtie
#' @importFrom BiocParallel bpworkers bplapply
#' @importFrom GEOquery gunzip
#' @importFrom Rsamtools mergeBam sortBam countBam asBam scanBamFlag ScanBamParam filterBam scanBamWhat BamFile asSam
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom rtracklayer export
#' @importFrom Biostrings BStringSet width intersect
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb genome seqinfo
#' @importFrom methods as
#' @importFrom InteractionSet GInteractions
#' @importFrom rbamtools bamReader readerToFastq
#' @importFrom IRanges narrow
#' @importClassesFrom Rsamtools BamFile
#' @importClassesFrom GenomicAlignments GAlignments
#' @importClassesFrom GenomeInfoDb Seqinfo
#' @importClassesFrom GenomicRanges GRanges
############################################## Main functions for stage 1
#-------------
#-------------
# main function for aligning the reads from the fastq files and creating
# pairedEnd BAM:
Stage_1_Main_fun = function(SA_prefix, S1_fastq1_usable_dir, S1_fastq2_usable_dir,
    S1_image, S1_BAMStream, S1_makeSam, S1_genome, S1_RbowtieIndexBuild, S1_RbowtieIndexDir,
    S1_RbowtieIndexPrefix, S1_RbowtieRefDir, S1_AnalysisDir) {
    # Take time:
    Analysis.time.start = Sys.time()
    # create directory for the align results only:
    if (!dir.exists(S1_AnalysisDir))
        dir.create(S1_AnalysisDir)
    #----------------
    # Create index if needed:
    #----------------
    IndexesPath = BuildBowtieIndex_fun(S1_RbowtieIndexBuild = S1_RbowtieIndexBuild,
        S1_RbowtieRefDir = S1_RbowtieRefDir, S1_RbowtieIndexDir = S1_RbowtieIndexDir,
        S1_RbowtieIndexPrefix = S1_RbowtieIndexPrefix)
    #----------------
    # Map the usable fastq files with zero mismatches:
    #----------------
    SAMout12_M0 = Map_fastq_V0_main_fun(S1_AnalysisDir = S1_AnalysisDir, S1_fastq1_usable_dir = S1_fastq1_usable_dir,
        S1_fastq2_usable_dir = S1_fastq2_usable_dir, IndexesPath = IndexesPath, SA_prefix = SA_prefix)
    #----------------
    # Convert SAM to BAM, filter and convert the unmapped to fastq again for second
    # round
    #----------------
    BAMstats12V0 = Convert_Filter_CreateFastq_main_fun(S1_AnalysisDir = S1_AnalysisDir,
        SA_prefix = SA_prefix, SAMout12_M0 = SAMout12_M0)
    #----------------
    # Map the unmapped with at most one mismatch
    #----------------
    SAMout12_M1 = Map_fastq_V1_main_fun(S1_AnalysisDir = S1_AnalysisDir, IndexesPath = IndexesPath,
        SA_prefix = SA_prefix)
    #----------------
    # Convert SAM to BAM again, filter out the unmapped and merge the BAM files V1,V0
    #----------------
    BAM12 = SAMtoBAM_convert_fun(S1_AnalysisDir = S1_AnalysisDir, SAMout12_M1 = SAMout12_M1,
        SA_prefix = SA_prefix, S1_image = S1_image, BAMstats12V0 = BAMstats12V0)
    #----------------
    # Merge BAM files from the two reads
    #----------------
    MergedBAM = MergeBAMfiles_fun(S1_AnalysisDir = S1_AnalysisDir, BAM12 = BAM12,
        SA_prefix = SA_prefix)
    #----------------
    # sort by Qname
    #----------------
    SortBAMQname_fun(MergedBAM = MergedBAM)
    #----------------
    # fix mates
    #----------------
    PairedEndBAMpath = FixMates_main_fun(MergedBAM = MergedBAM, S1_AnalysisDir = S1_AnalysisDir,
        S1_BAMStream = S1_BAMStream, SA_prefix = SA_prefix, S1_image = S1_image,
        S1_genome = S1_genome, CalledFromConvToPE_BAM = FALSE)
    #----------------
    # If they need the sam files, convert PairedEndBAM to two SAM files:
    #----------------
    if (S1_makeSam) {
        GetSAMFiles_fun(PairedEndBAMpath = PairedEndBAMpath, S1_AnalysisDir = S1_AnalysisDir,
            SA_prefix = SA_prefix)
    }
    #----------------
    # print:
    #----------------
    futile.logger::flog.info("=====================================", name = "SA_LogFile",
        capture = FALSE)
    futile.logger::flog.info("Stage 1 is done!", name = "SA_LogFile", capture = FALSE)
    futile.logger::flog.info(paste("Analysis results for stage 1 are in:\n", S1_AnalysisDir),
        name = "SA_LogFile", capture = FALSE)
    # save time:
    Analysis.time.end = Sys.time()
    Total.Time = Analysis.time.end - Analysis.time.start
    LogFile = paste("Total stage 1 time:", Total.Time, " ", units(Total.Time))
    futile.logger::flog.info(LogFile, name = "SA_LogFile", capture = FALSE)
}
# done
#----------------
#----------------
# function for building the bowtie index iff needed.
BuildBowtieIndex_fun = function(S1_RbowtieIndexBuild, S1_RbowtieRefDir, S1_RbowtieIndexDir,
    S1_RbowtieIndexPrefix) {
    if (S1_RbowtieIndexBuild) {
        # create folder for the index:
        cat("Building bowtie index...")
        if (!dir.exists(S1_RbowtieIndexDir))
            dir.create(S1_RbowtieIndexDir)
        # create the index, C=FALSE for colorspace index
        Rbowtie::bowtie_build(references = S1_RbowtieRefDir, outdir = S1_RbowtieIndexDir,
            prefix = S1_RbowtieIndexPrefix, force = TRUE, strict = TRUE, C = FALSE)
        cat("Done\n")
    }
    # indexes path:
    IndexesPath = file.path(S1_RbowtieIndexDir, S1_RbowtieIndexPrefix)
    return(IndexesPath)
}
# done
#----------------
#----------------
# main function for mapping the fastq files with 0 mismatch
Map_fastq_V0_main_fun = function(S1_AnalysisDir, S1_fastq1_usable_dir, S1_fastq2_usable_dir,
    IndexesPath, SA_prefix) {
    # -------align Notes: c=FALSE not cmd files input. C=FALSE not colorspace, v=0
    # max mismatches, ignore qualities unless ties.  m=1 for keeping unique mapped.
    # k=1 report one hit per read
    #----------------
    # map first reads
    #----------------
    cat("========>Mapping first reads with zero mismatch...\n")
    # output names:
    namefastqgz1 = basename(S1_fastq1_usable_dir)
    namefastq1 = paste(SA_prefix, "_usable_1_MR1.fastq", sep = "")
    namesam1 = paste(SA_prefix, "_usable_1_MR1.sam", sep = "")
    # run:
    SAMout1 = Map_fastq_V0_sub_fun(S1_AnalysisDir = S1_AnalysisDir, fastq_usable_dir = S1_fastq1_usable_dir,
        namefastqgz = namefastqgz1, namefastq = namefastq1, namesam = namesam1, IndexesPath = IndexesPath)
    cat("Done\n")
    #----------------
    # map second reads
    #----------------
    cat("========>Mapping second reads with zero mismatch...\n")
    # output names:
    namefastqgz2 = basename(S1_fastq2_usable_dir)
    namefastq2 = paste(SA_prefix, "_usable_2_MR1.fastq", sep = "")
    namesam2 = paste(SA_prefix, "_usable_2_MR1.sam", sep = "")
    # run:
    SAMout2 = Map_fastq_V0_sub_fun(S1_AnalysisDir = S1_AnalysisDir, fastq_usable_dir = S1_fastq2_usable_dir,
        namefastqgz = namefastqgz2, namefastq = namefastq2, namesam = namesam2, IndexesPath = IndexesPath)
    cat("Done\n")
    return(list(SAMout1 = SAMout1, SAMout2 = SAMout2))
}
# done
#----------------
#----------------
# function for mapping each fastq with 0 mismatch:
Map_fastq_V0_sub_fun = function(S1_AnalysisDir, fastq_usable_dir, namefastqgz, namefastq,
    namesam, IndexesPath) {
    #----------------
    # unzip usable fastq
    #----------------
    cat("Unziping reads for mapping ", namefastqgz, "...", sep = "")
    fastq_usable_dir_ungz = file.path(S1_AnalysisDir, namefastq)
    suppressMessages(GEOquery::gunzip(filename = fastq_usable_dir, overwrite = TRUE,
        remove = FALSE, destname = fastq_usable_dir_ungz))
    cat("OK\n")
    #----------------
    # Map usable fastq
    #----------------
    cat("Mapping reads in ", namefastq, "...", sep = "")
    SAMout = file.path(S1_AnalysisDir, namesam)
    Rbowtie::bowtie(sequences = fastq_usable_dir_ungz, index = IndexesPath, type = "single",
        outfile = SAMout, force = TRUE, strict = TRUE, c = FALSE, C = FALSE, trim5 = 0,
        trim3 = 0, quiet = TRUE, sam = TRUE, threads = BiocParallel::bpworkers(),
        verbose = FALSE, chunkmbs = 500, v = 0, k = 1, m = 1)
    cat("Done\n")
    #----------------
    # remove unziped fastq
    #----------------
    cat("Removing unnecessary files...\n")
    unlink(x = fastq_usable_dir_ungz, recursive = TRUE, force = TRUE)
    return(SAMout)
}
#----------------
#----------------
# main function for converting to bam, filtering, creating fastq for the
# unmapped.
Convert_Filter_CreateFastq_main_fun = function(S1_AnalysisDir, SA_prefix, SAMout12_M0) {
    cat("==================================================\n")
    cat("Preparing files for mapping with at most one mismatch...")
    #----------------
    # Create list of inputs:
    #----------------
    Converting = list()
    Converting[[1]] = list(SAMoutMR1 = SAMout12_M0[[1]], BAMoutMR1 = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_1_MR1", sep = "")), BAMoutbaiMR1 = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_1_MR1.bam.bai", sep = "")), BAMoutV0 = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_1_V0.bam", sep = "")), BAMoutMR2 = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_1_MR2.bam", sep = "")), BAMoutbaiMR2 = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_1_MR2.bam.bai", sep = "")), FastqMR2 = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_1_MR2.fastq.gz", sep = "")))
    Converting[[2]] = list(SAMoutMR1 = SAMout12_M0[[2]], BAMoutMR1 = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_2_MR1", sep = "")), BAMoutbaiMR1 = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_2_MR1.bam.bai", sep = "")), BAMoutV0 = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_2_V0.bam", sep = "")), BAMoutMR2 = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_2_MR2.bam", sep = "")), BAMoutbaiMR2 = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_2_MR2.bam.bai", sep = "")), FastqMR2 = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_2_MR2.fastq.gz", sep = "")))
    # run in parallel:
    BAMstats12V0 = BiocParallel::bplapply(X = Converting, FUN = Convert_Filter_CreateFastq_sub_fun)
    cat("Done\n")
    return(BAMstats12V0)
}
# DOne
#----------------
#----------------
# function for splitting and converting files for mapping round 2:
Convert_Filter_CreateFastq_sub_fun = function(ConvertingL) {
    #----------------
    # Convert to BAM:
    #----------------
    # suppress warnings for merging
    suppressWarnings(Rsamtools::asBam(file = ConvertingL$SAMoutMR1, destination = ConvertingL$BAMoutMR1,
        overwrite = TRUE, indexDestination = TRUE))
    ConvertingL$BAMoutMR1 = paste(ConvertingL$BAMoutMR1, ".bam", sep = "")
    # Get the total counts:
    TotReads = Rsamtools::countBam(ConvertingL$BAMoutMR1)$records
    #----------------
    # Filter the bam file to mapped with zero mismatch and unmapped:
    #----------------
    #-------------------------filter mapped:
    # make flag for keeping the mapped only:
    MappedFlag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
        isNotPassingQualityControls = FALSE, isDuplicate = FALSE)
    # make ScanBamParam:
    SBparam = Rsamtools::ScanBamParam(flag = MappedFlag)
    # filter, dont need index since they will be merged:
    Rsamtools::filterBam(file = ConvertingL$BAMoutMR1, index = ConvertingL$BAMoutbaiMR1,
        destination = ConvertingL$BAMoutV0, param = SBparam, indexDestination = FALSE)
    TotmappedV0 = Rsamtools::countBam(ConvertingL$BAMoutV0)$records
    #-------------------------filter umapped:
    # make flag for keeping the unmapped only:
    UnMappedFlag = Rsamtools::scanBamFlag(isUnmappedQuery = TRUE, isSecondaryAlignment = FALSE,
        isNotPassingQualityControls = FALSE, isDuplicate = FALSE)
    # make ScanBamParam:
    SBparam = Rsamtools::ScanBamParam(flag = UnMappedFlag)
    # filter, dont need index since they will be merged:
    Rsamtools::filterBam(file = ConvertingL$BAMoutMR1, index = ConvertingL$BAMoutbaiMR1,
        destination = ConvertingL$BAMoutMR2, param = SBparam, indexDestination = TRUE)
    #----------------
    # Convert the umapped to fastq:
    #----------------
    BAMreader = rbamtools::bamReader(filename = ConvertingL$BAMoutMR2, indexname = ConvertingL$BAMoutbaiMR2)
    rbamtools::readerToFastq(object = BAMreader, filename = ConvertingL$FastqMR2)
    #----------------
    # delete unnecessary files and return:
    #----------------
    unlink(x = c(ConvertingL$SAMoutMR1, ConvertingL$BAMoutMR1, ConvertingL$BAMoutbaiMR1,
        ConvertingL$BAMoutMR2, ConvertingL$BAMoutbaiMR2), recursive = TRUE, force = TRUE)
    return(list(TotReads = TotReads, TotmappedV0 = TotmappedV0))
}
#----------------
#----------------
# main function for mapping the fastq files with 1 mismatch
Map_fastq_V1_main_fun = function(S1_AnalysisDir, IndexesPath, SA_prefix) {
    # -------align Notes: c=FALSE not cmd files input. C=FALSE not colorspace, v=1
    # max mismatches, ignore qualities unless ties.  m=1 for keeping unique mapped.
    # k=1 report one hit per read
    #----------------
    # map first reads
    #----------------
    cat("========>Mapping first reads with one mismatch...\n")
    # output names:
    S1_fastq1_usable_dir = file.path(S1_AnalysisDir, paste(SA_prefix, "_usable_1_MR2.fastq.gz",
        sep = ""))
    namefastqgz1 = basename(S1_fastq1_usable_dir)
    namefastq1 = paste(SA_prefix, "_usable_1_MR2.fastq", sep = "")
    namesam1 = paste(SA_prefix, "_usable_1_MR2.sam", sep = "")
    # run:
    SAMout1 = Map_fastq_V1_sub_fun(S1_AnalysisDir = S1_AnalysisDir, fastq_usable_dir = S1_fastq1_usable_dir,
        namefastqgz = namefastqgz1, namefastq = namefastq1, namesam = namesam1, IndexesPath = IndexesPath)
    cat("Done\n")
    #----------------
    # map second reads
    #----------------
    cat("========>Mapping second reads with one mismatch...\n")
    # output names:
    S1_fastq2_usable_dir = file.path(S1_AnalysisDir, paste(SA_prefix, "_usable_2_MR2.fastq.gz",
        sep = ""))
    namefastqgz2 = basename(S1_fastq2_usable_dir)
    namefastq2 = paste(SA_prefix, "_usable_2_MR2.fastq", sep = "")
    namesam2 = paste(SA_prefix, "_usable_2_MR2.sam", sep = "")
    # run:
    SAMout2 = Map_fastq_V1_sub_fun(S1_AnalysisDir = S1_AnalysisDir, fastq_usable_dir = S1_fastq2_usable_dir,
        namefastqgz = namefastqgz2, namefastq = namefastq2, namesam = namesam2, IndexesPath = IndexesPath)
    cat("Done\n")
    return(list(SAMout1 = SAMout1, SAMout2 = SAMout2))
}
# done
#----------------
#----------------
# function for mapping each fastq with 1 mismatch:
Map_fastq_V1_sub_fun = function(S1_AnalysisDir, fastq_usable_dir, namefastqgz, namefastq,
    namesam, IndexesPath) {
    #----------------
    # unzip usable fastq
    #----------------
    cat("Unziping reads for mapping ", namefastqgz, "...", sep = "")
    fastq_usable_dir_ungz = file.path(S1_AnalysisDir, namefastq)
    suppressMessages(GEOquery::gunzip(filename = fastq_usable_dir, overwrite = TRUE,
        remove = FALSE, destname = fastq_usable_dir_ungz))
    cat("OK\n")
    #----------------
    # Map usable fastq
    #----------------
    cat("Mapping reads in ", namefastq, "...", sep = "")
    SAMout = file.path(S1_AnalysisDir, namesam)
    Rbowtie::bowtie(sequences = fastq_usable_dir_ungz, index = IndexesPath, type = "single",
        outfile = SAMout, force = TRUE, strict = TRUE, c = FALSE, C = FALSE, trim5 = 0,
        trim3 = 0, quiet = TRUE, sam = TRUE, threads = BiocParallel::bpworkers(),
        verbose = FALSE, chunkmbs = 500, v = 1, k = 1, m = 1)
    cat("Done\n")
    #----------------
    # remove unziped fastq
    #----------------
    cat("Removing unnecessary files...\n")
    unlink(x = c(fastq_usable_dir_ungz, fastq_usable_dir), recursive = TRUE, force = TRUE)
    return(SAMout)
}
# Done
#----------------
#----------------
# function to convert the two sam files to bam for V1, and merging the BAM v1 and
# v2
SAMtoBAM_convert_fun = function(S1_AnalysisDir, SAMout12_M1, SA_prefix, S1_image,
    BAMstats12V0) {
    cat("==================================================\n")
    cat("Preparing files for merging...")
    # output/input files:
    bpconvert = list()
    bpconvert[[1]] = list(SAMin = SAMout12_M1$SAMout1, BAMMR2 = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_1_MR2", sep = "")), BAMMR2bai = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_1_MR2.bam.bai", sep = "")), BAMout = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_1.bam", sep = "")), BAMV0 = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_1_V0.bam", sep = "")), BAMV1 = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_1_V1.bam", sep = "")))
    bpconvert[[2]] = list(SAMin = SAMout12_M1$SAMout2, BAMMR2 = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_2_MR2", sep = "")), BAMMR2bai = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_2_MR2.bam.bai", sep = "")), BAMout = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_2.bam", sep = "")), BAMV0 = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_2_V0.bam", sep = "")), BAMV1 = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_2_V1.bam", sep = "")))
    #---------
    # convert to BAM, merge and filter etc:
    #---------
    BAM12res = BiocParallel::bplapply(X = bpconvert, FUN = function(y) {
        # suppress warnings for merging
        suppressWarnings(Rsamtools::asBam(file = y$SAMin, destination = y$BAMMR2,
            overwrite = TRUE, indexDestination = TRUE))
        # update file names:
        y$BAMMR2 = paste(y$BAMMR2, ".bam", sep = "")
        #----------------
        # Filter the bam file to mapped with 1 mismatch
        #----------------
        # make flag for keeping the mapped only:
        MappedFlag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
            isNotPassingQualityControls = FALSE, isDuplicate = FALSE)
        # make ScanBamParam:
        SBparam = Rsamtools::ScanBamParam(flag = MappedFlag)
        # filter, dont need index since they will be merged:
        Rsamtools::filterBam(file = y$BAMMR2, index = y$BAMMR2bai, destination = y$BAMV1,
            param = SBparam, indexDestination = FALSE)
        TotmappedV1 = Rsamtools::countBam(y$BAMV1)$records
        #----------------
        # merge the V1 and V0 files:
        #----------------
        suppressWarnings(Rsamtools::mergeBam(files = c(y$BAMV0, y$BAMV1), destination = y$BAMout,
            overwrite = TRUE, byQname = FALSE, indexDestination = FALSE))
        # unlink:
        unlink(x = c(y$SAMin, y$BAMMR2, y$BAMMR2bai, y$BAMV0, y$BAMV1), recursive = TRUE,
            force = TRUE)
        return(list(TotmappedV1 = TotmappedV1, BAMout = y$BAMout))
    })
    cat("Done\n")
    #---------
    # print and save statistics:
    #---------
    # calculate the stats:
    TotReads_1 = BAMstats12V0[[1]]$TotReads
    TotmappedV0_1 = BAMstats12V0[[1]]$TotmappedV0
    TotmappedV0100_1 = TotmappedV0_1/TotReads_1 * 100
    TotmappedV1_1 = BAM12res[[1]]$TotmappedV1
    TotmappedV1100_1 = TotmappedV1_1/TotReads_1 * 100
    Totunmapped_1 = TotReads_1 - TotmappedV0_1 - TotmappedV1_1
    Totunmapped100_1 = Totunmapped_1/TotReads_1 * 100
    TotReads_2 = BAMstats12V0[[2]]$TotReads
    TotmappedV0_2 = BAMstats12V0[[2]]$TotmappedV0
    TotmappedV0100_2 = TotmappedV0_2/TotReads_2 * 100
    TotmappedV1_2 = BAM12res[[2]]$TotmappedV1
    TotmappedV1100_2 = TotmappedV1_2/TotReads_2 * 100
    Totunmapped_2 = TotReads_2 - TotmappedV0_2 - TotmappedV1_2
    Totunmapped100_2 = Totunmapped_2/TotReads_2 * 100
    # print:
    LogFile = list()
    LogFile[[1]] = "=========>Mapping statistics<========"
    LogFile[[2]] = "==>Statistics for usable_1.bam:"
    LogFile[[3]] = paste("Total reads processed:", TotReads_1)
    LogFile[[4]] = paste("Total reads mapped with zero mismatch:", TotmappedV0_1,
        "(", TotmappedV0100_1, "%)")
    LogFile[[5]] = paste("Total reads mapped with one mismatch:", TotmappedV1_1,
        "(", TotmappedV1100_1, "%)")
    LogFile[[6]] = paste("Total reads unmapped:", Totunmapped_1, "(", Totunmapped100_1,
        "%)")
    LogFile[[7]] = "==>Statistics for usable_2.bam:"
    LogFile[[8]] = paste("Total reads processed:", TotReads_2)
    LogFile[[9]] = paste("Total reads mapped with zero mismatch:", TotmappedV0_2,
        "(", TotmappedV0100_2, "%)")
    LogFile[[10]] = paste("Total reads mapped with one mismatch:", TotmappedV1_2,
        "(", TotmappedV1100_2, "%)")
    LogFile[[11]] = paste("Total reads unmapped:", Totunmapped_2, "(", Totunmapped100_2,
        "%)")
    for (lf in seq_len(11)) futile.logger::flog.info(LogFile[[lf]], name = "SA_LogFile",
        capture = FALSE)
    #---------
    # print and save statistics:
    #---------
    if (S1_image) {
        Get_image_S1_P1_fun(S1_AnalysisDir = S1_AnalysisDir, SA_prefix = SA_prefix,
            TotmappedV0100_1 = TotmappedV0100_1, TotmappedV1100_1 = TotmappedV1100_1,
            Totunmapped100_1 = Totunmapped100_1, TotmappedV0100_2 = TotmappedV0100_2,
            TotmappedV1100_2 = TotmappedV1100_2, Totunmapped100_2 = Totunmapped100_2)
    }
    # return:
    BAM12 = list(BAM1 = BAM12res[[1]]$BAMout, BAM2 = BAM12res[[2]]$BAMout)
    return(BAM12)
}
# done
#-------------
#-------------
# function for plotting for stage 1 part 1: mapped/unmapped reads
Get_image_S1_P1_fun = function(S1_AnalysisDir, SA_prefix, TotmappedV0100_1, TotmappedV1100_1,
    Totunmapped100_1, TotmappedV0100_2, TotmappedV1100_2, Totunmapped100_2) {
    # Rcheck:
    Value = NULL
    Kind = NULL
    # image dir:
    S1_P1_image_dir = file.path(S1_AnalysisDir, paste(SA_prefix, "_stage_1_p1_image.jpg",
        sep = ""))
    #-------------
    # create data:
    #-------------
    S1_imagedata_1 = data.frame(Kind = c(paste("Mapped reads with zero mismatch (",
        round(TotmappedV0100_1, digits = 1), "%)", sep = ""), paste("Mapped reads with one mismatch (",
        round(TotmappedV1100_1, digits = 1), "%)", sep = ""), paste("Unmapped reads (",
        round(Totunmapped100_1, digits = 1), "%)", sep = ""), paste("Mapped reads with zero mismatch (",
        round(TotmappedV0100_2, digits = 1), "%)", sep = ""), paste("Mapped reads with one mismatch  (",
        round(TotmappedV1100_2, digits = 1), "%)", sep = ""), paste("Unmapped reads (",
        round(Totunmapped100_2, digits = 1), "%)", sep = "")), Value = c(round(TotmappedV0100_1),
        round(TotmappedV1100_1), round(Totunmapped100_1), round(TotmappedV0100_2),
        round(TotmappedV1100_2), round(Totunmapped100_2)), facet = c(paste(SA_prefix,
        "_usable_1.bam", sep = ""), paste(SA_prefix, "_usable_1.bam", sep = ""),
        paste(SA_prefix, "_usable_1.bam", sep = ""), paste(SA_prefix, "_usable_2.bam",
            sep = ""), paste(SA_prefix, "_usable_2.bam", sep = ""), paste(SA_prefix,
            "_usable_2.bam", sep = "")))
    #-------------
    # plot the split:
    #-------------
    S1_image_p1 = ggplot2::ggplot(S1_imagedata_1, ggplot2::aes(x = "", y = Value,
        fill = factor(Kind))) + ggplot2::geom_bar(width = 1, stat = "identity") +
        ggplot2::facet_wrap(~facet) + ggplot2::coord_polar(theta = "y") + ggplot2::theme(axis.title = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 20, color = "black"), legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 17), axis.text = ggplot2::element_blank(),
        legend.position = "bottom", legend.direction = "vertical", axis.ticks = ggplot2::element_blank()) +
        ggplot2::ggtitle("Pie chart for the mapping occurence of bam files") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::scale_fill_brewer(palette = "Dark2")
    # save:
    ggplot2::ggsave(plot = S1_image_p1, file = S1_P1_image_dir, scale = 2)
}
# done
#----------------
#----------------
# function for merging the two mapped bam files: (in)
MergeBAMfiles_fun = function(S1_AnalysisDir, BAM12, SA_prefix) {
    cat("Merging ", basename(BAM12[[1]]), ", ", basename(BAM12[[2]]), " files...",
        sep = "")
    # output:
    MergedBAMpath = file.path(S1_AnalysisDir, paste(SA_prefix, "_usable_merged.bam",
        sep = ""))
    suppressWarnings(Rsamtools::mergeBam(files = unlist(BAM12), destination = MergedBAMpath,
        overwrite = TRUE, byQname = FALSE, indexDestination = TRUE))
    cat("Done\n")
    cat("Removing unnecessary files...\n")
    unlink(x = unlist(BAM12), recursive = TRUE, force = TRUE)
    # return:
    MergedBAM = list(BAM = MergedBAMpath, BAMbai = paste(MergedBAMpath, ".bai", sep = ""))
    return(MergedBAM)
}
# done
#----------------
#----------------
# function for sorting the merged bam by Qname (in)
SortBAMQname_fun = function(MergedBAM) {
    cat("Sorting ", basename(MergedBAM$BAM), " file by Qname...", sep = "")
    # output:
    SortBAM = unlist(strsplit(MergedBAM$BAM, ".bam"))
    suppressWarnings(Rsamtools::sortBam(file = MergedBAM$BAM, destination = SortBAM,
        byQname = TRUE))
    cat("Done\n")
}
# done
#----------------
#----------------
# main function for fixing mate-pairs (in)
FixMates_main_fun = function(MergedBAM, S1_AnalysisDir, S1_BAMStream, SA_prefix,
    S1_image, S1_genome, CalledFromConvToPE_BAM) {
    #----------------
    # initiate:
    #----------------
    cat("Pairing reads in ", basename(MergedBAM$BAM), " file...\n", sep = "")
    # counts:
    TotReadsScanned = 0
    TotPairsFound = 0
    TotBAMReads = Rsamtools::countBam(MergedBAM$BAM)$records
    SavedLastRead = NULL
    # what to read:
    SBparam = Rsamtools::ScanBamParam(what = Rsamtools::scanBamWhat())
    # open bam:
    BAMmain = Rsamtools::BamFile(file = MergedBAM$BAM, index = MergedBAM$BAMbai,
        yieldSize = S1_BAMStream, obeyQname = TRUE, asMates = FALSE)
    open(BAMmain)
    BAMcounter = 1
    BAM_path_yield_list = list()
    GIntScanned = NULL
    #----------------
    # scan and mate::
    #----------------
    repeat {
        #----------------
        # scan yield:
        #----------------
        BAMyield = GenomicAlignments::readGAlignments(file = BAMmain, use.names = TRUE,
            param = SBparam)
        # break repeat if yield empty:
        if (length(BAMyield) == 0)
            break
        # take yield length:
        Nyield = length(BAMyield)
        TotReadsScanned = TotReadsScanned + Nyield
        # add the read 1 on to iff it exists:
        if (!is.null(SavedLastRead)) {
            BAMyield = c(SavedLastRead, BAMyield)
            # update Nyield:
            Nyield = Nyield + 1
        }
        #----------------
        # remove dulicates because they cause problems:
        #----------------
        DupliRM = which(duplicated(names(BAMyield)))
        if (length(DupliRM) != 0) {
            BAMyield = BAMyield[-DupliRM]
            Nyield = Nyield - length(DupliRM)
        }
        #----------------
        # take names and get suffixes remove and keep the last read if it is first read:
        #----------------
        PairedNamesRes = Get_Paired_names_fun(BAMyield = BAMyield, Nyield = Nyield)
        #----------------
        # check last read and extract:
        #----------------
        if (PairedNamesRes$LastRead1) {
            SavedLastRead = BAMyield[Nyield]
        } else {
            SavedLastRead = NULL
        }
        # check if pairs found:
        if (PairedNamesRes$NyieldPairs == 0) {
            next
            # go to next yield
        } else {
            TotPairsFound = TotPairsFound + PairedNamesRes$NyieldPairs
            #----------------
            # split the yield in read1 and 2
            #----------------
            BAMyieldSplit = Get_Paired_Yields_fun(BAMyield = BAMyield, PairedNamesRes = PairedNamesRes)
            #----------------
            # add flags to the reads 1:
            #----------------
            FlagRes = FixMatesFlags_fun(BAMyieldSplit = BAMyieldSplit, GIntScanned = GIntScanned)
            GIntScanned = FlagRes$GIntScanned
            BAMyieldPaired = FlagRes$BAMyieldPaired
        }
        #----------------
        # save to BAM format:
        #----------------
        # add genome:
        GenomeInfoDb::genome(GenomeInfoDb::seqinfo(BAMyieldPaired)) = S1_genome
        # save:
        BAM_path_yield = file.path(S1_AnalysisDir, paste(SA_prefix, "_Paired_end_",
            BAMcounter, ".bam", sep = ""))
        rtracklayer::export(object = BAMyieldPaired, con = Rsamtools::BamFile(BAM_path_yield),
            format = "bam")
        BAM_path_yield_list[[BAMcounter]] = BAM_path_yield
        BAMcounter = BAMcounter + 1
        # print:
        TotReadsScanned100 = TotReadsScanned/TotBAMReads * 100
        TotPairsFound100 = 2 * TotPairsFound/TotReadsScanned * 100
        cat("||Total lines scanned: ", TotReadsScanned, "(", TotReadsScanned100,
            "%)|| ", "||Total Pairs registered: ", TotPairsFound, "(", TotPairsFound100,
            "% of the scanned lines)||\r", sep = "")
    }
    # close connection:
    close(BAMmain)
    #----------------
    # write and cat:
    #----------------
    LogFile = list()
    LogFile[[1]] = "=========>Pairing statistics<========"
    LogFile[[2]] = paste("Total reads processed: ", TotReadsScanned, "(", TotReadsScanned100,
        "%)")
    LogFile[[3]] = paste("Total pairs registered:", TotPairsFound, "(", TotPairsFound100,
        "% of the scanned lines)")
    for (lf in seq_len(3)) futile.logger::flog.info(LogFile[[lf]], name = "SA_LogFile",
        capture = FALSE)
    #----------------
    # plot:
    #----------------
    if (S1_image) {
        Get_image_S1_P2_fun(S1_AnalysisDir = S1_AnalysisDir, SA_prefix = SA_prefix,
            TotPairsFound100 = TotPairsFound100)
    }
    #----------------
    # Merge BAM files:
    #----------------
    cat("Merging bam files in ", paste(SA_prefix, "_Paired_end.bam", sep = ""), "...")
    PairedEndBAMpath = file.path(S1_AnalysisDir, paste(SA_prefix, "_Paired_end.bam",
        sep = ""))
    FilesToMerge = unlist(BAM_path_yield_list)
    if (length(FilesToMerge) > 1) {
        PairedEndBAMpath = Rsamtools::mergeBam(files = FilesToMerge, destination = PairedEndBAMpath,
            overwrite = TRUE, byQname = FALSE, indexDestination = TRUE)
    } else if (length(FilesToMerge) == 1) {
        # rename it:
        invisible(file.rename(from = FilesToMerge, to = PairedEndBAMpath))
        # bai:
        PairedEndBAMpathbai = paste(PairedEndBAMpath, ".bai", sep = "")
        FilesToMergebai = paste(FilesToMerge, ".bai", sep = "")
        invisible(file.rename(from = FilesToMergebai, to = PairedEndBAMpathbai))
    } else {
        futile.logger::flog.warn(paste("WARNING: No pairs found in:\n", MergedBAM$BAM),
            name = "SA_LogFile", capture = FALSE)
    }
    cat("Done\n")
    #----------------
    # delete those in BAM_path_yield_list and the usable_merged.bam/bai:
    #----------------
    if (!CalledFromConvToPE_BAM) {
        cat("Removing unnecessary files...\n")
        unlink(x = c(unlist(BAM_path_yield_list), unlist(MergedBAM)), recursive = TRUE,
            force = TRUE)
        unlink(x = paste(unlist(BAM_path_yield_list), ".bai", sep = ""), recursive = TRUE,
            force = TRUE)
    }
    # else they are deleted at ConvertToPE_BAM function.
    return(PairedEndBAMpath)
}
# done
#----------------
#----------------
# function for getting the paired-ends of the reads
Get_Paired_names_fun = function(BAMyield, Nyield) {
    #----------------
    # use BStringSet because it is faster:
    #----------------
    Namesyield = Biostrings::BStringSet(x = names(BAMyield), start = 1)
    NamesyieldWidth = Biostrings::width(Namesyield)
    NamesyieldPref = IRanges::narrow(Namesyield, start = 1, end = NamesyieldWidth -
        2)
    NamesyieldSuf = IRanges::narrow(Namesyield, start = NamesyieldWidth, end = NamesyieldWidth)
    NamesyieldSuf = as.character(NamesyieldSuf)
    # get read1 and read2 positions:
    Which_r1 = which(NamesyieldSuf == "1")
    Which_r2 = which(NamesyieldSuf == "2")
    # split prefixes of read 1 and 2:
    NamesyieldPref_r1 = NamesyieldPref[Which_r1]
    NamesyieldPref_r2 = NamesyieldPref[Which_r2]
    #----------------
    # find intersection of the names (those are your paires for the yield): if any of
    # the two above is empty then NamesyieldPref_r12 will be empty
    # length(NamesyieldPref_r12)==0
    NamesyieldPref_r12 = Biostrings::intersect(NamesyieldPref_r1, NamesyieldPref_r2)
    #----------------
    # check if the last element is read1, if it is then keep it for the next yield.
    LastRead1 = FALSE
    if (Nyield %in% Which_r1)
        LastRead1 = TRUE
    # return:
    return(list(NamesyieldPref_r12 = NamesyieldPref_r12, Which_r1 = Which_r1, Which_r2 = Which_r2,
        NamesyieldPref_r1 = NamesyieldPref_r1, NamesyieldPref_r2 = NamesyieldPref_r2,
        LastRead1 = LastRead1, NyieldPairs = length(NamesyieldPref_r12)))
}
# done
#----------------
#----------------
# function for separating the BAMyield into two pairs
Get_Paired_Yields_fun = function(BAMyield, PairedNamesRes) {
    #----------------
    # get names of read 1 in the intersection:
    #----------------
    read1names = which(as.character(PairedNamesRes$NamesyieldPref_r1) %in% as.character(PairedNamesRes$NamesyieldPref_r12))
    read1_id = PairedNamesRes$Which_r1[read1names]
    BAMyield_r1 = BAMyield[read1_id]
    # replace names without the prefix:
    names(BAMyield_r1) = as.character(PairedNamesRes$NamesyieldPref_r1[read1names])
    #----------------
    # get names of read 2 in the intersection:
    #----------------
    read2names = which(as.character(PairedNamesRes$NamesyieldPref_r2) %in% as.character(PairedNamesRes$NamesyieldPref_r12))
    read2_id = PairedNamesRes$Which_r2[read2names]
    BAMyield_r2 = BAMyield[read2_id]
    # replace names without the prefix:
    names(BAMyield_r2) = as.character(PairedNamesRes$NamesyieldPref_r2[read2names])
    # the orders should be identical but check:
    if (!identical(order(names(BAMyield_r1)), order(names(BAMyield_r2)))) {
        stop("Pairing BAM files failed. The names of the two BAM files are not sorted!",
            call. = FALSE)
    }
    # return:
    return(list(BAMyield_r1 = BAMyield_r1, BAMyield_r2 = BAMyield_r2))
}
# done
#----------------
#----------------
# function for fixing the mates flags and returning the final BAMyieldPaired
FixMatesFlags_fun = function(BAMyieldSplit, GIntScanned) {
    # split:
    BAMyield_r1 = BAMyieldSplit$BAMyield_r1
    BAMyield_r2 = BAMyieldSplit$BAMyield_r2
    # the only info I have is their strand. Take the info before adding other flags:
    Read1RevStrand = which(S4Vectors::mcols(BAMyield_r1)$flag == 16)
    Read2RevStrand = which(S4Vectors::mcols(BAMyield_r2)$flag == 16)
    #----------------
    # Add +1 on all flags to make them paired: Add +64 for first in pair and +128 for
    # second in pair: So add 65 on the first and 129 on the second
    #----------------
    S4Vectors::mcols(BAMyield_r1)$flag = S4Vectors::mcols(BAMyield_r1)$flag + 65
    S4Vectors::mcols(BAMyield_r2)$flag = S4Vectors::mcols(BAMyield_r2)$flag + 129
    #----------------
    # Add +32 for the mates in reverse strand:
    #----------------
    if (length(Read1RevStrand) != 0) {
        S4Vectors::mcols(BAMyield_r2)$flag[Read1RevStrand] = S4Vectors::mcols(BAMyield_r2)$flag[Read1RevStrand] +
            32
    }
    if (length(Read2RevStrand) != 0) {
        S4Vectors::mcols(BAMyield_r1)$flag[Read2RevStrand] = S4Vectors::mcols(BAMyield_r1)$flag[Read2RevStrand] +
            32
    }
    #----------------
    # Add mate positions:
    #----------------
    Pos_r1 = S4Vectors::mcols(BAMyield_r1)$pos
    Pos_r2 = S4Vectors::mcols(BAMyield_r2)$pos
    S4Vectors::mcols(BAMyield_r1)$mpos = Pos_r2
    S4Vectors::mcols(BAMyield_r2)$mpos = Pos_r1
    #----------------
    # Add mate chromosome:
    #----------------
    Chr_r1 = S4Vectors::mcols(BAMyield_r1)$rname
    Chr_r2 = S4Vectors::mcols(BAMyield_r2)$rname
    S4Vectors::mcols(BAMyield_r1)$mrnm = Chr_r2
    S4Vectors::mcols(BAMyield_r2)$mrnm = Chr_r1
    #----------------
    # Check dublicates:
    #----------------
    if (is.null(GIntScanned)) {
        # create the first instance:
        GRscanned1 = methods::as(BAMyield_r1, "GRanges")
        S4Vectors::mcols(GRscanned1) = NULL
        names(GRscanned1) = NULL
        GRscanned2 = methods::as(BAMyield_r2, "GRanges")
        S4Vectors::mcols(GRscanned2) = NULL
        names(GRscanned2) = NULL
        GIntScanned = InteractionSet::GInteractions(anchor1 = GRscanned1, anchor2 = GRscanned2)
        NGIntScanned = 0  #to begin with
    } else {
        # create new instance:
        GRscanned1_new = methods::as(BAMyield_r1, "GRanges")
        S4Vectors::mcols(GRscanned1_new) = NULL
        names(GRscanned1_new) = NULL
        GRscanned2_new = methods::as(BAMyield_r2, "GRanges")
        S4Vectors::mcols(GRscanned2_new) = NULL
        names(GRscanned2_new) = NULL
        GIntScanned_new = InteractionSet::GInteractions(anchor1 = GRscanned1_new,
            anchor2 = GRscanned2_new)
        # merge with the old (take length of the old first)
        NGIntScanned = length(GIntScanned)
        GIntScanned = c(GIntScanned, GIntScanned_new)
    }
    # ckeck duplicates, it gives the second instance as duplicated:
    Duplicated = which(duplicated(GIntScanned) == TRUE)
    # add flags:
    if (length(Duplicated) != 0) {
        # get correct positions:
        DuplPos = Duplicated - NGIntScanned
        # add flag:
        S4Vectors::mcols(BAMyield_r1[DuplPos])$flag = S4Vectors::mcols(BAMyield_r1[DuplPos])$flag +
            1024
        S4Vectors::mcols(BAMyield_r2[DuplPos])$flag = S4Vectors::mcols(BAMyield_r2[DuplPos])$flag +
            1024
        # remove duplicates using Duplicated here:
        GIntScanned = GIntScanned[-Duplicated]
    }
    #----------------
    # Merge in one and return:
    #----------------
    BAMyieldPaired = c(BAMyield_r1, BAMyield_r2)
    return(list(BAMyieldPaired = BAMyieldPaired, GIntScanned = GIntScanned))
}
# done
#-------------
#-------------
# function for plotting for stage 1 part 1: paired/unpaired reads:
Get_image_S1_P2_fun = function(S1_AnalysisDir, SA_prefix, TotPairsFound100) {
    # Rcheck:
    Value = NULL
    Kind = NULL
    # image dir:
    S1_P2_image_dir = file.path(S1_AnalysisDir, paste(SA_prefix, "_stage_1_p2_image.jpg",
        sep = ""))
    #-------------
    # create data:
    #-------------
    S1_imagedata_p2 = data.frame(Kind = c(paste("Paired Reads (", round(TotPairsFound100,
        digits = 1), "%)", sep = ""), paste("Unpaired Reads (", round(100 - TotPairsFound100,
        digits = 1), "%)", sep = "")), Value = c(round(TotPairsFound100, digits = 1),
        round(100 - TotPairsFound100, digits = 1)))
    #-------------
    # plot the split:
    #-------------
    S1_image_p2 = ggplot2::ggplot(S1_imagedata_p2, ggplot2::aes(x = "", y = Value,
        fill = factor(Kind))) + ggplot2::geom_bar(width = 1, stat = "identity") +
        ggplot2::coord_polar(theta = "y") + ggplot2::theme(axis.title = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 20, color = "black"), legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 17), axis.text = ggplot2::element_blank(),
        legend.position = "bottom", legend.direction = "vertical", axis.ticks = ggplot2::element_blank()) +
        ggplot2::ggtitle(paste("Pie chart for ", SA_prefix, "_usable_1.sam file",
            sep = "")) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::scale_fill_brewer(palette = "Dark2")
    # save:
    ggplot2::ggsave(plot = S1_image_p2, file = S1_P2_image_dir, scale = 2)
}
# done
#-------------
#-------------
# function for splitting the paired-end bam file into two sam files
GetSAMFiles_fun = function(PairedEndBAMpath, S1_AnalysisDir, SA_prefix) {
    cat("Creating SAM files from paired-end BAM...")
    #----------------
    # make mates list input:
    #----------------
    FilterMATE = list()
    FilterMATE[[1]] = list(PairedBAM = PairedEndBAMpath, PairedBAMbai = paste(PairedEndBAMpath,
        ".bai", sep = ""), firstmate = TRUE, BAMread = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_1", sep = "")), BAMreadIndex = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_1", sep = "")), SAMout = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_1.sam", sep = "")))
    FilterMATE[[2]] = list(PairedBAM = PairedEndBAMpath, PairedBAMbai = paste(PairedEndBAMpath,
        ".bai", sep = ""), firstmate = FALSE, BAMread = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_2", sep = "")), BAMreadIndex = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_2", sep = "")), SAMout = file.path(S1_AnalysisDir,
        paste(SA_prefix, "_usable_2.sam", sep = "")))
    #----------------
    # split:
    #----------------
    BiocParallel::bplapply(X = FilterMATE, FUN = function(y) {
        # make flag:
        MateFlag = Rsamtools::scanBamFlag(isFirstMateRead = y$firstmate, isSecondMateRead = !y$firstmate,
            isDuplicate = FALSE)
        # make ScanBamParam:
        SBparam = Rsamtools::ScanBamParam(flag = MateFlag)
        # filter, dont need index since they will be merged:
        Rsamtools::filterBam(file = y$PairedBAM, index = y$PairedBAMbai, destination = paste(y$BAMread,
            ".bam", sep = ""), param = SBparam, indexDestination = TRUE)
        # sort by Qname, done already
        Rsamtools::sortBam(file = paste(y$BAMread, ".bam", sep = ""), destination = y$BAMread,
            byQname = TRUE)
        # convert to sam:
        Rsamtools::asSam(file = paste(y$BAMread, ".bam", sep = ""), overwrite = TRUE,
            denstination = y$SAMout)
        # delete files:
        unlink(x = paste(y$BAMread, c(".bam", ".bam.bai"), sep = ""), recursive = TRUE,
            force = TRUE)
    })
    cat("Done\n")
}
# done
