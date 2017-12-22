#----------#' make title:
#' @title Convert two BAM files into one paired-end BAM file.
#----------#' author:
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#----------#' make reference on the article:
#' @references
#' Vardaxis I, Drabløs F, Rye M and Lindqvist BH (2018). \emph{Model-based Analysis for ChIA-PET (MACPET)}.
#' To be published.
#'
#----------#' for the current function to be exported on the NAMESPACE:
#' @export
#----------#' make discription:
#' @description Stage 2 in \code{\link{MACPETUlt}} needs a paired-end BAM file to run. This can be
#'created in Stage 1 using the usable_1 and usable_2 fastq files
#'created in Stage 0. However the user might have two single-end
#'BAM files already created but not paired (by filtering with another way than that in Stage 0 or
#'mapping using another algorithm than that in Stage 1)
#' and only needs to run Stages 2 and 3 in \code{\link{MACPETUlt}}.
#'  \code{ConvertToPE_BAM} can be used on the two BAM files for pairing them, and the
#'  resulted paired-end BAM file can then be used in Stage 2 in \code{\link{MACPETUlt}}.
#'
#----------#' make details sections
#'@details The BAM files \code{BAM_file_1} and \code{BAM_file_2} do not need to
#'be sorted, but their Qnames have to end with /1 and /2 respectively.
#'Furthermore, the BAM files have to include the header section.
#----------#'
#' @seealso
#' \code{\link{MACPETUlt}}, \code{\link{SampleChIAPETDataRead_1.bam}},
#' \code{\link{SampleChIAPETDataRead_2.bam}}
#----------#'output
#' @return A paired-end BAM file named \code{SA_prefix_MACPET_Paired_end.bam} and
#' its index, saved in \code{S1_AnalysisDir}.
#'
#----------#' function parameters:
#'@param S1_AnalysisDir The directory where the resulted paired-end BAM file will be saved.
#'@param SA_prefix see \code{\link{MACPETUlt}}.
#'@param S1_BAMStream see \code{\link{MACPETUlt}}.
#'@param S1_image see \code{\link{MACPETUlt}}.
#'@param S1_genome see \code{\link{MACPETUlt}}.
#'@param BAM_file_1 The directory of the BAM file with the first reads. Their Qnames have to end with /1.
#'@param BAM_file_2 The directory of the BAM file with the second reads. Their Qnames have to end with /2.
#'@param S1_makeSam see \code{\link{MACPETUlt}}.
#'
#'@importFrom methods is
#'@importFrom Rsamtools scanBamHeader sortBam mergeBam filterBam
# ----------#'
#' @examples
#' requireNamespace('ggplot2')
#'
#' #Create a temporary forder, or anywhere you want:
#' S1_AnalysisDir=file.path(tempdir(),'MACPETtest')
#' dir.create(S1_AnalysisDir)#where you will save the results
#'
#' #directories of the BAM files:
#' BAM_file_1=system.file('extdata', 'SampleChIAPETDataRead_1.bam', package = 'MACPET')
#' BAM_file_2=system.file('extdata', 'SampleChIAPETDataRead_2.bam', package = 'MACPET')
#' SA_prefix='MACPET'
#'
#' #convert to paired-end BAM:
#' ConvertToPE_BAM(S1_AnalysisDir=S1_AnalysisDir,
#'                 SA_prefix=SA_prefix,
#'                 S1_BAMStream=2000000,
#'                 S1_image=TRUE,
#'                 S1_genome='hg19',
#'                 BAM_file_1=BAM_file_1,
#'                 BAM_file_2=BAM_file_2)
#'
#' #test if the resulted BAM is paired-end:
#' PairedBAM=file.path(S1_AnalysisDir,paste(SA_prefix,'_Paired_end.bam',sep=''))
#' Rsamtools::testPairedEndBam(file = PairedBAM, index = PairedBAM)
#'
#' bamfile = Rsamtools::BamFile(file = PairedBAM,asMates = TRUE)
#' GenomicAlignments::readGAlignmentPairs(file = bamfile,use.names = FALSE,
#'                                        with.which_label = FALSE,
#'                                        strandMode = 1)
#'
#' #-----delete test directory:
#' unlink(S1_AnalysisDir,recursive=TRUE)
#'
#'
#'
ConvertToPE_BAM = function(S1_AnalysisDir = "", SA_prefix = "MACPET", S1_BAMStream = 2e+06, 
    S1_image = TRUE, S1_genome = "hg19", BAM_file_1 = "", BAM_file_2 = "", S1_makeSam = FALSE) {
    cat("Checking inputs...")
    #------------
    # check directory:
    #------------
    if (!methods::is(S1_AnalysisDir, "character")) {
        stop("S1_AnalysisDir:", S1_AnalysisDir, " is not a directory!", call. = FALSE)
    } else if (!dir.exists(S1_AnalysisDir)) {
        stop("S1_AnalysisDir:", S1_AnalysisDir, " directory does not exist!", call. = FALSE)
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
    # check the logical:
    #------------
    if (!methods::is(S1_image, "logical")) {
        stop("S1_image: ", S1_image, " has to be logical!", call. = FALSE)
    } else if (S1_image) {
        if (!requireNamespace("ggplot2", quietly = TRUE)) {
            stop("ggplot2 is needed if S1_image==TRUE. Please install it.", call. = FALSE)
        }
    }
    if (!methods::is(S1_makeSam, "logical")) {
        stop("S1_makeSam: ", S1_makeSam, " has to be logical!", call. = FALSE)
    }
    #------------
    # check S1_BAMStream:
    #------------
    if (!methods::is(S1_BAMStream, "numeric")) {
        stop("S1_BAMStream: ", S1_BAMStream, " has to be a numeric!", call. = FALSE)
    } else if (S1_BAMStream <= 0) {
        stop("S1_BAMStream: ", S1_BAMStream, " has to be a positive numeric!", call. = FALSE)
    } else {
        S1_BAMStream = ceiling(S1_BAMStream)
    }
    #------------
    # check S1_genome:
    #------------
    if (!methods::is(S1_genome, "character")) {
        stop("S1_genome: ", S1_genome, " has to be a character!", call. = FALSE)
    } else if (!S1_genome %in% names(sysdata)) {
        LogFile = paste("S1_genome: ", S1_genome, ", is not a part of ", paste(names(sysdata), 
            collapse = "/"), ". If S2_BlackList==TRUE at stage 2, then no black-listed regions will be removed from the data.\n", 
            sep = "")
        warning(LogFile)
    }
    #------------
    # check BAM_file_1:
    #------------
    if (!methods::is(BAM_file_1, "character")) {
        stop("BAM_file_1: ", BAM_file_1, " variable has to be a file directory!", 
            call. = FALSE)
    } else if (!file.exists(BAM_file_1)) {
        stop("BAM_file_1: ", BAM_file_1, " file does not exist!", call. = FALSE)
    }
    # check format:
    BAM_file_1_Name = basename(BAM_file_1)
    BAM_file_1_Name = strsplit(BAM_file_1_Name, ".", fixed = TRUE)
    BAM_file_1_Name = unlist(BAM_file_1_Name)
    Format_1 = BAM_file_1_Name[length(BAM_file_1_Name)]
    if (all(!c("bam") %in% Format_1)) {
        stop("BAM_file_1: ", BAM_file_1, " has to be of .bam format!", call. = FALSE)
    } else {
        # check header:
        Header_1 = Rsamtools::scanBamHeader(files = BAM_file_1, what = "text")
        Header_1 = Header_1[[1]]$text
        Header_1 = Header_1[which(names(Header_1) %in% "@SQ")]
        if (length(Header_1) == 0) {
            stop("The BAM file ", BAM_file_1, " is missing the header section!", 
                call. = FALSE)
        }
        Header_1 = do.call(rbind, Header_1)
        # check if any header is missing:
        SNpos = which(grepl("SN:", Header_1[1, ]))
        LNpos = which(grepl("LN:", Header_1[1, ]))
        if (length(SNpos) == 0) {
            stop("BAM file ", BAM_file_1, " header is missing the SN entry!\n", call. = FALSE)
        }
        if (length(LNpos) == 0) {
            stop("BAM file ", BAM_file_1, " header is missing the LN entry!\n", call. = FALSE)
        }
    }
    #------------
    # check BAM_file_2:
    #------------
    if (!methods::is(BAM_file_2, "character")) {
        stop("BAM_file_2: ", BAM_file_2, " variable has to be a file directory!", 
            call. = FALSE)
    } else if (!file.exists(BAM_file_2)) {
        stop("BAM_file_2: ", BAM_file_2, " file does not exist!", call. = FALSE)
    }
    # check format:
    BAM_file_2_Name = basename(BAM_file_2)
    BAM_file_2_Name = strsplit(BAM_file_2_Name, ".", fixed = TRUE)
    BAM_file_2_Name = unlist(BAM_file_2_Name)
    Format_2 = BAM_file_2_Name[length(BAM_file_2_Name)]
    if (all(!c("bam") %in% Format_1)) {
        stop("BAM_file_2: ", BAM_file_2, " has to be of .bam format!", call. = FALSE)
    } else {
        # check header:
        Header_2 = Rsamtools::scanBamHeader(files = BAM_file_2, what = "text")
        Header_2 = Header_2[[1]]$text
        Header_2 = Header_2[which(names(Header_2) %in% "@SQ")]
        if (length(Header_2) == 0) {
            stop("The BAM file ", BAM_file_2, " is missing the header section!", 
                call. = FALSE)
        }
        Header_2 = do.call(rbind, Header_2)
        # check if any header is missing:
        SNpos = which(grepl("SN:", Header_2[1, ]))
        LNpos = which(grepl("LN:", Header_2[1, ]))
        if (length(SNpos) == 0) {
            stop("BAM file ", BAM_file_2, " header is missing the SN entry!\n", call. = FALSE)
        }
        if (length(LNpos) == 0) {
            stop("BAM file ", BAM_file_2, " header is missing the LN entry!\n", call. = FALSE)
        }
    }
    #------------
    # check Header_1 and Header_2:
    #------------
    if (!identical(Header_1, Header_2)) {
        stop("The headers of the two BAM files are not identical!", call. = FALSE)
    }
    cat("OK\n")
    #----------------
    # Keep track of files to delete:
    #----------------
    FileToDelete = c()
    #----------------
    # sort bam files by coordinates for index creation:
    #----------------
    # first bam:
    BAM_file_1.bai = paste(BAM_file_1, ".bai", sep = "")
    if (!file.exists(BAM_file_1.bai)) {
        cat("Sorting ", basename(BAM_file_1), " for index creation...")
        # create output bam:
        BAM_file_1.sorted = file.path(S1_AnalysisDir, paste(SA_prefix, "_BAM_1_sorted", 
            sep = ""))
        # sort:
        suppressWarnings(Rsamtools::sortBam(file = BAM_file_1, destination = BAM_file_1.sorted, 
            byQname = FALSE))
        cat("Done\n")
        BAM_file_1.sorted = paste(BAM_file_1.sorted, ".bam", sep = "")
        # create index:
        cat("Creating BAM index...")
        Rsamtools::indexBam(file = BAM_file_1.sorted)
        cat("Done\n")
        BAM_file_1 = BAM_file_1.sorted
        BAM_file_1.bai = paste(BAM_file_1, ".bai", sep = "")
        # save files to delete:
        FileToDelete = c(FileToDelete, BAM_file_1, BAM_file_1.bai)
    }
    # second bam:
    BAM_file_2.bai = paste(BAM_file_2, ".bai", sep = "")
    if (!file.exists(BAM_file_2.bai)) {
        cat("Sorting ", basename(BAM_file_2), " for index creation...")
        # create output bam:
        BAM_file_2.sorted = file.path(S1_AnalysisDir, paste(SA_prefix, "_BAM_2_sorted", 
            sep = ""))
        # sort:
        suppressWarnings(Rsamtools::sortBam(file = BAM_file_2, destination = BAM_file_2.sorted, 
            byQname = FALSE))
        cat("Done\n")
        BAM_file_2.sorted = paste(BAM_file_2.sorted, ".bam", sep = "")
        # create index:
        cat("Creating BAM index...")
        Rsamtools::indexBam(file = BAM_file_2.sorted)
        cat("Done\n")
        BAM_file_2 = BAM_file_2.sorted
        BAM_file_2.bai = paste(BAM_file_2, ".bai", sep = "")
        FileToDelete = c(FileToDelete, BAM_file_2, BAM_file_2.bai)
    }
    #----------------
    # Remove unmapped from the BAM files
    #----------------
    # destinations:
    BAMfilt1 = file.path(S1_AnalysisDir, paste(SA_prefix, "_usable_1_filt.bam", sep = ""))
    BAMfilt2 = file.path(S1_AnalysisDir, paste(SA_prefix, "_usable_2_filt.bam", sep = ""))
    # make flag for keeping the mapped only:
    MappedFlag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE, 
        isNotPassingQualityControls = FALSE, isDuplicate = FALSE)
    # make ScanBamParam:
    SBparam = Rsamtools::ScanBamParam(flag = MappedFlag)
    # split bam 1:
    cat("Filtering ", basename(BAM_file_1), "...")
    suppressWarnings(Rsamtools::filterBam(file = BAM_file_1, index = BAM_file_1.bai, 
        destination = BAMfilt1, param = SBparam, indexDestination = FALSE))
    cat("Done\n")
    # split bam 2:
    cat("Filtering ", basename(BAM_file_2), "...")
    suppressWarnings(Rsamtools::filterBam(file = BAM_file_2, index = BAM_file_2.bai, 
        destination = BAMfilt2, param = SBparam, indexDestination = FALSE))
    cat("Done\n")
    # add to FileToDelete:
    FileToDelete = c(FileToDelete, BAMfilt1, BAMfilt2)
    #----------------
    # Merge bam files
    #----------------
    cat("Merging ", basename(BAMfilt1), ", ", basename(BAMfilt2), " files...", sep = "")
    # output:
    MergedBAMpath = file.path(S1_AnalysisDir, paste(SA_prefix, "_usable_merged.bam", 
        sep = ""))
    suppressWarnings(Rsamtools::mergeBam(files = c(BAMfilt1, BAMfilt2), destination = MergedBAMpath, 
        overwrite = TRUE, byQname = FALSE, indexDestination = TRUE))
    cat("Done\n")
    MergedBAMpath.bai = paste(MergedBAMpath, ".bai", sep = "")
    # add to FileToDelete:
    FileToDelete = c(FileToDelete, MergedBAMpath, MergedBAMpath.bai)
    #----------------
    # sort by Qname:
    #----------------
    cat("Sorting ", basename(MergedBAMpath), " file by Qname...", sep = "")
    suppressWarnings(Rsamtools::sortBam(file = MergedBAMpath, destination = unlist(strsplit(MergedBAMpath, 
        ".bam")), byQname = TRUE))
    cat("Done\n")
    #----------------
    # fix mates
    #----------------
    MergedBAM = list(BAM = MergedBAMpath, BAMbai = MergedBAMpath.bai)
    PairedEndBAMpath = FixMates_main_fun(MergedBAM = MergedBAM, S1_AnalysisDir = S1_AnalysisDir, 
        S1_BAMStream = S1_BAMStream, CalledFromConvToPE_BAM = TRUE, SA_prefix = SA_prefix, 
        S1_image = S1_image, S1_genome = S1_genome)
    cat("Deleting unnecessary files.")
    unlink(x = FileToDelete, recursive = TRUE, force = TRUE)
    #----------------
    # If they need the sam files, convert PairedEndBAM to two SAM files:
    #----------------
    if (S1_makeSam) {
        GetSAMFiles_fun(PairedEndBAMpath = PairedEndBAMpath, S1_AnalysisDir = S1_AnalysisDir, 
            SA_prefix = SA_prefix)
    }
    cat("The paired-end BAM is in: \n", PairedEndBAMpath, "\n")
}
