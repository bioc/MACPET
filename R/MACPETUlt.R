#----------#' make title:
#' @title Paired-end Tag (PET) Analysis Function.
#----------#' author:
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#----------#' make reference on the article:
#' @references
#' Vardaxis I, Drabløs F, Rye M and Lindqvist BH (2018). \emph{MACPET: Model-based Analysis for ChIA-PET}.
#' To be published.
#'
#' Consortium EP (2012) \emph{An integrated encyclopedia of DNA elements in the human genome.}.
#'  Nature, 489(7414), pp. 57–74. \url{http://dx.doi.org/10.1038/nature11247}.
#'
#----------#' for the current function to be exported on the NAMESPACE:
#' @export
#----------#' make discription:
#' @description \code{MACPETUlt} is used for running analysis based on paired-end DNA data, including stages
#' for linker removal, mapping to the reference genome, PET classification and binding site identification.
#'
#----------#' make details sections
#'@details
#'
#' Every stage has parameters associated with it. Parameters with prefix \code{SA} correspond to all
#' stages, \code{S0} to Stage 0, \code{S1} to Stage 1 etc. Parameters with \code{SA} prefix are mandatory
#' for every stage.
#'
#'If \code{SA_stages} parameter is given as vector, then the vector has to be continuous,
#'that is for example c(0:3) or c(2:3), not c(0,2,3). In general the best practice is to run
#'all the stages at once.
#'
#'The fastq files in \code{S0_fastq1} and \code{S0_fastq2} have to be of same length and be sorted
#'by their ID. Furthermore, the IDs in \code{S0_fastq1} have to end with /1 and the ones in
#'\code{S0_fastq2} with /2, representing the 5- and 3-end tags respectively. In other words,
#'for the same line in  \code{S0_fastq1} and  \code{S0_fastq2}, their IDs have to be
#'identical, except form their suffixes /1 and /2 respectively. Moreover, the "/" symbol
#'can be replaced with any other symbol, this will not cause any problems.
#'
#'\code{S0_LinkerOccurence} parameter defines the linker-occurence mode and separates the
#'usable from the ambiguous PETs. PETs with both reads including linkers are
#' not affeted by \code{S0_LinkerOccurence}. Also, reads which do not meet the
#'   \code{S0_MaxReadLength}/\code{S0_MinReadLength} lengths, are moved to
#' ambiguous anyway. The four values of \code{S0_LinkerOccurence} are:
#' \describe{
#' \item{\code{Mode 0:}}{ Both reads have to include a linker in order to
#'   be checked as usable or chimeric, if they dont, they are moved to ambiguous.}
#' \item{\code{Mode 1:}}{ If read 1 is not matching any linker, but read 2 does,
#' then the PET will be moved to usable.}
#' \item{\code{Mode 2:}}{ If read 2
#'  is not matching any linker, but read 1 does,  then the PET will be moved to usable.}
#' \item{\code{Mode 3:}}{ If any of the reads does not match any linker then the PET they will be moved to
#' usable.}
#' \item{\code{Mode 4:}}{ If both reads do not match any of the linkers, then the PET will be moved to usable. }
#' }
#'
#'
#'
#' \code{S0_MaxReadLength} has to be greater than \code{S0_MinReadLength}. The user should leave
#' those two at default unless the PET data is produced by tagmentation.
#'
#' \code{S1_fastq1_usable_dir} and \code{S1_fastq2_usable_dir} are not mandatory if Stage 0 is
#' run right before Stage 1 (\code{SA_stages}=c(0,1)). Those two are only mandatory if
#' Stage 1 is run separately. Then those parameters assume to have the usable reads only. The same
#' fastq specifications apply as those for \code{S0_fastq1} and \code{S0_fastq2}.
#'
#'The parameter \code{S1_genome} is very important. First the genome name given in \code{S1_genome}
#'should be the same as the one used for building the bowtie index for mapping.
#'This parameter will add an 'AS' column to the paired-end bam file with the genome information.
#'In Stage 2, this header will be used for identifying which kind on black-listed regions
#'to use if \code{S2_BlackList==TRUE}.
#'
#'If \code{S1_RbowtieIndexBuild==FALSE} then the bowtie index is assumed to be already built and saved in
#' \code{S1_RbowtieIndexDir}. Then the \code{S1_RbowtieIndexDir} folder should include the following
#' files: \code{S1_RbowtieIndexPrefix.1.ebwt}, \code{S1_RbowtieIndexPrefix.2.ebwt},
#'\code{S1_RbowtieIndexPrefix.3.ebwt}, \code{S1_RbowtieIndexPrefix.4.ebwt},
#'\code{S1_RbowtieIndexPrefix.rev.1.ebwt} and \code{S1_RbowtieIndexPrefix.rev.2.ebwt},
#'or with .ebwtl. Where \code{S1_RbowtieIndexPrefix} is also given as input.
#'
#'If \code{S1_RbowtieIndexBuild==TRUE} then the bowtie index will be build using the
#' \code{\link[Rbowtie:Rbowtie]{bowtie_build}} function. This function will need the .fa files
#' which should be given as input in the \code{S1_RbowtieRefDir} vector. This is a character vector
#' with the directories of the .fa files to use. The output index will be saved in
#' \code{S1_RbowtieIndexDir}. if  \code{S1_RbowtieIndexBuild==FALSE} then  \code{S1_RbowtieRefDir}
#' can be an empty string.
#'
#' The parameter \code{S2_PairedEndBAMpath} has to be specified only if Stage 2 is run without running
#' Stage 1 right before (SA_Stages=c(2) or c(2,3), not c(1,2) or c(0,1,2) for example).
#' If this is the case, the \code{S2_PairedEndBAMpath} has to be the path to the BAM/SAM paired-end file.
#' The file has to include the header with the 'SN', 'LN' and 'AS' columns. Moreover the mate flags
#' of the file have to be correct and also the duplicated PETs must be flagged too. Stage 2 will
#' upload the whole data in R using \code{\link[GenomicAlignments:readGAlignments]{readGAlignmentPairs}} function
#' with flags \code{isDuplicate=FALSE} and \code{isPaired=TRUE}. So if duplicated PETs are not
#' flagged, they will be used in the analysis. If the previous
#' stages are run in sequence, then  \code{S2_PairedEndBAMpath} will be overwritten with
#' the newly created BAM file, which will have the correct flags.
#'
#' If \code{S2_BlackList==TRUE} then which genome black-list is going to be used
#' is decided by the 'AS' column in the \code{S2_PairedEndBAMpath} file, which is
#' specified by the \code{S1_genome} if Stage 1 is also run. The black-listed regions cover
#' the following genomes: 'hg19', 'ce10', 'dm3', 'hg38', 'mm9', 'mm10'.
#' If the 'AS' header column is missing from the \code{S2_PairedEndBAMpath} file, or
#' if the  \code{S1_genome} is not matching any of the above named genomes, then
#' a warning will be produced saying that no black-listed regions will be removed.
#' Alternatily, the user can provide its own black-listed regions as a  \code{\linkS4class{GRanges}}
#' object.
#'
#' The parameter \code{S3_fileSelfDir} is not mandatory if the stages are run in sequence,
#' if Stage 2 is run right before stage 3. If this is the case then \code{S3_fileSelfDir}
#' will be overwritten with the data produced in Stage 2. If Stage 3 is run
#' separately, then \code{S3_fileSelfDir} has to be provided. It should be
#' a \code{\linkS4class{PSelf}} object and both the name of the object in the directory
#' and the one uploaded in R should be  \code{SA_prefix_pselfData}.
#'
#----------#' make stages description:
#' @section Stages description: \code{MACPETUlt} runs a complete or partial analysis
#' for PET data, depending on the stages of the analysis the user wants to run.
#' The stages of the analysis are the following:
#' \describe{
#' \item{Stage 0:}{
#' Linker identification stage: This stage uses the two fastq files for the 5- and 3-end tags and identifies
#'  which tags contain any of the linkers. Based on the
#'  linker combinations it classifies the PETs as usable (linkers A/A or B/B),
#'   chimeric (linkers A/B or B/A)
#'   and ambiguous (linkers non/A, non/B, A/non, B/non unless chosen otherwise by \code{S0_LinkerOccurence},
#'   or be smaller/bigger than the \code{S0_MinReadLength}/\code{S0_MaxReadLength} after the linker removal, respectively).
#'    Only usable PETs are considered in the subsequent steps.
#' }
#' \item{Stage 1:}{
#' PET mapping stage: This stage uses the usable PETs identified by stage 0. It maps them separately to the
#' reference genome using  the \code{\link[Rbowtie:Rbowtie]{bowtie}} function with no mismatch per read,
#' and keeps the uniquely mapped reads only. It then maps the unmapped reads to the reference genome
#' with at most one mismatch and keeps the uniquely mapped reads. Uniquely mapped reads with
#' zero or one mismatch are then merged and paired, their duplicates are marked and a paired-end bam file
#'  is created which is used in State 2.
#' }
#' \item{Stage 2:}{
#' PET classification stage: This stage takes the BAM paired-end file from stage 1 and
#' classifies the PETs as:
#' Inter-chromosomal PETs (which connect two different chromosomes),
#' Intra-chromosomal PETs (which connect regions of the same chromosome) and
#' Self-ligated PETs (which are used for binding site analysis).
#' Self-ligated PETs are used for finding the protein binding sites (peaks),
#' while Intra- and Inter-chromosomal are used for interactions between
#' the peaks. The algorithm uses the elbow-method to seperate the
#' Self-ligated from the Intra-chromosomal population.  Note that loading the data
#' into R might take a while depending on the size of the data.
#'
#' }
#' \item{Stage 3:}{ Peak calling stage: This stage uses the Self-ligated PETs and it runs the
#' EM algorithm to find clusters which represent candidate peaks/binding sites in
#' 2 dimentional space using skewed generalized students-t distributions (SGT).
#'  After the peak-calling analysis is done, the algorithm assesses the
#' significance of the candidate peaks using a local Poisson model.}
#' }
#'
#' @section Parallel:
#' All stages can be run in parallel using the
#' \code{\link[BiocParallel:register]{register}}
#' function. The user has to register a parallel backhead before starting the
#' function.
#----------#'
#' @seealso
#' \linkS4class{PSelf}, \linkS4class{PIntra},
#' \linkS4class{PInter},\code{\link{summary}},
#' \code{\link{AnalysisStatistics}}, \code{\link{plot}}
#' \code{\link{BiocParallel}},\code{\link{ConvertToPSelf}},
#' \code{\link{exportPeaks}},\code{\link{TagsToGInteractions}},
#'  \code{\link{PeaksToGRanges}},  \code{\link{PeaksToNarrowPeak}},
#'  \code{\link{ConvertToPE_BAM}}
#----------#'output
#' @return All outputs are saved at the \code{SA_AnalysisDir}. The output depents of the stages run:
#' \describe{
#' \item{Stage 0: (outputs saved in a folder named \code{S0_results} in  \code{SA_AnalysisDir})}{
#' \describe{
#' \item{\code{SA_prefix_usable_1.fastq.gz: }}{fastq.gz files with the usable 5-end tags. To be used in Stage 1.}
#' \item{\code{SA_prefix_usable_2.fastq.gz: }}{fastq.gz files with the usable 3-end tags. To be used in Stage 1.}
#' \item{\code{SA_prefix_chimeric_1.fastq.gz: }}{fastq.gz files with the chimeric 5-end tags.}
#' \item{\code{SA_prefix_chimeric_2.fastq.gz: }}{fastq.gz files with the chimeric 3-end tags.}
#' \item{\code{SA_prefix_ambiguous_1.fastq.gz: }}{fastq.gz files with the ambiguous 5-end tags.}
#' \item{\code{SA_prefix_ambiguous_2.fastq.gz: }}{fastq.gz files with the ambiguous 3-end tags.}
#' \item{\code{SA_prefix_stage_0_image.jpg: }}{Pie chart image with the split of two fastq files used as input (if \code{S0_image==TRUE}).}
#' }
#' }
#' \item{Stage 1: (outputs saved in a folder named \code{S1_results} in  \code{SA_AnalysisDir})}{
#' \describe{
#' \item{\code{SA_prefix_usable_1.sam: }}{sam file with the mapped 5-end reads (if \code{S1_makeSam==FALSE}).}
#' \item{\code{SA_prefix_usable_2.sam: }}{sam file with the mapped 3-end reads (if \code{S1_makeSam==FALSE}).}
#' \item{\code{SA_prefix_Paired_end.bam: }}{paired-end bam file with the mapped PETs. To be used in Stage 2}
#' \item{\code{SA_prefix_Paired_end.bam.bai: }}{.bai file for \code{SA_prefix_Paired_end.bam}. To be used in Stage 2.}
#' \item{\code{SA_prefix_stage_1_p1_image.jpg: }}{Pie-chart for the mapped/unmapped reads from \code{SA_prefix_usable_1.sam}, \code{SA_prefix_usable_2.sam} (if \code{S1_image==TRUE}).}
#' \item{\code{SA_prefix_stage_1_p2_image.jpg: }}{Pie-chart for the paired/unpaired reads of \code{SA_prefix_Paired_end.bam} (if \code{S1_image==TRUE}).}
#' }
#' }
#' \item{Stage 2: (outputs saved in a folder named \code{S2_results} in  \code{SA_AnalysisDir})}{
#' \describe{
#' \item{\code{SA_prefix_pselfData: }}{An object of \code{\linkS4class{PSelf}} class with
#' the Self-ligated PETs. To be used in Stage 3.}
#' \item{\code{SA_prefix_pintraData: }}{An object of \code{\linkS4class{PIntra}} class with
#' the Intra-chromosomal PETs.}
#' \item{\code{SA_prefix_pinterData: }}{An object of \code{\linkS4class{PInter}} class with
#' the Inter-chromosomal PETs.}
#' \item{\code{SA_prefix_stage_2_p1_image.jpg: }}{Pie-chart reliable/dublicated/black-listed PETs of \code{SA_prefix_Paired_end.bam} (if \code{S2_image==TRUE}).}
#' \item{\code{SA_prefix_stage_2_p2_image.jpg: }}{Histogram with the self-ligated/intra-chromosomal cut-off for \code{SA_prefix_Paired_end.bam} (if \code{S2_image==TRUE}).}
#' \item{\code{SA_prefix_stage_2_p3_image.jpg: }}{Pie-chart for the self-ligated/intra-chromosomal/inter-chromosomal PETs of \code{SA_prefix_Paired_end.bam} (if \code{S2_image==TRUE}).}
#' }
#' }
#' \item{Stage 3: (outputs saved in a folder named \code{S3_results} in  \code{SA_AnalysisDir})}{
#' \describe{
#' \item{\code{SA_prefix_psfitData: }}{An object of \code{\linkS4class{PSFit}} class with the peak information.}
#' \item{\code{SA_prefix_stage_3_p1_image.jpg: }}{Sizes of the upstream vs downstream peaks of each binding site given the binding site's FDR (if \code{S3_image==TRUE}).}
#' \item{\code{SA_prefix_stage_3_p2_image.jpg: }}{FDR of the binding sites. The horizontal red line is at FDR=0.05 (if \code{S3_image==TRUE}).}
#' \item{\code{SA_prefix_stage_3_p3_image.jpg: }}{Comparison of binding site sizes given their FDR (if \code{S3_image==TRUE}).}
#' \item{\code{SA_prefix_stage_3_p3_image.jpg: }}{FDR for the upstream/donwstream peaks of the binding sites given the binding sites FDR (if \code{S3_image==TRUE}).}
#' }
#' }
#' \item{Stage 0:3 :}{ All the above outputs. Furthermore, a log file named \code{SA_prefix_analysis.log} is always
#' created in \code{SA_AnalysisDir} with information about the process.}
#' }
#'
#'
#----------#' function parameters:
#'@param SA_AnalysisDir A directory were all the ouput is to be saved.
#'This parameter is mandatory for every stage.
#'@param SA_stages Numeric vector or integer (if stages are run separately).
#'This parameter is mandatory for every stage (see details).
#'@param SA_prefix A string which is going to be used as prefix for the outputs (default: 'MACPET').
#'This parameter is mandatory for every stage.
#'
#'@param S0_fastq1 A string with the directory of the 5-end fastq (or fastq.gz) file.
#'This parameter is mandatory if Stage 0 is run (see details).
#'@param S0_fastq2 A string with the directory of the 3-end fastq (or fastq.gz) file.
#'This parameter is mandatory if Stage 0 is run (see details).
#'@param S0_LinkerA A string with the first linker sequence (default 'GTTGGATAAG').
#'This parameter is mandatory if Stage 0 is run (see details).
#'@param S0_LinkerB A string with the second linker sequence (default 'GTTGGAATGT').
#'This parameter is mandatory if Stage 0 is run (see details).
#'@param S0_MinReadLength A positive integer with the minimum read
#'length after linker trimming (default: 18).
#'This parameter is mandatory if Stage 0 is run (see details).
#'@param S0_MaxReadLength  A positive integer with the maximum read
#' length after linker trimming (default: 50).
#'This parameter is mandatory if Stage 0 is run (see details).
#'@param S0_LinkerOccurence One of the following: 0, 1, 2, 3, 4. This parameter
#'defines the linker-occurence mode (see details). Default 0.
#'@param S0_image Logical, indicating if a pie-chart image for the fastq files
#' classification will be produced (default=TRUE).
#'This parameter is mandatory if Stage 0 is run.
#'@param S0_fastqStream Positive integer for total lines of fastq files to be
#' loaded in R (best to leave it at default because it might cause memory crash).
#'This parameter is mandatory if Stage 0 is run.
#'
#'@param S1_fastq1_usable_dir String with the directory of the 5-end usable fastq (or fastq.gz) files.
#'This parameter might not be mandatory (see details).
#'@param S1_fastq2_usable_dir String with the directory of the 3-end usable fastq (or fastq.gz) files.
#'This parameter might not be mandatory (see details).
#'@param S1_image Logical indicating if images for the mapping percentage and
#' the pairing percentage will be produced (default=TRUE).
#'This parameter is mandatory if Stage 1 is run.
#'@param S1_BAMStream Positive integer for the total number of bam file lines
#'to be loaded in R in a loop for pairing  (best to leave it at default because it might cause memory crash).
#'This parameter is mandatory if Stage 1 is run.
#'@param S1_makeSam Logical indicating whether the resulted paired-end BAM file will be splitted to two
#'SAM files (one for each read). The output SAM files can be used as input in the MANGO algorithm (default=TRUE).
#'Note, that the user has to remove the SAM header before running MANGO.
#'This parameter is mandatory if Stage 1 is run.
#'@param S1_genome String with the genome to be used in the bam file header (default='hg19').
#'This parameter is mandatory if Stage 1 is run (see details).
#'@param S1_RbowtieIndexBuild Logical indicating whether you want to build
#'the bowtie index or not (default=FALSE).
#'This parameter is mandatory if Stage 1 is run (see details).
#'@param S1_RbowtieIndexDir String with the directory of the bowtie
#'index (if \code{S1_RbowtieIndexBuild==FALSE})
#'or with the directory where the bowtie index will be
#'saved (if \code{S1_RbowtieIndexBuild==TRUE}).
#'This parameter is mandatory if Stage 1 is run (see details).
#'@param S1_RbowtieIndexPrefix String with the prefix for the bowtie
#' indeces in \code{S1_RbowtieIndexDir} (see details).
#'This parameter is mandatory if Stage 1 is run (see details).
#'@param S1_RbowtieRefDir A vector with the directories of the .fa files,
#' used if \code{S1_RbowtieIndexBuild==TRUE}.
#'This parameter is mandatory if Stage 1 is run and \code{S1_RbowtieIndexBuild==TRUE} (see details).
#'
#'@param S2_PairedEndBAMpath A string with the directory of the paired-end bam file (or paired-end sam file).
#' This parameter might not be mandatory (see details).
#'@param S2_image Logical indicating whether images for the s
#'elf-ligated/intra-chromosomal cut-off as well as pie-charts
#'for the PET classification will be produced (default=TRUE).
#' This parameter is mandatory if Stage 1 is run.
#'@param S2_BlackList Logical indicating whether black-listed
#'regions will be removed from the data
#' based on the \code{S1_genome} parameter (see details).
#'Alternatively a \code{\linkS4class{GRanges}} object with the user specified regions.
#'This parameter is mandatory if Stage 2 is run.
#'
#'@param S3_fileSelfDir A string with the directory of the of the object of
#'class \code{\linkS4class{PSelf}}.
#' This parameter might not be mandatory (see details).
#'@param S3_image Logical indicating whether images for the binding site's FDR,
#' sizes of the binding sites, sizes of binding site's upstream/downstream peaks
#' will be created.
#' This parameter is mandatory if Stage 3 is run.
#'@param S3_method String with the FDR method used for finding
#'p-values of significant peaks in the data.
#'See  \code{\link[stats:p.adjust]{p.adjust.methods}} (default= 'BH').
#'This parameter is mandatory if Stage 3 is run.
# ----------#'
#' @examples
#'
#' #Create a temporary forder, or anywhere you want:
#' SA_AnalysisDir=file.path(tempdir(),'MACPETtest')
#' dir.create(SA_AnalysisDir)#where you will save the results
#' #give directory of the BAM file:
#' S2_PairedEndBAMpath=system.file('extdata', 'SampleChIAPETData.bam', package = 'MACPET')
#'
#' #give prefix name:
#' SA_prefix='MACPET'
#'
#' #parallel backhead can be created using the BiocParallel package
#' #parallel backhead can be created using the BiocParallel package
#' #requireNamespace('BiocParallel')
#' #snow <- BiocParallel::SnowParam(workers = 4, type = 'SOCK', progressbar=FALSE)
#' #BiocParallel::register(snow, default=TRUE)
#'
#' #-run for the whole binding site analysis:
#' MACPETUlt(SA_AnalysisDir=SA_AnalysisDir,
#'           SA_stages=c(2:3),
#'           SA_prefix=SA_prefix,
#'           S2_PairedEndBAMpath=S2_PairedEndBAMpath,
#'           S2_image=TRUE,
#'           S2_BlackList=TRUE,
#'           S3_image=TRUE)
#'
#'
#'
#' #load results:
#' SelfObject=paste(SA_prefix,'_pselfData',sep='')
#' load(file.path(SA_AnalysisDir,'S2_results',SelfObject))
#' SelfObject=get(SelfObject)
#' class(SelfObject) # see methods for this class
#'
#' IntraObject=paste(SA_prefix,'_pintraData',sep='')
#' load(file.path(SA_AnalysisDir,'S2_results',IntraObject))
#' IntraObject=get(IntraObject)
#' class(IntraObject) # see methods for this class
#'
#' InterObject=paste(SA_prefix,'_pinterData',sep='')
#' load(file.path(SA_AnalysisDir,'S2_results',InterObject))
#' InterObject=get(InterObject)
#' class(InterObject) # see methods for this class
#'
#' SelfFitObject=paste(SA_prefix,'_psfitData',sep='')
#' load(file.path(SA_AnalysisDir,'S3_results',SelfFitObject))
#' SelfFitObject=get(SelfFitObject)
#' class(SelfFitObject) # see methods for this class
#'
#' #-----delete test directory:
#' unlink(SA_AnalysisDir,recursive=TRUE)
MACPETUlt = function(SA_AnalysisDir = "", SA_stages = c(0:3), SA_prefix = "MACPET",
    S0_fastq1 = "", S0_fastq2 = "", S0_LinkerA = "GTTGGATAAG", S0_LinkerB = "GTTGGAATGT",
    S0_MinReadLength = 18, S0_MaxReadLength = 50, S0_LinkerOccurence = 0, S0_image = TRUE,
    S0_fastqStream = 2e+06, S1_fastq1_usable_dir = "", S1_fastq2_usable_dir = "",
    S1_image = TRUE, S1_BAMStream = 2e+06, S1_makeSam = TRUE, S1_genome = "hg19",
    S1_RbowtieIndexBuild = FALSE, S1_RbowtieIndexDir = "", S1_RbowtieIndexPrefix = "",
    S1_RbowtieRefDir = "", S2_PairedEndBAMpath = "", S2_image = TRUE, S2_BlackList = TRUE,
    S3_fileSelfDir = "", S3_image = TRUE, S3_method = "BH") {
    #--------------------------------------------
    #---------------Take and reorder Input:
    #--------------------------------------------
    # Take time:
    Analysis.time.start = Sys.time()
    # get arguments:
    InArg = list(SA_AnalysisDir = SA_AnalysisDir, SA_stages = SA_stages, SA_prefix = SA_prefix,
        S0_fastq1 = S0_fastq1, S0_fastq2 = S0_fastq2, S0_LinkerA = S0_LinkerA, S0_LinkerB = S0_LinkerB,
        S0_MinReadLength = S0_MinReadLength, S0_MaxReadLength = S0_MaxReadLength,
        S0_LinkerOccurence = S0_LinkerOccurence, S0_image = S0_image, S0_fastqStream = S0_fastqStream,
        S1_fastq1_usable_dir = S1_fastq1_usable_dir, S1_fastq2_usable_dir = S1_fastq2_usable_dir,
        S1_image = S1_image, S1_BAMStream = S1_BAMStream, S1_makeSam = S1_makeSam, S1_genome = S1_genome,
        S1_RbowtieIndexBuild = S1_RbowtieIndexBuild, S1_RbowtieIndexDir = S1_RbowtieIndexDir,
        S1_RbowtieIndexPrefix = S1_RbowtieIndexPrefix, S1_RbowtieRefDir = S1_RbowtieRefDir,
        S2_PairedEndBAMpath = S2_PairedEndBAMpath, S2_image = S2_image, S2_BlackList = S2_BlackList,
        S3_fileSelfDir = S3_fileSelfDir, S3_image = S3_image, S3_method = S3_method)
    #--------------------------------------------
    #---------------Check input is correct:
    #--------------------------------------------
    InArg = InputCheckMACPETUlt(InArg = InArg)
    #--------------------------------------------
    #---------------Run stages:
    #--------------------------------------------
    if (c(0) %in% SA_stages) {
        #------------------------------------------------------------------#
        #---------------------   fastq separation    ----------------------#
        #------------------------------------------------------------------#
        # take arguments:
        InArgS0 = InArg$InArgS0
        # for the log file.
        LogFile = list()
        LogFile[1] = "|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|"
        LogFile[2] = "|-------------------Filtering Linkers-------------------|"
        LogFile[3] = "|-----------------------Stage 0-------------------------|"
        LogFile[4] = "|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|"
        for (lf in seq_len(4))  futile.logger::flog.info(LogFile[[lf]], name="SA_LogFile",capture=FALSE)
        # call Stage 0:
        do.call(what = Stage_0_Main_fun, args = InArgS0)
        SA_LogFile.dir = InArgS0$SA_LogFile.dir  #used in the at the end
    }
    if (c(1) %in% SA_stages) {
        #------------------------------------------------------------------#
        #------------------- Mapping And building Paired-end BAM-----------#
        #------------------------------------------------------------------#
        # take arguments:
        InArgS1 = InArg$InArgS1
        # for the log file.
        LogFile = list()
        LogFile[1] = "|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|"
        LogFile[2] = "|--Mapping to Ref. Genome And building paired-end BAM---|"
        LogFile[3] = "|-----------------------Stage 1-------------------------|"
        LogFile[4] = "|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|"
        for (lf in seq_len(4)) futile.logger::flog.info(LogFile[[lf]], name="SA_LogFile",capture=FALSE)
        # call Stage 1:
        do.call(what = Stage_1_Main_fun, args = InArgS1)
        SA_LogFile.dir = InArgS1$SA_LogFile.dir  #used in the at the end
    }
    if (c(2) %in% SA_stages) {
        #------------------------------------------------------------------#
        #--------------------- Run PETClassification-----------------------#
        #------------------------------------------------------------------#
        # take arguments:
        InArgS2 = InArg$InArgS2
        # for the log file.
        LogFile = list()
        LogFile[1] = "|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|"
        LogFile[2] = "|--------------PET Classification Analysis--------------|"
        LogFile[3] = "|-----------------------Stage 2-------------------------|"
        LogFile[4] = "|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|"
        for (lf in seq_len(4)) futile.logger::flog.info(LogFile[[lf]], name="SA_LogFile",capture=FALSE)
        # load data, it is saved from stage 1 else it is already loaded:
        if (c(1) %in% SA_stages) {
            S2_PairedData = LoadBAM_FromMACPETUlt_fun(S2_PairedEndBAMpath = InArgS2$S2_PairedEndBAMpath,
                S2_BlackList = InArgS2$S2_BlackList,
                S2_image = InArgS2$S2_image, S2_AnalysisDir = InArgS2$S2_AnalysisDir,
                SA_prefix = InArgS2$SA_prefix)
            InArgS2$S2_PairedData = S2_PairedData
            InArgS2$S2_BlackList = NULL
            InArgS2$S2_PairedEndBAMpath = NULL
        }
        # call Stage 2:
        do.call(what = Stage_2_Main_fun, args = InArgS2)
        SA_LogFile.dir = InArgS2$SA_LogFile.dir  #used in the at the end
    }
    if (c(3) %in% SA_stages) {
        #------------------------------------------------------------------#
        #------------------------  Run PeakFinder -------------------------#
        #------------------------------------------------------------------#
        # take arguments:
        InArgS3 = InArg$InArgS3
        LogFile = list()  #for the log file.
        LogFile[1] = "|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|"
        LogFile[2] = "|----------------Binding Site Analysis------------------|"
        LogFile[3] = "|-----------------------Stage 3-------------------------|"
        LogFile[4] = "|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|"
        for (lf in seq_len(4)) futile.logger::flog.info(LogFile[[lf]], name="SA_LogFile",capture=FALSE)
        # if stage 2 is run then the pself is saved, so load it
        if (c(2) %in% SA_stages) {
            # then load the data:
            load(InArgS3$S3_fileSelfDir)
            InArgS3$S3_Selfobject = get(paste(InArgS3$SA_prefix, "_pselfData", sep = ""))
            InArgS3$S3_fileSelfDir = NULL
        }
        # call stage 3:
        do.call(what = Stage_3_Main_fun, args = InArgS3)
        SA_LogFile.dir = InArgS3$SA_LogFile.dir  #used in the at the end
    }
    #------------------------
    # finallize:
    #------------------------
    futile.logger::flog.info("|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|", name="SA_LogFile",capture=FALSE)
    futile.logger::flog.info("Global analysis is done!", name="SA_LogFile",capture=FALSE)
    # take time:
    Analysis.time.end = Sys.time()
    Total.Time = Analysis.time.end - Analysis.time.start
    LogFile = paste("Global analysis time:", Total.Time, " ", units(Total.Time))
    futile.logger::flog.info(LogFile, name="SA_LogFile",capture=FALSE)
}
# done
