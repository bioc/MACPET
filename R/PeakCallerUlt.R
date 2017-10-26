#--#' make title:
#' @title Complete Binding Site Analysis Function
#--#' author:
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#--#' make reference on the article:
#' @references
#' \insertRef{macpetcite}{MACPET}
#'
#' \insertRef{ENCODE_1}{MACPET}
#--#' for the current function to be exported on the NAMESPACE:
#' @export
#--#' make discription:
#' @description \code{PeakCallerUlt} is used for running a complete binding site
#' analysis of ChIA-PET data.
#'
#--#' make details:
#' @details \code{PeakCallerUlt} runs a complete or partial binding site analysis
#' for ChIA-PET data, depending on the stages of the analysis the user wants to run.
#' The stages of the analysis are the following:
#' \describe{
#' \item{Stage 0:}{
#' PET classification stage. This stage takes the BAM/SAM ChIA-PET file, or an object of \code{\link[GenomicAlignments]{GAlignmentPairs}} class  as input
#' and it seperates the PETs into three categories:
#' Inter-chromosomal PETs (which connect two different chromosomes),
#' Intra-chromosomal PETs (which connect regions of the same chromosome) and
#' Self-ligated PETs (which are used for binding site analysis).
#' Self-ligated PETs are used for finding the protein binding sites (peaks),
#' while Intra- and Inter-chromosomal are used for interactions between
#' the peaks. Furthermore, it removes identically mapped PETs for reducing noise created by
#' amplification procedures. The algorithm uses the elbow-method to seperate the
#' Self-ligated from the Intra-chromosomal population.  Note that loading the data
#' into R might take a while depending on the size of the data.
#'
#' }
#' \item{Stage 1:}{ Peak calling stage. This stage uses the Self-ligated PETs and it runs the
#' EM algorithm to find clusters which represent candidate peaks/binding sites in
#' 2 dimentional space using skewed generalized students-t distributions (SGT).
#'  After the peak-calling analysis is done, the algorithm assesses the
#' significance of the candidate peaks using a local Poisson model.}
#' }
#'
#' @section Parallel:
#' The function can be run in parallel using the
#' \code{\link[BiocParallel]{register}}
#' function. The user has to register a parallel backhead before starting the
#' function. However the algorithm is mainly written in C++ and it is therefore fast,
#' therefore no need of parallel
#' backhead is necessary.
#'
#' @seealso
#' \linkS4class{PSelf}, \linkS4class{PIntra},
#' \linkS4class{PInter},\code{\link{summary}},
#' \code{\link{AnalysisStatistics}}, \code{\link{plot}}
#' \code{\link{BiocParallel}},\code{\link{ConvertToPSelf}},
#' \code{\link{exportPeaks}},\code{\link{TagsToGInteractions}},
#'  \code{\link{PeaksToGRanges}},  \code{\link{PeaksToNarrowPeak}}
#'
#' @return All outputs are saved at the \code{AnalysisDir}. The output depents of the stages run:
#' \describe{
#' \item{Stage 0:}{
#' \describe{
#' \item{\code{fileSelf: }}{An object of \code{\linkS4class{PSelf}} class with
#' the Self-ligated PETs. The name is specified by the \code{fileSelf} argument.
#'  This output will be further used in stage 1.}
#' \item{\code{fileIntra: }}{An object of \code{\linkS4class{PIntra}} class with
#' the Intra-chromosomal PETs. The name is specified by the \code{fileIntra}
#' argument.}
#' \item{\code{fileInter: }}{An object of \code{\linkS4class{PInter}} class with
#' the Inter-chromosomal PETs. The name is specified by the \code{fileInter}
#' argument.}
#' \item{\code{Image: }}{If \code{PopImage=TRUE}, an image for the elbow method of the
#'  Self-ligated population and its cut-off from the Intra-ligated population.}
#' }
#' }
#' \item{Stage 1:}{
#' \describe{
#' \item{\code{fileSelfFit: }}{An object of \code{\linkS4class{PSFit}} class with
#' the peak information. The name is specified by the \code{fileSelfFit}
#' argument.}
#' }
#' }
#' \item{Stage 0:1:}{ All the above outputs.}
#' }
#' The function also creates a file named \code{MACPET_analysis.log}
#' with the progress of the analysis.
#'
#'
#--#' function parameters:
#' @param DataDir A string with the directory of the BAM/SAM file (stage 0 parameter). In case \code{DataFile} is an object of
#' \code{\link[GenomicAlignments]{GAlignmentPairs}} class, \code{DataDir} will not be taken into account.
#' @param DataFile  A string with the name of the BAM/SAM file. It should end
#'  with .bam or .sam accordingly (stage 0 parameter). Alternatively, an object of \code{\link[GenomicAlignments]{GAlignmentPairs}}
#'  class with the paired-end tags to be analyzed.
#' @param AnalysisDir  A string with the directory where you want the output
#' of the algorithm to be saved (stage 0,1 parameter).
#' @param GenomePkg A string with the name of the genome as defined by the
#' \code{\link[BSgenome]{available.genomes}}. Note that the genome has to be
#' installed before running the algorithm (stage 0 parameter).
#' @param fileSelf  A string for the name of the Self-ligated PETs produced
#' as output (stage 0), and taken as input (stage 1).
#' @param fileIntra A string for the name of the Intra-chromosomal PET sproduced
#' as output (stage 0 parameter).
#' @param fileInter A string  for the name of the Inter-chromosomal PETs produced
#' as output (stage 0 parameter).
#' @param PopImage TRUE or FALSE indicating whether you want to create a figure
#'  for the Self-ligated borders/population or not (stage 0 parameter).
#'@param BlackList TRUE or FALSE depending on whether the user wants to remove
#' peaks which overlap with black listed regions before inference. Currently
#'only available for \code{hg="hg19"}, however the user can provide his own as a
#'data.frame with the following columns (stage 0 parameter):
#'\describe{
#'\item{\code{Chrom}}{A character vector with the chromosome names ("chr1" etc)}
#'\item{\code{Region.Start}}{An integer vector with the start of the
#'black listed regions.}
#'\item{\code{Region.End}}{An integer vector with the end of the black
#'listed regions.}
#'}
#'
#' @param fileSelfFit A string for the name of the fitted Self-ligated PETs produced
#' as output (stage 1 parameter).
#'
#' @param method Which method to use for adjusting for FDR for the candidate peaks found,
#' using the  \code{\link[stats]{p.adjust}} function (stage 1 parameter).
#'
#' @param Stages A numeric (0 or 1) or a vector (c(0,1)), indicating which stage(s) the
#' algorithm should run. If the vector is given, both stages will be run at once.
#'
#'
#' @examples
#' #Create a test forder on the desktop, or anywhere you want:
#' AnalysisDir=file.path(path.expand('~'),'Desktop')
#' dir.create(file.path(AnalysisDir,"MACPET.test"))
#' AnalysisDir=file.path(AnalysisDir,"MACPET.test")#where you will save
#' #the results
#'
#' #load sample data to use in the algorithm and give inputs:
#' DataDir=system.file("extdata", package = "MACPET") #data path
#' DataFile="SampleChIAPETData.bam" #data name
#' PopImage=TRUE #Sample data, not very good classifiction results.
#' GenomePkg="BSgenome.Hsapiens.UCSC.hg19" #genome of the data.
#' fileSelf="pselfData" #name for Self-ligated
#' fileIntra="pintraData" #name for Intra-chromosomal
#' fileInter="pinterData" #name for Inter-chromosomal
#' BlackList=TRUE #remove PETs in black listed regions
#' Stages=c(0)#run single stage for the example.
#'
#' #parallel backhead can be created using the BiocParallel package
#' #parallel backhead can be created using the BiocParallel package
#' #requireNamespace("BiocParallel")
#' #snow <- BiocParallel::SnowParam(workers = 1, type = "SOCK", progressbar=FALSE)
#' #BiocParallel::register(snow, default=TRUE)
#'
#' #-run for the whole binding site analysis:
#' PeakCallerUlt(DataDir=DataDir,
#'               DataFile=DataFile,
#'               AnalysisDir=AnalysisDir,
#'               GenomePkg=GenomePkg,
#'               fileSelf=fileSelf,
#'               fileIntra=fileIntra,
#'               fileInter=fileInter,
#'               PopImage=PopImage,
#'               BlackList=BlackList,
#'               Stages=Stages)
#'
#'
#'
#' #load results:
#' load(file.path(AnalysisDir,fileSelf))
#' class(pselfData) # see methods for this class
#'
#' load(file.path(AnalysisDir,fileIntra))
#' class(pintraData) # see methods for this class
#'
#' load(file.path(AnalysisDir,fileInter))
#' class(pinterData) # see methods for this class
#'
#' #-----delete test directory:
#' unlink(AnalysisDir,recursive=TRUE)

PeakCallerUlt=function(DataDir="",DataFile=NULL,AnalysisDir="",GenomePkg=NULL,
                       fileSelf=NULL,fileIntra=NULL,fileInter=NULL,PopImage=TRUE,
                       BlackList=TRUE,fileSelfFit=NULL,method="BH",Stages=c(0:1)){
    #--------------------------------------------
    #---------------Take and reorder Input:
    #--------------------------------------------
    # Take time:
    Analysis.time.start=Sys.time()
    # get arguments:
    InArg=list(DataDir=DataDir,DataFile=DataFile,AnalysisDir=AnalysisDir,
             GenomePkg=GenomePkg,fileSelf=fileSelf,fileIntra=fileIntra,
             fileInter=fileInter,PopImage=PopImage,fileSelfFit=fileSelfFit,
             method=method,BlackList=BlackList,Stages=Stages)
    #--------------------------------------------
    #---------------Check input is correct:
    #--------------------------------------------
    InArg=InputCheckPeakCallerUlt(InArg=InArg)
    #--------------------------------------------
    #---------------Run stages:
    #--------------------------------------------
    if(0%in%InArg$Stages){
        #------------------------------------------------------------------#
        #--------------------- Run PETClassification-----------------------#
        #------------------------------------------------------------------#
        LogFile=list()#for the log file.
        LogFile[1]="|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|\n"
        LogFile[2]="|------------Starting classification process------------|\n"
        LogFile[3]="|-----------------------Stage 0-------------------------|\n"
        LogFile[4]="|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|\n"
        for(lf in 1:4) cat(LogFile[[lf]])
        for(lf in 1:4) write(LogFile[[lf]],file=InArg$LogFile.dir,append=TRUE)
        #create input:
        InArg_S0=which(names(InArg)%in%c("DataDir","DataFile","AnalysisDir",
                                        "BlackList","fileSelf",
                                        "fileIntra","fileInter","PopImage",
                                        "LogFile.dir","Format",
                                        "GenInfo","ChromLengths"))
        InArg_S0=InArg[InArg_S0]
        #call Stage 0:
        do.call(what=PETClassification_fun,args=InArg_S0)

    }
    if(1%in%InArg$Stages){
        #------------------------------------------------------------------#
        #------------------------  Run PeakFinder -------------------------#
        #------------------------------------------------------------------#
        LogFile=list()#for the log file.
        LogFile[1]="|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|\n"
        LogFile[2]="|-------------Starting Binding Site Analysis------------|\n"
        LogFile[3]="|-----------------------Stage 1-------------------------|\n"
        LogFile[4]="|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|\n"
        for(lf in 1:4) cat(LogFile[[lf]])
        for(lf in 1:4) write(LogFile[[lf]],file=InArg$LogFile.dir,append=TRUE)
        #load self.data:
        if(!c(0)%in%InArg$Stages){
            # then the data is loaded:
            Selfobject=InArg$Selfobject
        }else{
            # then step 0 is run so load it
            load(file.path(InArg$AnalysisDir,InArg$fileSelf))
            Selfobject=get(InArg$fileSelf)
        }
        #create input:
        InArg_S1=list(AnalysisDir=InArg$AnalysisDir,Selfobject=Selfobject,
                    fileSelfFit=InArg$fileSelfFit,method=InArg$method,
                    LogFile.dir=InArg$LogFile.dir)
        #call stage 1:
        do.call(what=PeakFinder_fun,args=InArg_S1)
    }
    # finallize:
    Analysis.time.end=Sys.time()
    Total.Time=Analysis.time.end-Analysis.time.start
    LogFile[5]=paste("Total analysis time: ",Total.Time," ",
                   units(Total.Time),"\n",sep="")
    cat(LogFile[[5]])
    write(LogFile[[5]],file=InArg$LogFile.dir,append=TRUE)#write in log file.
    LogFile[6]="Global Analysis in done!\n"
    write(LogFile[[6]],file=InArg$LogFile.dir,append=TRUE)
    #remove chuncks:
    rm(list=ls())
    return(cat("Global Analysis in done!\n"))
}


