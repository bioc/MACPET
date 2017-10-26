#-----------------------------------------------------------------------------#
#----------------------Input Checking Functions ------------------------------#
#-----------------------------------------------------------------------------#
#' @importFrom GenomeInfoDb genome seqlengths
#' @importFrom BSgenome available.genomes installed.genomes getBSgenome
#' @importFrom stats p.adjust.methods

#------------------------------------
#input check for PeakCallerUlt.R
#------------------------------------
InputCheckPeakCallerUlt=function(InArg){
    #global variables for Rcheck:
    pkgname=NULL
    #--------------------------
    LogFile=list()#for the log file.
    LogFile[1]="|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|\n"
    LogFile[2]="|-------------MACPET analysis input checking------------|\n"
    LogFile[3]="|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|\n"
    for(lf in 1:3) cat(LogFile[[lf]])
    #--------------------------
    #--Keep needed input:
    #--------------------------
    DefaultInputs=c("DataDir","DataFile","AnalysisDir","GenomePkg","PopImage",
                    "fileSelf","fileIntra","fileInter","BlackList",
                    "fileSelfFit","method","Stages")
    InArg=InArg[which(names(InArg)%in%DefaultInputs)]
    #--------------------------
    # analysis directory-all stages
    #--------------------------
    if(!dir.exists(InArg$AnalysisDir)){
        stop("AnalysisDir does not exist!",call. = FALSE)
    }else{
        # create log file dir:
        LogFile.dir=file.path(InArg$AnalysisDir,"MACPET_analysis.log")
        if(file.exists(LogFile.dir)) file.remove(LogFile.dir)
        # write in log file:
        for(lf in 1:3) write(LogFile[[lf]],file=LogFile.dir,append=TRUE)
    }
    #--------------------------
    # check stages-all stages
    #--------------------------
    if(!any(InArg$Stages%in%c(0:1))){
        stop("Stages argument should take one of the following elements: 0,1,c(0,1).",call. = FALSE)
    }else{
        LogFile1=paste("Stages chosen to run:",paste(InArg$Stages,collapse=" "),"\n")
        cat(LogFile1)
        write(LogFile1,file=LogFile.dir,append=TRUE)
    }
    #--------------------------
    # image inputs-stage 0, ALL
    #--------------------------
    if(c(0)%in%InArg$Stages){
        if(class(InArg$PopImage)!="logical"){
            stop("PopImage has to be logical!",call.=FALSE)
        }else{
            if(!requireNamespace("ggplot2",quietly=TRUE)){
                stop("ggplot2 needed for this function to work if PopImage==T. Please install it.",call.=FALSE)
            }
        }
    }
    #--------------------------
    # DataDir and data input:- stages 0 and ALL
    #--------------------------
    if(c(0)%in%InArg$Stages){
        InArg$Format=NA#initiate, not with NULL
        if(class(InArg$DataFile)=="character"){
            if(!dir.exists(InArg$DataDir)){
                stop("DataDir does not exist!",call.=FALSE)
            }else if(file.exists(file.path(InArg$DataDir,InArg$DataFile))){
                Format=strsplit(InArg$DataFile,".",fixed=TRUE)
                Format=unlist(Format)
                if(length(Format)<2){
                    stop("DataFile has given wrong input!",call.=FALSE)
                }else if(!Format[length(Format)]%in%c("bam","sam")){
                    stop("DataFile is neither a BAM nor a SAM format!",
                         call.=FALSE)
                }else{
                    InArg$Format=Format
                    LogFile2=paste("Format detected:",Format[length(Format)],"\n")
                    cat(LogFile2)
                    write(LogFile2,file=LogFile.dir,append=TRUE)
                }
            }else{
                stop("DataFile does not exist in the DataDir directory!",
                     call.=FALSE)
            }
        }else if(class(InArg$DataFile)!="GAlignmentPairs"){
            stop("DataFile has to be a BAM/SAM file name, or an object of class GAlignmentPairs.",call.=FALSE)
        }else{
            # then a GAlignmentPairs object is given
            LogFile2=paste("DataFile in a GAlignmentPairs object.\n")
            cat(LogFile2)
            write(LogFile2,file=LogFile.dir,append=TRUE)
        }

    }
    #--------------------------
    # Test the genome: stages 0 and ALL
    #--------------------------
    if(c(0)%in%InArg$Stages){
        if(!is.character(InArg$GenomePkg)){
            stop("GenomePkg has to be of chracter class!",call. = FALSE)
        }else if(!InArg$GenomePkg%in%BSgenome::available.genomes()){
            stop("GenomePkg is not a part of the available.genomes! See ??BSgenome::available.genomes",call.=FALSE)
        }else if(!InArg$GenomePkg%in%BSgenome::installed.genomes()){
            stop("GenomePkg is not installed! See ??BSgenome::installed.genomes",call.=FALSE)
        }else{
            GenInfo=subset(BSgenome::installed.genomes(splitNameParts=TRUE),
                           pkgname==InArg$GenomePkg)
            ChromLengths=BSgenome::getBSgenome(InArg$GenomePkg)
            ChromLengths=GenomeInfoDb::seqlengths(ChromLengths)
            ChromLengths=data.frame(Chrom=names(ChromLengths),size=ChromLengths,
                                    stringsAsFactors=FALSE)
            rownames(ChromLengths)=NULL
            #save:
            InArg$GenInfo=GenInfo
            InArg$ChromLengths=ChromLengths
        }
    }
    #--------------------------
    # load black list:-stages 0 and ALL
    #--------------------------
    if(c(0)%in%InArg$Stages){
        if(is.data.frame(InArg$BlackList)){
            if(!c("Chrom","Region.Start","Region.End")%in%colnames(InArg$BlackList)){
                stop("Give correct colnames to BlackList if it is given as data.frame!",call.=FALSE)
            }
        }else if(!is.logical(InArg$BlackList)){
            InArg$BlackList=NULL
        }else if(InArg$BlackList){
            hg=GenInfo$provider_version
            if(hg=="hg19"){
                InArg$BlackList=sysdata$Black_list_hg19
            }else{
                InArg$BlackList=NULL
            }
        }else{
            InArg$BlackList=NULL
        }
    }
    #--------------------------
    #check method: stage 1 and ALL
    #--------------------------
    if(c(1)%in%InArg$Stages){
        if(!InArg$method%in%stats::p.adjust.methods){
            stop("method value is wrong!",call. = FALSE)
        }
    }
    #--------------------------
    #file names inputs: stage 0 and ALL
    #--------------------------
    if(c(0)%in%InArg$Stages){
        if(!is.character(InArg$fileSelf)|!is.character(InArg$fileIntra)|
           !is.character(InArg$fileInter)){
            stop("fileSelf, fileIntra and fileInter must be characters!",call. = FALSE)
        }
    }
    #--------------------------
    #file names inputs: stage 1 and ALL
    #--------------------------
    if(c(1)%in%InArg$Stages){
        # check if file name is correct
        if(!is.character(InArg$fileSelfFit)){
            stop("fileSelfFit has to be character!",call. = FALSE)
        }
        # if stage 0 or ALL is NOT in the stages, then the fileSelf has to be
        # loaded and it has to be PSelf class!
        if(!c(0)%in%InArg$Stages){
            # means only stage 1 is to be run.
            # check if file name correct:
            if(!is.character(InArg$fileSelf)){
                stop("fileSelf must be character!",call. = FALSE)
            }
            # load self object:
            LoadPSelf=file.path(InArg$AnalysisDir,InArg$fileSelf)
            if(!file.exists(LoadPSelf)){
                # file does not exist:
                stop("fileSelf does not exist!",call. = FALSE)
            }
            load(LoadPSelf)#load file
            Selfobject=get(InArg$fileSelf)#get object
            if(class(Selfobject)!="PSelf"){
                stop("fileSelf is not PSelf class!",call. = FALSE)
            }else{
                InArg$Selfobject=Selfobject#save it
                LogFile4="PSelf Class object loaded.\n"
                cat(LogFile4)
                write(LogFile4,file=LogFile.dir,append=TRUE)
            }
        }
    }
    #--------------------------
    # Finallize
    #--------------------------
    LogFile5="All inputs correct! Starting MACPET analysis...\n"
    cat(LogFile5)
    write(LogFile5,file=LogFile.dir,append=TRUE)
    # return:
    InArg$LogFile.dir=LogFile.dir
    return(InArg)
}
#---------------------
#---------------------

