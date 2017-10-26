# Helping functions for PeakFinder.R
#' @importFrom plyr dlply . ldply llply ddply
#' @importFrom intervals Intervals clusters interval_included size
#' @importFrom stats ppois p.adjust
#' @importFrom S4Vectors metadata
#' @importFrom GenomeInfoDb seqlengths seqinfo
#' @importFrom BiocParallel bplapply
#'
###############################################################################
# Main function for stage 1:
###############################################################################
#-------------
#-------------
PeakFinder_fun=function(AnalysisDir,Selfobject,fileSelfFit,method,LogFile.dir){
    #global variables for Rcheck:
    #--------------------------------------------
    #---------------Take and reorder Input:
    #--------------------------------------------
    # Take time:
    Analysis.time.start=Sys.time()
    #--------------------------------------------
    #---------------Keep data you need only:
    #--------------------------------------------
    FitDataInfo=TakeFitAnalysisData_fun(x=Selfobject)
    PETsData=FitDataInfo$PETsData#all the PETs data, used in inference too.
    ChromInf=FitDataInfo$ChromInf#used in inference and Fitting.
    #--------------------------------------------
    #---------------Segment to regions:
    #--------------------------------------------
    SegmRes=FindRegions_fun(x=PETsData,ChromInf=ChromInf,
                            LogFile.dir=LogFile.dir)
    SegPETS=SegmRes$x#segmented PETs data, used in fitting.
    RegionCounts=SegmRes$RegionCounts#take the counts
    #------------
    #-----save the region counts on main data:
    #------------
    Ordermatch=match(RegionCounts$Chrom,
                     S4Vectors::metadata(Selfobject)$Self_info$Chrom)
    S4Vectors::metadata(Selfobject)$Self_info$Region.counts[Ordermatch]=
        RegionCounts$Counts
    #--------------------------------------------
    #---------------Fit:
    #--------------------------------------------
    GlobalPeakCallRes=FitCallGlobal_fun(SegPETS=SegPETS,ChromInf=ChromInf,
                                        LogFile.dir=LogFile.dir)
    #--------------------------------------------
    #---------------update the results
    #--------------------------------------------
    #------classification information
    S4Vectors::metadata(Selfobject)$Classification.Info=
        GlobalPeakCallRes$Classification.Info
    #------Peak counts
    S4Vectors::metadata(Selfobject)$Self_info$Peak.counts=0
    Match=match(GlobalPeakCallRes$Peak.counts$Chrom,
                S4Vectors::metadata(Selfobject)$Self_info$Chrom)
    S4Vectors::metadata(Selfobject)$Self_info$Peak.counts[Match]=
        GlobalPeakCallRes$Peak.counts$V1
    #--------------------------------------------
    #---------------Run inference:
    #--------------------------------------------
    Peaks.Info=GlobalPeakCallRes$Peaks.Info#peaks found by fit
    #run inference
    InfRes=SignificanceCall_fun(Peaks.Info=Peaks.Info,PETsData=PETsData,
                                ChromInf=ChromInf,method=method,
                                LogFile.dir=LogFile.dir)
    #---add Peak information for peaks found in data with their FDR etc
    S4Vectors::metadata(Selfobject)$Peaks.Info=InfRes
    #-----update class and save:
    class(Selfobject)="PSFit"
    assign(fileSelfFit,Selfobject)#assign value.
    SavePath=file.path(AnalysisDir,fileSelfFit)
    save(list=fileSelfFit,file=SavePath)
    # save log:
    Analysis.time.end=Sys.time()
    Total.Time=Analysis.time.end-Analysis.time.start
    LogFile1=paste("Total peak-calling time: ",Total.Time," ",
                   units(Total.Time),"\n",sep="")
    cat(LogFile1)
    write(LogFile1,file=LogFile.dir,append=TRUE)#write in log file.
    #print its done:
    return(cat("Binding Site Analysis is done!\n"))
}
#-------------
#-------------
###############################################################################
# Functions for inputs and data
###############################################################################
#-------------
#-------------
#function for converting the data in the analysis format:
TakeFitAnalysisData_fun=function(x){
    cat("Converting data for analysis...\n")
    #------------
    #take Chromosome information(some used in inference):
    #------------
    ChromInf=GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(x))
    ChromInf=data.frame(Chrom=names(ChromInf),
                        size=as.numeric(ChromInf),
                        stringsAsFactors=FALSE)
    Self_info=S4Vectors::metadata(x)$Self_info
    MatchChrom=match(Self_info$Chrom,ChromInf$Chrom)
    ChromInf$PET.counts[MatchChrom]=Self_info$PET.counts
    #------------
    #convert to data frame and Swap Anchors:
    #------------
    x=as.data.frame(x)
    x$UTag=(x$start1+x$end1)/2#use Tag mids
    x$DTag=(x$start2+x$end2)/2#use Tag mids
    Swap=which(x$UTag>x$DTag)
    if(length(Swap)!=0){
        UTagswap=x$UTag[Swap]
        DTagswap=x$DTag[Swap]
        x$UTag[Swap]=DTagswap
        x$DTag[Swap]=UTagswap
    }
    #------------
    #add mainindex
    #------------
    x$MainIndex=1:nrow(x)
    #------------
    #keep left part of downstream and right of upstream
    #------------
    x$Chrom=as.character(x$seqnames1)
    x=x[,c("Chrom","UTag","DTag","MainIndex")]
    #------------
    #return:
    #------------
    Rtrn=list(PETsData=x,ChromInf=ChromInf)
    return(Rtrn)
}
#-------------
#-------------
#####################################
# Functions for region segmentation:
#####################################
#-------------
#-------------
#main function for region segmentation:
FindRegions_fun=function(x,ChromInf,LogFile.dir){
    #rcheck:
    Chrom=Region=NULL
    #------------
    #------------
    #split by chromosome to list to pass to bplapply
    bppass=plyr::dlply(x,plyr::.(Chrom),function(x) x)
    #------------
    #------------
    #--Segment to regions:
    #------------
    cat("Segmenting into regions...\n")
    SegmRes=BiocParallel::bplapply(
        X=bppass,
        FUN=function(y,ChromInf){
            ChromCur=as.character(unique(y$Chrom))
            ChromInfCur=subset(ChromInf,Chrom==ChromCur)
            #------------
            #Give the region ids:
            #------------
            y=RegionIds_fun(PETinf=y,ChromCur=ChromCur,ChromInfCur=ChromInfCur)
            return(y)
        },
        ChromInf=ChromInf)
    #------------
    #---break and rbind:
    #------------
    #take the region info:
    PETinf=lapply(SegmRes,"[[",1)
    PETinf=plyr::ldply(PETinf,function(y)y)
    PETinf=subset(PETinf,!is.na(Region))#keep not NA
    if(nrow(PETinf)!=0){
        #update data with the regions:
        x$Region=NA
        x$Region[PETinf$MainIndex]=PETinf$Region
        x=subset(x,!is.na(Region))
    }else{
        stop("No regions found in data. The data is probably too sparse.")
    }
    #------------
    #--take total regions info:
    #------------
    RegionCounts=lapply(SegmRes,"[[",2)
    RegionCounts=plyr::ldply(RegionCounts,function(y)y)
    RegionCounts=RegionCounts[,c("Chrom","Counts")]
    LogFile1=paste("Total Regions found: ",sum(RegionCounts$Counts),"\n")
    cat(LogFile1)
    write(LogFile1,file=LogFile.dir,append=TRUE)#write to log file
    return(list(x=x,RegionCounts=RegionCounts))
}
#-------------
#-------------
#Fuction for creating region ids by segmenting into regions:
RegionIds_fun=function(PETinf,ChromCur,ChromInfCur){
    #global variables for Rcheck:
    #--------------------------
    #------------
    #initiate
    #------------
    PETinf$Region=NA
    TotRegions=data.frame(Chrom=ChromCur,Counts=0)
    #------------
    #take Pet Pmid to create global intervals:
    #------------
    PETInt=intervals::Intervals(PETinf[,c("UTag","DTag")],closed=c(TRUE,TRUE))
    #------------
    #Take clusters:
    #------------
    PETclst=intervals::clusters(PETInt,which=TRUE,w=0)
    TotRegions$Counts=length(PETclst)#total regions
    if(TotRegions$Counts!=0){
        #------------
        #region ids:
        #------------
        Region=rep(1:length(PETclst),lengths(PETclst))
        #------------
        #pets in the regions:
        #------------
        TagsIncluded=unlist(PETclst)
        #------------
        #add them according to the rownames of the PETinf
        #------------
        PETinf$Region[TagsIncluded]=Region
    }#else no region found
    PETinf=PETinf[,c("Chrom","MainIndex","Region")]
    return(list(PETinf=PETinf,TotRegions=TotRegions))
}
#######################################
# Functions for fitting the main model:
########################################
#-------------
#-------------
#main function for fitting: breaking each region
FitCallGlobal_fun=function(SegPETS,ChromInf,LogFile.dir){
    #---R-check
    Region=Chrom=NULL
    #----
    bppass=plyr::dlply(SegPETS,plyr::.(Chrom,Region),function(y) y)
    #---------------------------------------
    #--------call the function in parallel for fitting regions:
    #---------------------------------------
    cat("Running peak calling process...\n")
    #the following runs in parallel:
    GlobalFitRes=BiocParallel::bplapply(X=bppass,FUN=FitCallLocal_fun_Rcpp,
                                        ChromInf=ChromInf)
    cat("Fit completed!\n")
    #---------------------------
    #--------save Classification:
    #---------------------------
    Peak.Id.Inf=plyr::llply(GlobalFitRes,function(x) x[[1]])
    Peak.Id.Inf=do.call(rbind,Peak.Id.Inf)
    if(is.null(Peak.Id.Inf)){
        stop("No Peaks found by the algorithm!")
    }
    #---------------------------
    #--mach classes with SegPETS
    #---------------------------
    SegPETS$Peak.ID=NA
    Match=match(Peak.Id.Inf[,2],SegPETS$MainIndex)
    SegPETS$Peak.ID[Match]=Peak.Id.Inf[,1]
    Classification.Info=SegPETS[,c("MainIndex","Region","Peak.ID")]#to return
    #---------------------------
    #--------save Peaks:
    #---------------------------
    Peaks.Info=plyr::ldply(GlobalFitRes,function(x) x[[2]])
    if(is.null(Peaks.Info)){
        stop("No Peaks found by the algorithm!")
    }
    Peaks.Info=Peaks.Info[,which(colnames(Peaks.Info)!=".id")]
    LogFile1=paste("Total ",nrow(Peaks.Info)," candidate peaks found in data.\n")
    cat(LogFile1)
    write(LogFile1,file=LogFile.dir,append=TRUE)#write to log file
    #---------------------------
    #--------Find counts:
    #---------------------------
    Peak.counts=plyr::ddply(Peaks.Info,plyr::.(Chrom),nrow)

    return(list(Peaks.Info=Peaks.Info,Classification.Info=Classification.Info,
                Peak.counts=Peak.counts))
}
#-------------
#-------------
####################################################
# Functions for the Inference:
####################################################
#-------------
#-------------
#Main significance call function:
SignificanceCall_fun=function(Peaks.Info,PETsData,ChromInf,method,LogFile.dir){
    #global variables for Rcheck:
    Chrom=NULL
    #--------------------------
    #--------------------------------------------
    #---------------Take information for the analysis:
    #--------------------------------------------
    Peaks.Info$Chrom=as.character(Peaks.Info$Chrom)
    PETsData$MainIndex=NULL#remove index not needed
    #--------------------------------------------
    #---------------run inference:
    #--------------------------------------------
    cat("Splitting data by chromosome for inference...\n")
    #split by chromosome all the three inputs:
    # then merge them by chromosome, pass a three way list in c++
    bppass=Get_SignificanceCall_Data_fun(Peaks.Info=Peaks.Info,
                                         PETsData=PETsData,
                                         ChromInf=ChromInf)
    #------------
    #Run inference(call c++):
    #------------
    LogFile1="Computing p-values...\n"
    cat(LogFile1)
    write(LogFile1,file=LogFile.dir,append=TRUE)#write to log file

    InfRes=BiocParallel::bplapply(X=bppass,FUN=PoissonLocalInference_fun,
                                  windows=c(10,15))
    #------------
    #merge and return:
    #------------
    InfRes=do.call(rbind.data.frame,InfRes)
    rownames(InfRes)=NULL
    #------
    #FDR at 0.05:
    #------
    cat("FDR adjusting p-values...\n")
    InfRes$FDRUp=stats::p.adjust(p=InfRes$p.valueUp,method=method)
    InfRes$FDRDown=stats::p.adjust(p=InfRes$p.valueDown,method=method)
    InfRes$FDR=stats::p.adjust(p=InfRes$p.value,method=method)
    #------
    #---return:
    #------
    cat("Inference is done!\n")
    return(InfRes)
}
#-------------
#-------------
# Function for breaking the inputs for significance calling:
Get_SignificanceCall_Data_fun=function(Peaks.Info,PETsData,ChromInf){
    # Rcheck:
    Chrom=NULL
    #-------------
    # break Peaks.Info by chrom:
    #-------------
    Peaks.Info_list=plyr::dlply(Peaks.Info,plyr::.(Chrom),function(x) x)
    #-------------
    # break PETsData by chrom:
    #-------------
    PETsData$MainIndex=NULL#dont need index
    PETsData_list=plyr::dlply(PETsData,plyr::.(Chrom),function(x) x)
    #-------------
    # break ChromInf by chrom:
    #-------------
    ChromInf_list=plyr::dlply(ChromInf,plyr::.(Chrom),function(x) x)
    #-------------
    # create the bppass input:
    #-------------
    bppass=list()
    for(chri in 1:length(Peaks.Info_list)){
        # take Peaks_Info_x:
        Peaks_Info_x=Peaks.Info_list[[chri]]
        Chr_x=Peaks_Info_x$Chrom[1]
        # take PETsData_x:
        PETsData_x=PETsData_list[[which(names(PETsData_list)==Chr_x)]]
        PETsData_x$Chrom=NULL#dont need that anymore
        PETsData_x=as.matrix(PETsData_x)
        # take ChromInf_x:
        ChromInf_x=ChromInf_list[[which(names(ChromInf_list)==Chr_x)]]
        ChromInf_x$Chrom=NULL#dont need that

        bppass[[chri]]=list(Peaks_Info_x=Peaks_Info_x,PETsData_x=PETsData_x,
                          ChromInf_x=ChromInf_x)
    }
    return(bppass)
}
#-------------
#-------------
# Function for Inference:
PoissonLocalInference_fun=function(bppass_x,windows){
    #-------------
    # break data:
    #-------------
    Peaks_Info_x=bppass_x$Peaks_Info_x
    PETsData_x=bppass_x$PETsData_x
    # take UTag
    UTag_x=intervals::Intervals(PETsData_x[,c(1,1)],closed=c(TRUE,TRUE))
    # take DTag
    DTag_x=intervals::Intervals(PETsData_x[,c(2,2)],closed=c(TRUE,TRUE))
    # take chromosome info:
    ChromSize=bppass_x$ChromInf_x$size
    ChromPETs=bppass_x$ChromInf_x$PET.counts
    #-------------
    # ------ Up window 1:
    #-------------
    UpW1=cbind(Peaks_Info_x$CIQ.Up.start-Peaks_Info_x$CIQ.Up.size*windows[1]/2.0,
               Peaks_Info_x$CIQ.Up.start+Peaks_Info_x$CIQ.Up.size*windows[1]/2.0)
    UpW1=intervals::Intervals(UpW1,closed=c(TRUE,TRUE))
    # observed pets in intervals
    UpW1_obs=lengths(intervals::interval_included(UpW1,UTag_x))
    lambdaUpW1=UpW1_obs*Peaks_Info_x$CIQ.Up.size/(intervals::size(UpW1)+1)
    #-------------
    # ------ Up window 2:
    #-------------
    UpW2=cbind(Peaks_Info_x$CIQ.Up.start-Peaks_Info_x$CIQ.Up.size*windows[2]/2.0,
               Peaks_Info_x$CIQ.Up.start+Peaks_Info_x$CIQ.Up.size*windows[2]/2.0)
    UpW2=intervals::Intervals(UpW2,closed=c(TRUE,TRUE))
    # observed pets in intervals
    UpW2_obs=lengths(intervals::interval_included(UpW2,UTag_x))
    lambdaUpW2=UpW2_obs*Peaks_Info_x$CIQ.Up.size/(intervals::size(UpW2)+1)
    #-------------
    # ------ Down window 1:
    #-------------
    DownW1=cbind(Peaks_Info_x$CIQ.Down.start-Peaks_Info_x$CIQ.Down.size*windows[1]/2.0,
                 Peaks_Info_x$CIQ.Down.start+Peaks_Info_x$CIQ.Down.size*windows[1]/2.0)
    DownW1=intervals::Intervals(DownW1,closed=c(TRUE,TRUE))
    # observed pets in intervals
    DownW1_obs=lengths(intervals::interval_included(DownW1,DTag_x))
    lambdaDownW1=DownW1_obs*Peaks_Info_x$CIQ.Down.size/(intervals::size(DownW1)+1)
    #-------------
    # ------ Down window 2:
    #-------------
    DownW2=cbind(Peaks_Info_x$CIQ.Down.start-Peaks_Info_x$CIQ.Down.size*windows[2]/2.0,
                 Peaks_Info_x$CIQ.Down.start+Peaks_Info_x$CIQ.Down.size*windows[2]/2.0)
    DownW2=intervals::Intervals(DownW2,closed=c(TRUE,TRUE))
    # observed pets in intervals
    DownW2_obs=lengths(intervals::interval_included(DownW2,DTag_x))
    lambdaDownW2=DownW2_obs*Peaks_Info_x$CIQ.Down.size/(intervals::size(DownW2)+1)
    #-------------
    # ------ genome backgrounds:
    #-------------
    lambdaBGC_Up=ChromPETs*Peaks_Info_x$CIQ.Up.size/ChromSize;
    lambdaBGC_Down=ChromPETs*Peaks_Info_x$CIQ.Down.size/ChromSize;
    #-------------
    # ------ find maximum:
    #-------------
    OptimalLambda=plyr::llply(1:nrow(Peaks_Info_x),function(i,lambdaUpW1,lambdaUpW2,lambdaBGC_Up,
                                        lambdaDownW1,lambdaDownW2,lambdaBGC_Down){

        UpMax=max(c(2,lambdaUpW1[i],lambdaUpW2[i],lambdaBGC_Up[i]))
        DownMax=max(c(2,lambdaDownW1[i],lambdaDownW2[i],lambdaBGC_Down[i]))

        return(c(UpMax,DownMax))

    },lambdaUpW1=lambdaUpW1,lambdaUpW2=lambdaUpW2,lambdaBGC_Up=lambdaBGC_Up,
    lambdaDownW1=lambdaDownW1,lambdaDownW2=lambdaDownW2,lambdaBGC_Down=lambdaBGC_Down)

    OptimalLambda=do.call(rbind,OptimalLambda)
    #-------------
    # ------ Find folds and p.values:
    #-------------
    # Upstream
    FoldEnrichUp=Peaks_Info_x$Pets/OptimalLambda[,1]
    pvalueUp=stats::ppois(q=Peaks_Info_x$Pets,lambda=OptimalLambda[,1],lower.tail=FALSE)
    # Downstream
    FoldEnrichDown=Peaks_Info_x$Pets/OptimalLambda[,2]
    pvalueDown=stats::ppois(q=Peaks_Info_x$Pets,lambda=OptimalLambda[,2],lower.tail=FALSE)
     # merged p-value:
    pvalue=pvalueUp*pvalueDown
    #-------------
    # ------Return
    #-------------
    Peaks_Info_x$lambdaUp=OptimalLambda[,1]
    Peaks_Info_x$FoldEnrichUp=FoldEnrichUp
    Peaks_Info_x$p.valueUp=pvalueUp
    Peaks_Info_x$lambdaDown=OptimalLambda[,2]
    Peaks_Info_x$FoldEnrichDown=FoldEnrichDown
    Peaks_Info_x$p.valueDown=pvalueDown
    Peaks_Info_x$p.value=pvalue

    return(Peaks_Info_x)

}
#-------------
#-------------
