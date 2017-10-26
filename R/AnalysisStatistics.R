#' @title Count Statistics for ChIA-PET.
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#' @references \insertRef{macpetcite}{MACPET}
#' @export
#' @description \code{AnalysisStatistics} prints and saves count statistics
#' for the current inputs of the peak-calling analysis.
#--#' function parameters:
#' @param x.self An object of class \code{\linkS4class{PSelf}} or
#'  \code{\linkS4class{PSFit}}.
#' @param x.intra An object of class \code{\linkS4class{PIntra}}
#' (optional).
#' @param x.inter An object of class \code{\linkS4class{PInter}}
#' (optional).
#' @param file.out A string with the name of the output to be saved to
#' \code{savedir}. If \code{NULL} the function will only print the output.
#' @param threshold A numeric indicating the FDR cut-off, used when
#' \code{class(x.self)=\linkS4class{PSFit}}. If NULL, no threshold is applied.
#' @param savedir A string with the directory to save the ouput file. If
#' \code{NULL} then the function will only print the output.
#--#' what the function returns:

#' @return Based on the inputs, \code{AnalysisStatistics} prints the total
#'  Self-ligated, Intra- and Inter-chromosomal PETs, as well as the total regions,
#'   total candidate peaks and total significant
#' peaks (if \code{threshold!=NULL} and
#' \code{class(x.self)=\linkS4class{PSFit}}). If \code{file.out} and
#' \code{savedir} are not \code{NULL} then it also saves the output to a csv
#' file in \code{savedir}.
#'
#' @seealso \code{\linkS4class{PSelf}},
#'  \code{\linkS4class{PSFit}},\code{\linkS4class{PIntra}},
#' \code{\linkS4class{PInter}}
#'
#--#' import packages to NAMESPACE:
#' @importFrom S4Vectors metadata
#' @importFrom plyr ddply .
#' @importFrom GenomeInfoDb genome
#' @importFrom knitr kable
#' @importFrom utils write.table
#'
#' @examples
#' #Create a test forder on the desktop, or anywhere you want:
#' savedir=file.path(path.expand('~'),'Desktop')
#' dir.create(file.path(savedir,"MACPET.test"))
#' savedir=file.path(savedir,"MACPET.test")#where you will save
#' #the results
#'
#' #load Inter-chromosomal data:
#' load(system.file("extdata", "pinterData.rda", package = "MACPET"))
#' class(pinterData)
#'
#' #load Intra-chromosomal data:
#' load(system.file("extdata", "pintraData.rda", package = "MACPET"))
#' class(pintraData)
#'
#' #################################################################
#' #load Self-ligated data: (class=PSelf)
#' load(system.file("extdata", "pselfData.rda", package = "MACPET"))
#' class(pselfData)
#'
#' #Print analysis:
#' AnalysisStatistics(x.self=pselfData,
#'                    x.intra=pintraData,
#'                    x.inter=pinterData,
#'                    file.out="AnalysisStats",
#'                    savedir=savedir)
#'
#' #################################################################
#' #load Self-ligated data: (class=PSFit)
#' load(system.file("extdata", "psfitData.rda", package = "MACPET"))
#' class(psfitData)
#'
#' #Print analysis:
#' AnalysisStatistics(x.self=psfitData,
#'                    x.intra=pintraData,
#'                    x.inter=pinterData,
#'                    file.out="AnalysisStats",
#'                    savedir=savedir,
#'                    threshold=1e-5)
#'
#' #-----delete test directory:
#' unlink(savedir,recursive=TRUE)


AnalysisStatistics=function(x.self,x.intra=NULL,x.inter=NULL,file.out=NULL,
                            threshold=1e-5,savedir=NULL){
    #global variables for Rcheck:
    FDR=Chrom=NULL
    #--------------------------
    #input Check:
    if(!class(x.self)%in%c("PSelf","PSFit")){
        stop("x.self has to be on of the following classes:
             PSelf or PSFit")
    }else if(class(x.self)=="PSFit"&!is.numeric(threshold)){
        threshold=NULL
    }
    #take intra:
    if(!is.null(x.intra)){
        if(class(x.intra)!="PIntra"){
            stop("x.intra has to be PIntra class!")
        }else{
            RES.intra=S4Vectors::metadata(x.intra)$InteractionCounts
            Nintra=sum(RES.intra$Counts)
        }
    }
    #take inter:
    if(!is.null(x.inter)){
        if(class(x.inter)!="PInter"){
            stop("x.inter has to be PInter class!")
        }else{
            RES.inter=S4Vectors::metadata(x.inter)$InteractionCounts
            RES.inter=colSums(RES.inter)
            RES.inter=data.frame(Chrom=names(RES.inter),
                                 Counts=as.numeric(RES.inter))
            Ninter=sum(RES.inter$Counts)
        }

    }
    #take self:
    Nself=length(x.self)
    SLmean=S4Vectors::metadata(x.self)$SLmean
    MaxSize=S4Vectors::metadata(x.self)$MaxSize
    MinSize=S4Vectors::metadata(x.self)$MinSize
    hg=unique(GenomeInfoDb::genome(GenomeInfoDb::seqinfo(x.self)))
    Organism=as.character(S4Vectors::metadata(x.self)$GenInfo$organism)
    Info=data.frame(SLmean=SLmean,hg=hg,Organism=Organism,
                    SelfBord=paste(MinSize,"/",MaxSize," bp",
                                   sep=""),Total.PETs=Nself)
    RES=S4Vectors::metadata(x.self)$Self_info
    if(class(x.self)=="PSelf"){
        colnames(RES)=c("Chrom","Self")
        colnames(Info)=c("Self-lig. mean size","Genome","Organism",
                         "Self Borders","Tot. Self")

    }else if(class(x.self)=="PSFit"&is.null(threshold)){
        colnames(RES)=c("Chrom","Self","Regions",
                        "Peaks")
        Info$Total.regions=sum(S4Vectors::metadata(x.self)$
                                   Self_info$Region.counts)
        Info$Total.bs=sum(S4Vectors::metadata(x.self)$Self_info$Peak.counts)
        colnames(Info)=c("Self-lig. mean size","Genome","Organism",
                         "Self Borders","Tot. Self","Regions","Peaks")
    }else if(class(x.self)=="PSFit"&!is.null(threshold)){
        colnames(RES)=c("Chrom","Self","Regions",
                        "Peaks")
        Info$Total.regions=sum(S4Vectors::metadata(x.self)$Self_info$
                                   Region.counts)
        Info$Total.bs=sum(S4Vectors::metadata(x.self)$Self_info$Peak.counts)
        Sig.bs=S4Vectors::metadata(x.self)$Peaks.Info
        Sig.bs=subset(Sig.bs,FDR<threshold)
        Sig.bs=plyr::ddply(Sig.bs,plyr::.(Chrom),nrow)
        RES$`Sign. Peaks`=0
        RES$`Sign. Peaks`[match(Sig.bs$Chrom,RES$Chrom)]=
            Sig.bs$V1
        Info$Total.sig=sum(RES$`Sign. Peaks`)
        colnames(Info)=c("Self-lig. mean size","Genome","Organism",
                         "Self Borders","Tot. Self","Regions",
                         "Peaks","Sign. Peaks")
    }
    #append:
    if(!is.null(x.intra)){
        RES[match(RES.intra$Chrom,RES$Chrom),
            c("Intra")]=RES.intra$Counts
        if(any(is.na(RES$Intra))){
            RES$Intra[which(is.na(RES$Intra))]=0
        }
        Info[,c("Tot. Intra")]=Nintra
    }

    if(!is.null(x.inter)){
        RES[match(RES.inter$Chrom,RES$Chrom),
            c("Inter")]=RES.inter$Counts
        if(any(is.na(RES$Inter))){
            RES$Inter[which(is.na(RES$Inter))]=0
        }
        Info[,c("Tot. Inter")]=Ninter
    }
    #print:
    cat(paste(rep("-",25),collapse=""))
    cat("\n PETs Counts Summary \n")
    cat(paste(rep("-",25),collapse=""))
    print(knitr::kable(RES,row.names = FALSE,
                       align=c("c"),format="markdown",padding=1,output=TRUE))

    #print:
    print(knitr::kable(Info,row.names = FALSE,
                       align=c("c"),format="rst",padding=1,output=TRUE))
    #save
    if(!is.null(file.out)&!is.null(savedir)){
        if(!dir.exists(savedir)){
            stop("savedir does not exist!",call. = FALSE)
        }
        Desc="#Count statistics for each chromosome:"
        writeLines(paste(Desc,sep = "\n"),
                   file.path(savedir,paste(file.out,".csv",sep="")))

        suppressWarnings(
            utils::write.table(RES,quote=FALSE,
                               file=file.path(savedir,
                                              paste(file.out,".csv",sep="")),
                                            sep=";",
                                            col.names=TRUE,row.names=FALSE,
                                            qmethod="double",append=TRUE))

        return(cat("The output has been saved at the savedir"))
    }
}
