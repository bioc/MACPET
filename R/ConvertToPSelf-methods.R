#' @title Convert GInteraction object to PSelf object
#' @author Ioannis Vardaxis, \email{ioannis.vardaxis@@ntnu.no}
#' @references \insertRef{macpetcite}{MACPET}
#'
#' @description \code{ConvertToPSelf} converts a GInteractions object to
#' class to \code{\linkS4class{PSelf}} object.
#'
#' @details \code{\link{PeakCallerUlt}} with Stage=0 is the first function the user can
#' use for separating Inter-chromosomal, Intra-chromosomal and Self-ligated PETs by giving
#' the whole paired-end BAM/SAM file as input. However the user might only have
#' Self-ligated data available and already separated from the Inter/Intra-chromosomal
#' PETs. \code{ConvertToPSelf} can then be used in the Self-ligated data to convert
#' a \code{\link[InteractionSet]{GInteractions}} object containing only the Self-ligated
#' PETs to a \code{\linkS4class{PSelf}} class for further analysis.
#'
#' @param object An object of \code{\link[InteractionSet]{GInteractions}} class.
#' @param GenomePkg See  \code{\link{PeakCallerUlt}}.
#' @param ... not used.
#' @param BlackList See  \code{\link{PeakCallerUlt}}.
#'
#' @seealso \code{\linkS4class{PSelf}}
#---#'define default:
#' @name ConvertToPSelf
#' @export
#' @include AllClasses.R
#'
#'
#' @examples
#' #load Self-ligated data: (class=PSelf)
#' load(system.file("extdata", "pselfData.rda", package = "MACPET"))
#' class(pselfData)
#'
#' object=pselfData
#' #--remove information and convert to GInteractions:
#' S4Vectors::metadata(object)=list(NULL)
#' class(object)="GInteractions"
#' GenomePkg="BSgenome.Hsapiens.UCSC.hg19" #genome of the data.
#' BlackList=TRUE
#'
#' object=ConvertToPSelf(object=object,GenomePkg=GenomePkg,BlackList=BlackList)
#' class(object)
#'

#default:
ConvertToPSelf <- function(object,...){
    UseMethod("ConvertToPSelf",object=object)
}

#' @rdname ConvertToPSelf
#' @method ConvertToPSelf default
#' @export
ConvertToPSelf.default = function(object,...) {
    stop(paste("No ConvertToPSelf method for class ",class(object),sep=""))
}

#' @rdname ConvertToPSelf
#' @method ConvertToPSelf GInteractions
#' @return An object of class \code{\linkS4class{PSelf}}.
#' @export
ConvertToPSelf.GInteractions=function(object,GenomePkg,BlackList,...){
    #----R-check:
    pkgname=NULL
    #----
    #-------------------------
    #---------Test the genome:
    #-------------------------
    if(!is.character(GenomePkg)){
        stop("GenomePkg has to be of chracter class!",call. = FALSE)
    }else if(!GenomePkg%in%BSgenome::available.genomes()){
        stop("GenomePkg is not a part of the available.genomes! See ??BSgenome::available.genomes",call. = FALSE)
    }else if(!GenomePkg%in%BSgenome::installed.genomes()){
        stop("GenomePkg is not installed! See ??BSgenome::installed.genomes",call. = FALSE)
    }else{
        GenInfo=subset(BSgenome::installed.genomes(splitNameParts=TRUE),
                       pkgname==GenomePkg)
        ChromLengths=BSgenome::getBSgenome(GenomePkg)
        ChromLengths=GenomeInfoDb::seqlengths(ChromLengths)
        ChromLengths=data.frame(chrom=names(ChromLengths),size=ChromLengths)
        rownames(ChromLengths)=NULL
    }
    #--------------------------
    #---------load black list:
    #--------------------------
    if(is.data.frame(BlackList)){
        if(!c("Chrom","Region.Start","Region.End")%in%colnames(BlackList)){
            stop("Give correct colnames to BlackList if it is given as data.frame!",call. = FALSE)
        }

    }else if(!is.logical(BlackList)){
        BlackList=NULL

    }else if(BlackList){
        hg=GenInfo$provider_version
        if(hg=="hg19"){
            BlackList=sysdata$Black_list_hg19
        }else{
            BlackList=NULL
        }
    }else{
        BlackList=NULL
    }
    #-------------------------
    #---------Trim anchors, add info etc:
    #-------------------------
    ArgPairedData=list(ChromLengths=ChromLengths,GenInfo=GenInfo,
                       BlackList=BlackList,LogFile.dir=NULL)
    object=GInteractionsCovnert_fun(ChIAPET=object,ArgPairedData=ArgPairedData)
    Nobject=length(object)#before removing any
    cat(paste("Total PETs found in data:",Nobject),"\n")
    #-------------------------
    #---------Check if Inter-chromosomal:
    #-------------------------
    object=subset(object,!is.na(object$Dist))
    Nreduced=length(object)
    if(Nobject!=Nreduced){
        cat(paste("Total Inter-Chromsomal removed from data:",
                  Nobject-Nreduced),"\n")
    }
    if(Nreduced==0){
        #only Inter-chromosomal:
        stop("object contained only Inter-chromosomal PETs!")
    }
    #-------------------------
    #-----Convert to PSelf class:
    #-------------------------
    SelfIndicator=1:Nreduced#all the data
    object=FindSelf_fun(x=object,SelfIndicator=SelfIndicator,
                        ArgFindSelf=ArgPairedData,PSelf.Convert=TRUE)
    return(object)
}

