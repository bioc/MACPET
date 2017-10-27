## ----style,eval=TRUE,echo=FALSE,results='hide'-----------------------------
BiocStyle::latex2

## ----eval=TRUE,echo=TRUE---------------------------------------------------
#Create a temporary test folder, or anywhere you want:
AnalysisDir=file.path(tempdir(),"MACPETtest")
dir.create(AnalysisDir)#where you will save the results.

## --------------------------------------------------------------------------
library(MACPET)

## --------------------------------------------------------------------------
load(system.file("extdata", "pselfData.rda", package = "MACPET"))
class(pselfData) #example name
pselfData #print method

## --------------------------------------------------------------------------
metadata(pselfData)

## --------------------------------------------------------------------------
seqinfo(pselfData)

## --------------------------------------------------------------------------
load(system.file("extdata", "psfitData.rda", package = "MACPET"))
class(psfitData) #example name
psfitData #print method

## --------------------------------------------------------------------------
head(metadata(psfitData)$Peaks.Info)

## --------------------------------------------------------------------------
load(system.file("extdata", "pinterData.rda", package = "MACPET"))
class(pinterData) #example name
pinterData #print method

## ----eval=TRUE-------------------------------------------------------------
metadata(pinterData)

## --------------------------------------------------------------------------
load(system.file("extdata", "pintraData.rda", package = "MACPET"))
class(pintraData)#example name
pintraData#print method

## --------------------------------------------------------------------------
metadata(pintraData)

## --------------------------------------------------------------------------
class(pselfData)
summary(pselfData)

## --------------------------------------------------------------------------
class(psfitData)
summary(psfitData)

## --------------------------------------------------------------------------
class(pintraData)
requireNamespace("ggplot2")
requireNamespace("reshape2")
summary(pintraData,heatmap=TRUE)

## --------------------------------------------------------------------------
class(pinterData)
requireNamespace("ggplot2")
requireNamespace("reshape2")
summary(pinterData,heatmap=TRUE)

## --------------------------------------------------------------------------
requireNamespace("ggplot2")
class(pselfData)
# PET counts plot
plot(pselfData)

## --------------------------------------------------------------------------
class(psfitData)
#binding site couts:
plot(psfitData,kind="PeakCounts")
# region example with binding sites:
plot(psfitData,kind="PeakPETs",RegIndex=1)

## --------------------------------------------------------------------------
class(pintraData)
#plot counts:
plot(pintraData)

## --------------------------------------------------------------------------
class(pinterData)
requireNamespace("igraph")
#network plot:
plot(pinterData)

## ----eval=TRUE,echo=TRUE---------------------------------------------------
class(psfitData)#PSFit class
exportPeaks(object=psfitData,file.out="Peaks",threshold=1e-5,savedir=AnalysisDir)

## ----eval=TRUE,echo=TRUE---------------------------------------------------
class(psfitData)#PSFit class
object=PeaksToGRanges(object=psfitData,threshold=1e-5)
object

## ----eval=TRUE,echo=TRUE---------------------------------------------------
class(psfitData)#PSFit class
TagsToGInteractions(object=psfitData,threshold=1e-5)


## ----eval=TRUE,echo=TRUE---------------------------------------------------
class(psfitData)#PSFit class
PeaksToNarrowPeak(object=psfitData,threshold=1e-5,file.out="MACPET_peaks.narrowPeak",savedir=AnalysisDir)

## ----eval=TRUE,echo=TRUE---------------------------------------------------
 #--remove information and convert to GInteractions:
object=pselfData
S4Vectors::metadata(object)=list(NULL)
class(object)="GInteractions"
GenomePkg="BSgenome.Hsapiens.UCSC.hg19" #genome of the data.
BlackList=TRUE
object=ConvertToPSelf(object=object,GenomePkg=GenomePkg,BlackList=BlackList)
class(object)

## ----echo=TRUE,eval=TRUE---------------------------------------------------
AnalysisStatistics(x.self=pselfData,#One of the self-ligated classes.
                    x.intra=pintraData,#NULL for not printing the class.
                    x.inter=pinterData,#NULL for not printing the class.
                    #specify a name for the output to be saved.
                    file.out="Statistics",
                    savedir=AnalysisDir,#Where to save the output.
                    threshold=1e-5)#Theshold for FDR


## ----echo=TRUE,eval=TRUE---------------------------------------------------
 #load sample data to use in the algorithm and give inputs:
 DataDir=system.file("extdata", package = "MACPET") #data path
 DataFile="SampleChIAPETData.bam" #data name
 PopImage=TRUE #Sample data, not very good classifiction results.
 GenomePkg="BSgenome.Hsapiens.UCSC.hg19" #genome of the data.
 fileSelf="pselfData" #name for Self-ligated
 fileIntra="pintraData" #name for Intra-chromosomal
 fileInter="pinterData" #name for Inter-chromosomal
 BlackList=TRUE #remove PETs in black listed regions
 fileSelfFit="psfitData"
 method="BH"
 Stages=c(0:1)

#parallel backhead can be created using the BiocParallel package
# requireNamespace("BiocParallel")
# snow <- BiocParallel::SnowParam(workers = 1, type = "SOCK", progressbar=FALSE)
# BiocParallel::register(snow, default=TRUE)

#-run for the whole binding site analysis:
PeakCallerUlt(DataDir=DataDir,
           DataFile=DataFile,
           AnalysisDir=AnalysisDir,
           GenomePkg=GenomePkg,
           fileSelf=fileSelf,
           fileIntra=fileIntra,
           fileInter=fileInter,
           PopImage=PopImage,
           fileSelfFit=fileSelfFit,
           method=method,
           BlackList=BlackList,
           Stages=Stages)
#load results:
load(file.path(AnalysisDir,fileSelf))
class(pselfData) # see methods for this class
load(file.path(AnalysisDir,fileIntra))
class(pintraData) # see methods for this class
load(file.path(AnalysisDir,fileInter))
class(pinterData) # see methods for this class
load(file.path(AnalysisDir,fileSelfFit))
class(psfitData) # see methods for this class

#-----delete test directory:
unlink(AnalysisDir,recursive=TRUE)

