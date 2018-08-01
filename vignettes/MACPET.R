## ----style,eval=TRUE,echo=FALSE,results='hide'-----------------------------
BiocStyle::latex

## ----eval=TRUE,echo=TRUE---------------------------------------------------
#Create a temporary test folder, or anywhere you want:
SA_AnalysisDir=file.path(tempdir(),"MACPETtest")
dir.create(SA_AnalysisDir)#where you will save the results.

## --------------------------------------------------------------------------
library(MACPET)

## --------------------------------------------------------------------------
load(system.file("extdata", "MACPET_pselfData.rda", package = "MACPET"))
class(MACPET_pselfData) #example name
MACPET_pselfData #print method

## --------------------------------------------------------------------------
metadata(MACPET_pselfData)

## --------------------------------------------------------------------------
seqinfo(MACPET_pselfData)

## --------------------------------------------------------------------------
load(system.file("extdata", "MACPET_psfitData.rda", package = "MACPET"))
class(MACPET_psfitData) #example name
MACPET_psfitData #print method

## --------------------------------------------------------------------------
head(metadata(MACPET_psfitData)$Peaks.Info)

## --------------------------------------------------------------------------
load(system.file("extdata", "MACPET_pinterData.rda", package = "MACPET"))
class(MACPET_pinterData) #example name
MACPET_pinterData #print method

## ----eval=TRUE-------------------------------------------------------------
metadata(MACPET_pinterData)

## --------------------------------------------------------------------------
load(system.file("extdata", "MACPET_pintraData.rda", package = "MACPET"))
class(MACPET_pintraData)#example name
MACPET_pintraData#print method

## --------------------------------------------------------------------------
metadata(MACPET_pintraData)

## --------------------------------------------------------------------------
load(system.file("extdata", "MACPET_GenomeMapData.rda", package = "MACPET"))
class(MACPET_GenomeMapData) #example name
MACPET_GenomeMapData #print method

## --------------------------------------------------------------------------
metadata(MACPET_GenomeMapData)

## --------------------------------------------------------------------------
class(MACPET_pselfData)
summary(MACPET_pselfData)

## --------------------------------------------------------------------------
class(MACPET_psfitData)
summary(MACPET_psfitData)

## --------------------------------------------------------------------------
class(MACPET_pintraData)
requireNamespace("ggplot2")
requireNamespace("reshape2")
summary(MACPET_pintraData,heatmap=TRUE)

## --------------------------------------------------------------------------
class(MACPET_pinterData)
requireNamespace("ggplot2")
requireNamespace("reshape2")
summary(MACPET_pinterData,heatmap=TRUE)

## --------------------------------------------------------------------------
class(MACPET_GenomeMapData)
summary(MACPET_GenomeMapData)

## --------------------------------------------------------------------------
requireNamespace("ggplot2")
class(MACPET_pselfData)
# PET counts plot
plot(MACPET_pselfData)

## --------------------------------------------------------------------------
class(MACPET_psfitData)
#binding site couts:
plot(MACPET_psfitData,kind="PeakCounts")

## --------------------------------------------------------------------------
# region example with binding sites:
plot(MACPET_psfitData,kind="PeakPETs",RegIndex=1)

## --------------------------------------------------------------------------
class(MACPET_pintraData)
#plot counts:
plot(MACPET_pintraData)

## --------------------------------------------------------------------------
class(MACPET_pinterData)
requireNamespace("igraph")
#network plot:
plot(MACPET_pinterData)

## --------------------------------------------------------------------------
class(MACPET_GenomeMapData)
requireNamespace("igraph")
#network plot:
plot(MACPET_GenomeMapData,Type='network-circle')

## ----eval=TRUE,echo=TRUE---------------------------------------------------
class(MACPET_psfitData)#PSFit class
exportPeaks(object=MACPET_psfitData,file.out="Peaks",threshold=1e-5,savedir=SA_AnalysisDir)

## ----eval=TRUE,echo=TRUE---------------------------------------------------
class(MACPET_psfitData)#PSFit class
object=PeaksToGRanges(object=MACPET_psfitData,threshold=1e-5)
object

## ----eval=TRUE,echo=TRUE---------------------------------------------------
class(MACPET_psfitData)#PSFit class
TagsToGInteractions(object=MACPET_psfitData,threshold=1e-5)


## ----eval=TRUE,echo=TRUE---------------------------------------------------
class(MACPET_psfitData)#PSFit class
PeaksToNarrowPeak(object=MACPET_psfitData,threshold=1e-5,
                  file.out="MACPET_peaks.narrowPeak",savedir=SA_AnalysisDir)

## ----eval=TRUE,echo=TRUE---------------------------------------------------
 #--remove information and convert to GInteractions:
object=MACPET_pselfData
#--remove information and convert to GInteractions:
S4Vectors::metadata(object)=list(NULL)
class(object)='GInteractions'
#----input parameters
S2_BlackList=TRUE
SA_prefix="MACPET"
S2_AnalysisDir=SA_AnalysisDir

ConvertToPSelf(object=object,
               S2_BlackList=S2_BlackList,
               SA_prefix=SA_prefix,
               S2_AnalysisDir=S2_AnalysisDir)
#load object:
rm(MACPET_pselfData)#old object
load(file.path(S2_AnalysisDir,"MACPET_pselfData"))
class(MACPET_pselfData)


## ----eval=TRUE,echo=TRUE---------------------------------------------------
class(MACPET_GenomeMapData)#GenomeMap class
GetSignInteractions(object=MACPET_GenomeMapData,
                     threshold = NULL,
                     ReturnedAs='GInteractions')

## ----echo=TRUE,eval=TRUE---------------------------------------------------
AnalysisStatistics(x.self=MACPET_psfitData,
                   x.intra=MACPET_pintraData,
                   x.inter=MACPET_pinterData,
                   file.out='AnalysisStats',
                   savedir=SA_AnalysisDir,
                   threshold=1e-5)


## ----echo=TRUE,eval=TRUE---------------------------------------------------
requireNamespace('ggplot2')

#Create a temporary forder, or anywhere you want:
S1_AnalysisDir=SA_AnalysisDir

#directories of the BAM files:
BAM_file_1=system.file('extdata', 'SampleChIAPETDataRead_1.bam', package = 'MACPET')
BAM_file_2=system.file('extdata', 'SampleChIAPETDataRead_2.bam', package = 'MACPET')
SA_prefix="MACPET"

#convert to paired-end BAM:
ConvertToPE_BAM(S1_AnalysisDir=S1_AnalysisDir,
             SA_prefix=SA_prefix,
             S1_BAMStream=2000000,S1_image=TRUE,
             S1_genome="hg19",BAM_file_1=BAM_file_1,
             BAM_file_2=BAM_file_2)

#test if the resulted BAM is paired-end:
PairedBAM=file.path(S1_AnalysisDir,paste(SA_prefix,"_Paired_end.bam",sep=""))
Rsamtools::testPairedEndBam(file = PairedBAM, index = PairedBAM)

bamfile = Rsamtools::BamFile(file = PairedBAM,asMates = TRUE)
GenomicAlignments::readGAlignmentPairs(file = bamfile,use.names = FALSE,
                                with.which_label = FALSE,
                                strandMode = 1)


## ----echo=TRUE,eval=TRUE---------------------------------------------------

#give directory of the BAM file:
S2_PairedEndBAMpath=system.file('extdata', 'SampleChIAPETData.bam', package = 'MACPET')

#give prefix name:
SA_prefix="MACPET"

#parallel backhead can be created using the BiocParallel package
#parallel backhead can be created using the BiocParallel package
#requireNamespace('BiocParallel')
#snow <- BiocParallel::SnowParam(workers = 4, type = 'SOCK', progressbar=FALSE)
#BiocParallel::register(snow, default=TRUE)

# packages for plotting:
requireNamespace('ggplot2')

#-run for the whole binding site analysis:
MACPETUlt(SA_AnalysisDir=SA_AnalysisDir,
       SA_stages=c(2:4),
       SA_prefix=SA_prefix,
       S2_PairedEndBAMpath=S2_PairedEndBAMpath,
       S2_image=TRUE,
       S2_BlackList=TRUE,
       S3_image=TRUE,
       S4_image=TRUE,
       S4_FDR_peak=1)# the data is small so use all the peaks found.

#load results:
SelfObject=paste(SA_prefix,"_pselfData",sep="")
load(file.path(SA_AnalysisDir,"S2_results",SelfObject))
SelfObject=get(SelfObject)
class(SelfObject) # see methods for this class

IntraObject=paste(SA_prefix,"_pintraData",sep="")
load(file.path(SA_AnalysisDir,"S2_results",IntraObject))
IntraObject=get(IntraObject)
class(IntraObject) # see methods for this class

InterObject=paste(SA_prefix,"_pinterData",sep="")
load(file.path(SA_AnalysisDir,"S2_results",InterObject))
InterObject=get(InterObject)
class(InterObject) # see methods for this class

SelfFitObject=paste(SA_prefix,"_psfitData",sep="")
load(file.path(SA_AnalysisDir,"S3_results",SelfFitObject))
SelfFitObject=get(SelfFitObject)
class(SelfFitObject) # see methods for this class

GenomeMapObject=paste(SA_prefix,"_GenomeMapData",sep="")
load(file.path(SA_AnalysisDir,"S4_results",GenomeMapObject))
GenomeMapObject=get(GenomeMapObject)
class(GenomeMapObject) # see methods for this class

#-----delete test directory:
unlink(SA_AnalysisDir,recursive=TRUE)


