#!/usr/bin/R
########################################################################
# Alexey Pindyurin, Anton Ivankin, September 12, 2014, DAM_count_statistics.R
#
# DESCRIPTION:
#   
#
# DATA:
#   input data is specified in the samples list file which is supplied as
#   argument to this runner script. File with samples to generate mannualy.
#
# OUTPUT:
#   
#
# VERSIONS:
#   140912: First revision!
#     
########################################################################
rm(list=ls())
library(gplots)
library(tools)
# Declare variables
###################
prefixDir <- "edge_inner" # directory for other experiments
workDir <- getwd()	# working directory (WD)
outputGff <- "gff"	# output folder for gff in WD
outputWig <- "wig"	# output folder for wig in WD
outputScttr <- "scatter_plots"	# output folder for scatter plots in WD
startCol <- 7	# the number of last column in GATCs file, default "7"
gatcFile <- paste(workDir, "/GATCs.txt", sep="")	# location you GATCs file
samplesListFile <- paste(workDir, "/source.csv", sep="")	# location you file with source data
needCombine <- F	# are you need to combine some columns into one (T or F)? we recommend set the "F"; if you select "T" - edit the "combine" vector on string #85 into this code
usePseudoCounts <- T	# are you need to add pseudo counts into source data (T or F)? Default "T"
pseudoCounts <- c(0.01)		# the vector of pseudo counts, default "c(1)"
corrMethod <- c("pearson", "spearman")	# the vector of type correlation method - please don't change this setting
heatmapColors <- greenred(200)	# check you color into heatmap
labelHeatmap <- c("B", "C")	# the vector of heatmap label - please don't change this setting
labelAcf <- c("A_ALL", "B_MA")	# the vector of acf label - please don't change this setting
writeTemp <- T		# use this option if you need to get the intermediate files, default "T"
needSomeFiles <- F	# if need calculate only some files from source list, default "F"
someFiles <- c(9:16)	# the region of files which need calculate
dir.create(file.path(workDir, prefixDir), showWarnings = FALSE)
dir.create(file.path(workDir, prefixDir, outputGff), showWarnings = FALSE)
dir.create(file.path(workDir, prefixDir, outputWig), showWarnings = FALSE)
dir.create(file.path(workDir, prefixDir, outputScttr), showWarnings = FALSE)

# Whether a script run earlier?
###############################
setwd(file.path(workDir, prefixDir))
alreadyRun <- list.files(prefixDir, "^.*Step_01.*")
if (length(alreadyRun) == 1) {
startCol <- 0
}
setwd(workDir)
############################################
################ FUNCTIONS #################
############################################

# Write intermediate files function
################################### 
WriteIntermediateFiles <- function(source, output.file) {
	if (writeTemp == T) {
		write.table(source, file=file.path(prefixDir, output.file), sep=";", row.names=F, col.names=T, quote=F, dec=".", append=F)
	}
}

# Scatter Plots on Averaged data function
#########################################
ScatterPlottingOnAveraged <- function(dataSet, tag) {
  for (j in 8:(ncol(dataSet))){
  	for (i in 8:(ncol(dataSet))){
  		if (j != i) {
			  bmp(filename=file.path(prefixDir, outputScttr, paste("scatter_", names(dataSet)[j], "_vs_", names(dataSet)[i], "_", "_", tag, ".bmp", sep="")), width=600, height=600, units = "px")
				Cor.P <- round(cor(dataSet[, j], dataSet[, i], method="pearson", use="pairwise.complete.obs"), digits=2)
				Cor.S <- round(cor(dataSet[, j], dataSet[, i], method="spearman", use="pairwise.complete.obs"), digits=2)
				titleName <- basename(file_path_sans_ext(samplesList[grep(names(DATA)[j], samplesList$id), 6]))
				plot(x=dataSet[, j], y=dataSet[, i], cex=0.3, xlab=names(dataSet[j]), ylab=names(dataSet[i]), main=titleName, text(x=min(dataSet[, j], na.rm=T) + 0.5, y=max(dataSet[, i], na.rm=T) - 0.5, labels=c(paste("r = ", Cor.P, "\n\n", sep=""), paste("s = ", Cor.S, sep=""))))
     		rm(Cor.P)
     		rm(Cor.S)
     		dev.off()
     	}
		}
	}
}

#######################################
################# END #################
#######################################
print("Run script")
# Load GATC counts in data frame
################################
if (startCol == 0) {
    step01 <- read.delim(alreadyRun, header=T, as.is=T, dec=".")
    startCol <- ncol(step01)
    gatcs <- step01
    } else {
        gatcs <- read.delim(gatcFile, header=T, as.is=T, dec=".")
      }
  samplesList <- read.delim(file=samplesListFile, header=T, dec=".", stringsAsFactors=F, as.is=T)
	if (needSomeFiles == T) {
	samplesList <- samplesList[someFiles, ]
	}
  gatcs <- cbind(gatcs, matrix(data=NA, nrow=nrow(gatcs), ncol=nrow(samplesList)))
    for (i in 1:nrow(samplesList)){
    colnames(gatcs)[startCol+i] <- samplesList$id[i]
    load(file=samplesList$path[i])
    if (all(gatcs$ID.il == reads2GATC$ID)) gatcs[, startCol + i] <- reads2GATC$count
    }
  rm(i)
  currentDate <- format(Sys.time(), "%d-%m-%Y")
  load.gatc.df <- paste("DF_Counts_Step_01_Raw_Counts_", currentDate, ".csv", sep="")
	WriteIntermediateFiles(source=gatcs, output.file=load.gatc.df)

# Remove unimportant data
#########################
  DATA <- gatcs[!(gatcs$chr %in% c("U", "M", "Uextra" )), ]
  use.chr.only <- paste("DF_Counts_Step_02_Useful_Chrs_Only_", currentDate, ".csv", sep="")
	WriteIntermediateFiles(source=DATA, output.file=use.chr.only)
# Scatter Plots on data
#######################
ScatterPlottingOnAveraged(dataSet=DATA, tag="")
