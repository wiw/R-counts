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
	outputScttr <- "scatter_plots"	# output folder for scatter plots in WD
	startCol <- 7	# the number of last column in GATCs file, default "7"
	gatcFile <- paste(workDir, "/GATCs.txt", sep="")	# location you GATCs file
	samplesListFile <- paste(workDir, "/source.csv", sep="")	# location you file with source data
	writeTemp <- T		# use this option if you need to get the intermediate files, default "T"
	dir.create(file.path(workDir, prefixDir), showWarnings = FALSE)
dir.create(file.path(workDir, prefixDir, outputScttr), showWarnings = FALSE)

# Whether a script run earlier?
###############################
setwd(file.path(workDir, prefixDir))
	alreadyRun <- list.files(prefixDir, "^.*Step_01.*")
setwd(workDir)
############################################
################ FUNCTIONS #################
############################################

# Write intermediate files function
################################### 
	WriteIntermediateFiles <- function(source, output.file) {
		if (writeTemp == T) {
			write.table(source, file=file.path(prefixDir, output.file), sep=";", row.names=F, col.names=T, quote=F, dec=",", append=F)
		}
	}

# Scatter Plots on Averaged data function
#########################################
ScatterPlottingOnAveraged <- function(dataSet, tag) {

	for (inner in 8:ncol(DATA)) { 
		if ( inner %% 2 != 0 ) {
			edge <- inner - 1
				Cor.P <- round(cor(dataSet[, inner], dataSet[, edge], method="pearson", use="pairwise.complete.obs"), digits=2)
				Cor.S <- round(cor(dataSet[, inner], dataSet[, edge], method="spearman", use="pairwise.complete.obs"), digits=2)
				rowCol <- grep(sub("([0-9a-zA-Z\\._]+)(_(edge|inner))(.*)", "\\1\\4", names(DATA)[inner], perl=T), stat$Data.set)
				stat[rowCol, 12] <- Cor.P
				stat[rowCol, 13] <- Cor.S
				bmp(filename=file.path(prefixDir, outputScttr, paste("scatter_", names(dataSet)[inner], "_vs_", names(dataSet)[edge], "_", tag, ".bmp", sep="")), width=600, height=600, units = "px")
				Cor.P <- round(cor(dataSet[, inner], dataSet[, edge], method="pearson", use="pairwise.complete.obs"), digits=2)
				Cor.S <- round(cor(dataSet[, inner], dataSet[, edge], method="spearman", use="pairwise.complete.obs"), digits=2)
				titleName <- sub("_inner_local_GATCcounts", "", basename(file_path_sans_ext(samplesList[grep(names(DATA)[inner], samplesList$id), 6])))

				plot(x=dataSet[, inner], y=dataSet[, edge], cex=0.3, xlab=names(dataSet[inner]), ylab=names(dataSet[edge]), main=titleName, text(x=min(dataSet[, inner], na.rm=T) + 0.5, y=max(dataSet[, edge], na.rm=T) - 0.5, labels=c(paste("r = ", Cor.P, "\n\n", sep=""), paste("s = ", Cor.S, sep=""))))
				dev.off()
		}

	}

	WriteIntermediateFiles(source=stat, output.file="statistics.csv")
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
gatcs <- cbind(gatcs, matrix(data=NA, nrow=nrow(gatcs), ncol=nrow(samplesList)))
	for (i in 1:nrow(samplesList)){
		colnames(gatcs)[startCol+i] <- samplesList$id[i]
			load(file=samplesList$path[i])
			if (all(gatcs$ID.il == reads2GATC$ID)) gatcs[, startCol + i] <- reads2GATC$count
	}
rm(i)

# Remove unimportant data
#########################
	DATA <- gatcs[!(gatcs$chr %in% c("U", "M", "Uextra" )), ]
# Scatter Plots on data
#######################
	stat <- read.csv("cutadapt_statistics.csv", header=T, sep=";", stringsAsFactors=F)
stat <- cbind(stat, as.data.frame(matrix(data=NA, nrow=nrow(stat), ncol=2)))
	names(stat)[12:13] <- c("pearson", "spearman")
	print("Run scatter")
	ScatterPlottingOnAveraged(dataSet=DATA, tag="25feb")
