#!/usr/bin/R
	rm(list=ls())
	library(ggplot2)

# Declare variables
###################
	prefixDir <- "RUN22-06-2015_PreTest" # output directory into working directory
	workDir <- getwd()	# working directory (WD)
	sourceDir <- "/home/anton/backup/output/RUN22-06-2015" # location your RData files. You can specify the highest folder as it is possible. Searching runs recursively.
	damIdLocation <- "/home/anton/data/DAM/RUN/damid_description.csv" # location your DamID-Description file
	startCol <- 7	# the number of last column in GATCs file, default "7"
	outputPreliminaryStat <- "pre_stat"
	gatcFile <- paste(workDir, "GATCs_mod.txt", sep="/")	# location you GATCs file

# Create folders
################
	dir.create(file.path(workDir, prefixDir), showWarnings = FALSE)
	dir.create(file.path(workDir, prefixDir, outputPreliminaryStat), showWarnings = FALSE)
# Whether a script run earlier?
###############################
setwd(workDir)

############################################
################ FUNCTIONS #################
############################################

# Make samples list file
##########################
MakeSamplesListFile <- function(SOURCE, DAMID) {
filePath <- list.files(path=SOURCE, pattern="*_local_GATCcounts.RData", full.names=T, recursive=T)
baseFile <- unique(sub("(.*)_(edge|inner).*", "\\1", basename(filePath), perl=T))
damIdDscrp <- read.delim(DAMID, header=T, sep="\t", stringsAsFactors=F)
for (i in baseFile) {
	if (length(grep("paired", i)) != 0){
	clearFileName <- sub("(.*)_paired", "\\1", i, perl=T)
	fileID <- subset(damIdDscrp, grepl(clearFileName, fastq.file))
	} else {
		fileID <- subset(damIdDscrp, grepl(paste(i, "\\..*", sep=""), fastq.file))
		if (nrow(fileID) == 1){
			fileID <- rbind(fileID, fileID)
		}
	}
	if (exists("damIdDscrpCut") == F){
		damIdDscrpCut <- fileID
	} else {
		damIdDscrpCut <- rbind(damIdDscrpCut, fileID)
	}
}
rm(fileID, i)
samplesList <<- as.data.frame(matrix(data=NA, nrow=length(filePath), ncol=6, dimnames=NULL))
names(samplesList) <<- c("id", "tissue", "protein", "conditions", "replicate", "path")
samplesList$path <<- filePath
for (colnumber in c(1:5)){
	for (i in c(1:nrow(samplesList))){
		ins <- sub("([0-9_.a-zA-Z-]+)_(edge|inner)(.*)", "\\2", basename(samplesList[i, 6]), perl=T)
		if (colnumber == 1){
		subst <- paste("\\1\\.\\2\\.\\3_", ins, "\\.\\4", sep="")
		} else if (colnumber %in% c(2,3,5)){
		subst <- paste("\\", colnumber - 1, sep="")
		} else {
		subst <- paste("\\3_", ins, sep="")
		}
		if (length(grep("paired", samplesList[i, 6])) != 0){
			samplesList[i, colnumber] <<- sub("^([a-zA-Z]+)\\.([a-zA-Z0-9]+)\\.([a-zA-Z0-9]+)_(?:R|F)\\.([0-9]+)", subst, damIdDscrpCut[i, 1], perl=TRUE)
		} else {
		samplesList[i, colnumber] <<- sub("^([a-zA-Z]+)\\.([a-zA-Z0-9]+)\\.([a-zA-Z0-9_]+)\\.([0-9]+)", subst, damIdDscrpCut[i, 1], perl=TRUE)
		}
		}
}
rm(i, colnumber, ins, subst, damIdDscrpCut)
write.table(samplesList, file="rdata_description.csv", sep=";", row.names=F, col.names=T, quote=F, dec=".", append=F)
}


# Make samples list file
##########################
MakeSamplesListFile(sourceDir, damIdLocation)

# Load GATC counts in data frame
################################

gatcs <- read.delim(gatcFile, header=T, as.is=T, dec=".")
gatcs <- cbind(gatcs, matrix(data=NA, nrow=nrow(gatcs), ncol=nrow(samplesList)))
	for (i in 1:nrow(samplesList)){
		colnames(gatcs)[startCol+i] <- samplesList$id[i]
		load(file=samplesList$path[i])
		if (all(gatcs$ID.il == reads2GATC$ID)) gatcs[, startCol + i] <- reads2GATC$count
	}

gatcs <- gatcs[!(gatcs$chr %in% c("U", "M", "Uextra" )), ]


for (i in 8:length(gatcs)) {
	if(i %% 2 == 0) {
		y <- i+1
		bmp(filename=file.path(prefixDir, outputPreliminaryStat, paste("Correlations_between_", names(gatcs)[i], "_and_", names(gatcs)[y], ".bmp", sep="")), width=800, height=800, units = "px")
		par(mai=c(1.5, 1.5, 0.5, 0.5))
		par(cex=1.3)
		Cor.P <- round(cor(gatcs[, i], gatcs[, y], method="pearson", use="pairwise.complete.obs"), digits=2)
		Cor.S <- round(cor(gatcs[, i], gatcs[, y], method="spearman", use="pairwise.complete.obs"), digits=2)
		print(ggplot(gatcs, aes(gatcs[, i], gatcs[, y]))+geom_point(alpha=1/10, colour="red", size=4) + xlab(names(gatcs)[i]) + ylab(names(gatcs)[y]) + geom_text(data = data.frame(), size = 4, hjust=0, aes(min(gatcs[, i], na.rm=T), max(gatcs[, y], na.rm=T)*0.75, label =c(paste("Pearson.Cor = ", Cor.P, "\n\n", sep=""), paste("Spearman.Cor = ", Cor.S, sep="")))) + theme_bw())
		dev.off()
	}
}
print("End preliminary correlation count.")


# plot(x=gatcs[, i], y=gatcs[, y], cex=0.3, xlab=names(gatcs)[i], ylab=names(gatcs)[y], las=1, bty="l", pch=".", text(x=min(gatcs[, i], na.rm=T), y=max(gatcs[, y], na.rm=T)*0.75, adj=0, labels=c(paste("pearson = ", Cor.P, "\n\n", sep=""), paste("spearman = ", Cor.S, sep=""))))

# scatter <- ggplot(gatcs, aes(gatcs[, i], gatcs[, y]))

# scatter + geom_point(alpha=1/10, colour="red", size=4) + xlab(names(gatcs)[i]) + ylab(names(gatcs)[y]) + geom_text(data = data.frame(), size = 4, hjust=0, aes(min(gatcs[, i], na.rm=T), max(gatcs[, y], na.rm=T)*0.75, label =c(paste("Pearson.Cor = ", Cor.P, "\n\n", sep=""), paste("Spearman.Cor = ", Cor.S, sep="")))) + theme_bw()