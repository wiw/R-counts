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
	library(ggplot2)
	library(snapCGH)
	library(cluster)
	library(limma)
	library(Vennerable)
	library(BSgenome.Dmelanogaster.UCSC.dm3)
	library(DESeq)
	library(Ringo)
	library(tools)
	library(plyr)
	options(scipen = 999)

	# libraries_string <- c("gplots", "ggplot2", "snapCGH", "cluster", "limma", "Vennerable", "BSgenome.Dmelanogaster.UCSC.dm3", "DESeq", "Ringo", "tools", "plyr")
	# library(libraries_string, character.only=T)

# Declare variables
###################
	prefixDir <- "RUN11-09-2015_ludo" # output directory into working directory
	onlyEdge <- F # use only edge reads to counts or not
	workDir <- "/home/anton/data/R-script/R-counts"	# working directory (WD)
	sourceDir <- "/home/anton/backup/ludo-output" # location your RData files. You can specify the highest folder as it is possible. Searching runs recursively.
	damIdLocation <- "/home/anton/data/DAM/RUN/damid_description.csv" # location your DamID-Description file
	genesFilePath <- "/home/anton/backup/2013-07-09_Assigning_GATCs_to_Genes_(FBgn_IDs)/Drosophila_melanogaster.BDGP5.25.64_Genes_AP130708.txt"
	HKGenesPath <- "/home/anton/backup/2013-07-09_Assigning_GATCs_to_Genes_(FBgn_IDs)/houskeeping_genes.csv"
	exprDataDir <- "/home/anton/backup/input/ExpressionData"
	outputGff <- "gff"	# output folder for gff in WD
	outputWig <- "wig"	# output folder for wig in WD
	outputScttr <- "scatter_plots"	# output folder for scatter plots in WD
	outputDomain <- "domains"
	outputCleanStat <- "clean_stat"
	outputBio <- "Bio"
	outputHeatmap <- "Heatmap"
	outputBioBoxplot <- "Boxplot"
	outputExpr <- "Expression"
	startCol <- 7	# the number of last column in GATCs file, default "7"
	gatcFile <- paste(workDir, "GATCs_mod.txt", sep="/")	# location you GATCs file
	gatcFile4Genes <- paste(workDir, "GATCs.txt", sep="/")
	needCombine <- F	# are you need to combine some columns into one (T or F)? we recommend set the "F"; if you select "T" - edit the "combine" vector on string #85 into this code
	usePseudoCounts <- F	# are you need to add pseudo counts into source data (T or F)? Default "T"
	pseudoCounts <- c(0.01)		# the vector of pseudo counts, default "c(1)"
	corrMethod <- c("pearson", "spearman")	# the vector of type correlation method - please don't change this setting
	heatmapColors <- greenred(200)	# check you color into heatmap
	labelHeatmap <- c("B", "C")	# the vector of heatmap label - please don't change this setting
	labelAcf <- c("A_ALL", "B_MA")	# the vector of acf label - please don't change this setting
	writeTemp <- T		# use this option if you need to get the intermediate files, default "T"
	needSomeFiles <- F	# if need calculate only some files from source list, default "F"
	useSomeFiles <- c(9:16)	# the region of files which need calculate
	tissue.bio.set <- c("BR", "FB", "Glia", "NRN", "Kc167") # Part of tissies from all
	protein.bio.set <- c("LAM", "HP1", "PC") # Part of protein from all
	conditions.bio.set <- c("m", "m_25mkM4HT", "mf_min", "std")
	tissueExprSet <- c("BR", "FB", "Glia", "NRN", "Kc167")
	tissueIntersectSet <- c("BR", "FB", "Glia", "NRN", "Kc167")
	proteinIntersectSet <- c("LAM", "HP1", "PC")


# Create folders
################
	dir.create(file.path(workDir, prefixDir), showWarnings = FALSE)
	dir.create(file.path(workDir, prefixDir, outputGff), showWarnings = FALSE)
	dir.create(file.path(workDir, prefixDir, outputWig), showWarnings = FALSE)
	dir.create(file.path(workDir, prefixDir, outputScttr), showWarnings = FALSE)
	dir.create(file.path(workDir, prefixDir, outputGff, outputDomain), showWarnings = FALSE)
	dir.create(file.path(workDir, prefixDir, outputCleanStat), showWarnings = FALSE)
	dir.create(file.path(workDir, prefixDir, outputHeatmap), showWarnings = FALSE)
	dir.create(file.path(workDir, prefixDir, outputBio), showWarnings = FALSE)
	dir.create(file.path(workDir, prefixDir, outputBio, outputBioBoxplot), showWarnings = FALSE)


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

# Make samples list file
##########################
MakeSamplesListFile <- function(SOURCE, DAMID) {
	filePath <- list.files(path=SOURCE, pattern="*_local_GATCcounts.RData", full.names=T, recursive=T)
	if (all(grepl("(edge|inner)", basename(filePath), perl=T))) {
	baseFile <- unique(sub("(.*)_(edge|inner).*", "\\1", basename(filePath), perl=T))
	ludoLabel <<- F
	} else {
		baseFile <- unique(sub("(.*)_local_GATCcounts.RData", "\\1", basename(filePath), perl=T))
		ludoLabel <<- T
	}
	damIdDscrp <- read.delim(DAMID, header=T, sep="\t", stringsAsFactors=F)
	for (i in baseFile) {
		if (length(grep("paired", i)) != 0){
		clearFileName <- sub("(.*)_paired", "\\1", i, perl=T)
		fileID <- subset(damIdDscrp, grepl(clearFileName, fastq.file))
		} else {
			fileID <- subset(damIdDscrp, grepl(paste(i, "\\..*", sep=""), fastq.file))
			if (nrow(fileID) == 1 & ludoLabel == F){
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
			if (ludoLabel == F) {
				ins <- sub("([0-9_.a-zA-Z-]+)(_edge|_inner)(.*)", "\\2", basename(samplesList[i, 6]), perl=T)
			} else {
				ins <- character()
			}
			if (colnumber == 1){
			subst <- paste("\\1\\.\\2\\.\\3", ins, "\\.\\4", sep="")
			} else if (colnumber %in% c(2,3,5)){
			subst <- paste("\\", colnumber - 1, sep="")
			} else {
			subst <- paste("\\3", ins, sep="")
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

# Write intermediate files function
################################### 
	WriteIntermediateFiles <- function(source, output.file) {
		if (writeTemp == T) {
			if (tolower(names(source)[1]) == "id") {
				names(source)[1] <- paste("GATC", names(source)[1], sep=".")
			}
			write.table(source, file=file.path(prefixDir, output.file), sep=";", row.names=F, col.names=T, quote=F, dec=".", append=F, eol="\r\n")
		}
	}

# Pearson and spearman correlations function
############################################
PearsonAndSpearmanCorrelations <- function(dataSet, use.method, use.opt) {
	corr <- matrix(data=NA, nrow=ncol(dataSet)-7, ncol=ncol(dataSet)-7, byrow=T)
		rownames(corr) <- names(dataSet)[8:(ncol(dataSet))]
		colnames(corr) <- names(dataSet)[8:(ncol(dataSet))]
		for (j in 1:(ncol(dataSet)-7)){
			for (i in 1:(ncol(dataSet)-7)){
				corr[i,j] <- round(cor(dataSet[,7+i], dataSet[,7+j], method=use.method, use=use.opt), digits=2)
			}
			rm(i)
		}
	rm(j)
		invisible(corr)
}
PearsonAndSpearmanCorrelationsHeatmapMod <- function(dataSet1, dataSet2, use.method, use.opt) {
	corr <- matrix(data=NA, nrow=ncol(dataSet1)-7, ncol=ncol(dataSet2)-7, byrow=T)
		rownames(corr) <- names(dataSet1)[8:(ncol(dataSet1))]
		colnames(corr) <- names(dataSet2)[8:(ncol(dataSet2))]
		for (j in 1:(ncol(dataSet2)-7)){
			for (i in 1:(ncol(dataSet1)-7)){
				corr[i,j] <- round(cor(dataSet1[,7+i], dataSet2[,7+j], method=use.method, use=use.opt), digits=2)
			}
			rm(i)
		}
	rm(j)
		invisible(corr)
}

# Combine samples function
##########################
CombineSamples <- function(dataFrame) {
	with(dataFrame, data.frame(DAM.1=`DAM-1.FCC4JPEACXX_L6_R1` + `DAM-1.FCC4JPEACXX_L6_R2`, DAM.2=`DAM-2.FCC4JPEACXX_L6_R1` + `DAM-2.FCC4JPEACXX_L6_R2`, LAM.1=`LAM-1.FCC4JPEACXX_L6_R1` + `LAM-1.FCC4JPEACXX_L6_R2`, LAM.2=`LAM-2.FCC4JPEACXX_L6_R1` + `LAM-2.FCC4JPEACXX_L6_R2`))
}

# Main correlations function
############################
MainCorrelations <- function(dataSet, corrMethod, labelHeatmap, use.opt="everything", suffixCSV, suffixPDF, corr.on.file, createPDF=T, counts=F) {  
	for (m in corrMethod) {
		corr.file <- assign(paste(m, ".cor", sep=""), as.data.frame(PearsonAndSpearmanCorrelations(dataSet, m, use.opt)))
			if (writeTemp == T) {
				write.table(corr.file, file=file.path(prefixDir, paste("DF_Counts_Step_", suffixCSV, "_", if (counts == T) {""} else {paste(name, "_", sep="")}, m, "_", currentDate, ".csv", sep="")), sep=";", row.names=T, col.names=T, quote=F, dec=".", append=F, eol="\r\n")
			}
		corr.file <- as.matrix(corr.file)
			if (createPDF == T){
				for (label in labelHeatmap) {
					options(warn=-1)
						pdf(file=file.path(prefixDir, paste("DF_Counts_Step_", suffixPDF, "_", label, "_", name, "_", m, "_", currentDate, ".pdf", sep="")), width=12, height=12)
						if (label == "B") {
							KEY <- F
								densityInfo <- "density"
						} else {
							KEY <- T
								densityInfo <- "none"
						}
					heatmap.2(x=corr.file, col=heatmapColors, breaks=seq(from=-1, to=1, by=0.01), Rowv=T, Colv=T, dendrogram="both", trace="none", cellnote=corr.file, notecol="white", notecex = 0.5, margins=c(7,7), main=paste(m, "'s correlation coefficients and hierarchical clustering on\n", "'", corr.on.file, "'", sep=""), cex.lab=1.1, cexRow=0.6, cexCol=0.6, lmat=matrix(c(4,2,3,1), ncol=2), lwid=c(0.1, 0.9), lhei=c(0.15, 0.85), key=KEY, density.info=densityInfo)
						options(warn=0)
						dev.off()
				}  
			}
	}
}

# ACF on data function
######################
AcfOnData <- function(dataSet, labelAcf, method, suffixPDF, ylab.val, na.data) {
	for (label in labelAcf) {
		pdf(file=file.path(prefixDir, paste("DF_Counts_Step_", suffixPDF, "_", label, "_", name, "_", currentDate, ".pdf", sep="")), width=11.69, height=8.27)
			par(mfrow=c(3, 4))
			par(mai=c(0.7, 0.7, 0.7, 0.5))
			if (label == "A_ALL") {
				DATAs.acf <- dataSet
					descr <- "all"
			} else {
				DATAs.acf <- dataSet[dataSet$presence.ma == 1, ]
					descr <- "ma"
			}
		acf.order.list <- order(names(DATAs.acf)[8:ncol(DATAs.acf)]) + 7
			if (method == "acf") {    
				for (i in acf.order.list) {
					if (na.data == T) {
						acf.na.data <- sum(!is.na(DATAs.acf[, i]))
					} else {
						acf.na.data <- sum(DATAs.acf[, i] != 0)
					}
					acf(DATAs.acf[, i], na.action=na.pass, main=paste(names(DATAs.acf[i]), "\n(", descr," GATCs: ", acf.na.data, "\nout of ", nrow(DATAs.acf), ")", sep=""), ylab=ylab.val)
				}
			} else {
				for (i in acf.order.list) {
					if (na.data == T) {
						acf.na.data <- sum(!is.na(DATAs.acf[, i]))
					} else {
						acf.na.data <- sum(DATAs.acf[, i] != 0)
					}
					plot(density(DATAs.acf[, i], na.rm=T), main=paste(names(DATAs.acf[i]), "\n(", descr," GATCs: ", acf.na.data, "\nout of ", nrow(DATAs.acf), ")", sep=""), ylab=ylab.val)
				}
			}		
		rm(i)
			dev.off()
	}
}

# DamID to WIG&GFF function
###########################
DamIdSeqToWigGff <- function(dataSet) {
	for (step in c(1:2)) {
		if (step == 1) {
			tag <- ""
			data.wg <- dataSet
		} else {
			tag <- ".ma"
			data.wg <- dataSet[dataSet$presence.ma == 1, ]
		}

	# Calculate WIG
		data.wg$start <- round((data.wg$start + data.wg$end)/2)
		data.wg$start <- sprintf("%d", data.wg$start)
		for (j in 8:(ncol(data.wg))) {
			wig.file <- file.path(prefixDir, outputWig, paste(names(data.wg)[j], tag, "_", name, ".wig", sep=""))
			chrs <- unique(data.wg$chr)
			for (i in 1:length(unique(data.wg$chr))){
				selected.chr <- data.wg[(data.wg$chr == chrs[i]), c(3, j)]
				selected.chr <- selected.chr[!is.na(selected.chr[, 2]), ]
				if (i == 1) {
					write.table(paste("variableStep chrom=", chrs[i], sep=""), file=wig.file, sep="\t", row.names=F, col.names=F, quote=F, dec=".", append=F)
				} else {
					write.table(paste("variableStep chrom=", chrs[i], sep=""), file=wig.file, sep="\t", row.names=F, col.names=F, quote=F, dec=".", append=T)
				}
				write.table(selected.chr, file=wig.file, sep=" ", row.names=F, col.names=F, quote=F, dec=".", append=T)
			}
		}
		rm(j, i)	

	# Calculate GFF
			for (j in 8:(ncol(data.wg))) {
				selected.set <- data.wg[, c(1:7, j)];
				selected.set <- cbind(selected.set, NA);
				selected.set[, 2] <- paste("chr", selected.set[, 2], sep="");
				selected.set <- selected.set[, c(2, 7, 5, 3, 4, 8, 7, 9, 1)];
				selected.set[, 2] <- ".";
				selected.set[, 3] <- paste(names(data.wg)[j], tag, sep="");
				selected.set[, 7] <- ".";
				selected.set[, 8] <- ".";
				selected.set <- selected.set[!is.na(selected.set[, 6]), ];
				gff.file <- file.path(prefixDir, outputGff, paste(names(data.wg)[j], tag, "_", name, ".gff", sep=""))
					write.table(selected.set, file=gff.file, sep="\t", row.names=F, col.names=F, quote=F, dec=".", append=F);
			}
			rm(j)
		}
}

# Scatter Plots on Averaged data function
#########################################
# ScatterPlotting <- function(dataSet, tag) {
# 	for (j in names(dataSet)[-c(1:7)]){
# 			i <- grep(sub("(.+)\\.[1-2]\\.norm", "\\1", j, perl=T), names(dataSet)[-c(1:7)], value=T)
# 			if (length(i) != 1){
# 					bmp(filename=file.path(prefixDir, outputScttr, paste("scatter_", sub("(.+)\\.norm", "\\1", i, perl=T)[1], "_vs_", sub("(.+)\\.norm", "\\1", i, perl=T)[2], "_", name, "_", tag, ".bmp", sep="")), width=600, height=600, units = "px")
# 					Cor.P <- round(cor(dataSet[[i[1]]], dataSet[[i[2]]], method="pearson", use="pairwise.complete.obs"), digits=2)
# 					Cor.S <- round(cor(dataSet[[i[1]]], dataSet[[i[2]]], method="spearman", use="pairwise.complete.obs"), digits=2)
# 					plot(x=dataSet[[i[1]]], y=dataSet[[i[2]]], cex=0.3, xlab=i[1], ylab=i[2], text(x=min(dataSet[[i[1]]], na.rm=T) + 0.5, y=max(dataSet[[i[2]]], na.rm=T) - 0.5, labels=c(paste("r = ", Cor.P, "\n\n", sep=""), paste("s = ", Cor.S, sep=""))))
# 					rm(Cor.P)
# 					rm(Cor.S)
# 					dev.off()
# 			} else {
# 				print(paste("Skip make the Scatter Plots from", j, sep=" "))
# 			}
# 	}
# }
# Scatter Plots on Averaged data function
#########################################
ScatterPlotting3D <- function(dataSet, tag) {
	for (j in unique(sub("([a-zA-Z_\\.]+)\\.[0-9].*", "\\1", names(dataSet)[-c(1:7)], perl=T))) {
		repSet <- sort(grep(j, names(dataSet), value=T))
		if (length(repSet) > 1) {
		bmp(filename=file.path(prefixDir, outputScttr, paste("scatter_on_", j, "_", tag, ".bmp", sep="")), width=600*length(repSet), height=600*(length(repSet)-1), units = "px")
		par(mfrow=c(length(repSet)-1, length(repSet)))
		par(mai=c(1.5, 1.5, 0.7, 0.5))
		par(cex=1.5)
			for (i in repSet) {
				for (x in repSet) {
					if (x != i) {
						Cor.P <- round(cor(dataSet[[i]], dataSet[[x]], method="pearson", use="pairwise.complete.obs"), digits=2)
						Cor.S <- round(cor(dataSet[[i]], dataSet[[x]], method="spearman", use="pairwise.complete.obs"), digits=2)
						plot(x=dataSet[[i]], y=dataSet[[x]], cex=0.3, xlab=i, ylab=x, text(x=min(dataSet[[i]], na.rm=T) + 1.5, y=max(dataSet[[i]], na.rm=T)*0.75, labels=c(paste("r = ", Cor.P, "\n\n", sep=""), paste("s = ", Cor.S, sep=""))))
						# print(ggplot(dataSet, aes(dataSet[[i]], dataSet[[x]]))+geom_point(alpha=1/10, colour="red", size=4) + xlab(i) + ylab(x) + geom_text(data = data.frame(), size = 4, hjust=0, aes(min(dataSet[, i], na.rm=T), max(dataSet[, x], na.rm=T)*0.75, label =c(paste("Pearson.Cor = ", Cor.P, "\n\n", sep=""), paste("Spearman.Cor = ", Cor.S, sep="")))) + theme_bw())
						rm(Cor.P)
						rm(Cor.S)
					}
				}
			}
		dev.off()
		} else {
			print(paste("Skip make the Scatter Plots from", j, sep=" "))
		}
	}
}
# FeatureCalls to GFF function
##############################
FeatureCalls.to.GFF.like <- function(start.coordinate, end.coordinate, feature.type) {
   e.ind <- cumsum(rle(feature.type)$lengths)
   s.ind <- c(1, e.ind[1:(length(e.ind)-1)]+1)
   gff <- data.frame(start=start.coordinate[s.ind], end=end.coordinate[e.ind], value=feature.type[s.ind])
   invisible(gff)
}

# Run biological Headen Mark Model function
###########################################
runBioHMM <- function (mval, datainfo, useCloneDists = TRUE, criteria = "AIC", delta = NA, var.fixed = FALSE, epsilon = 1e-06, numiter = 30000) {
	crit = TRUE
	if (criteria == "AIC") {
		aic = TRUE
		bic = FALSE
	} else if (criteria == "BIC") {
		bic = TRUE
		aic = FALSE
	} else crit = FALSE
	if ((crit == 1) || (crit == 2)) {
		if (criteria == "BIC") {
			if (is.na(delta)) {
				delta <- c(1)
			}
		}
		res <- try(fit.model(obs = mval, datainfo = datainfo, useCloneDists = useCloneDists, aic = aic, bic = bic, delta = delta, var.fixed = var.fixed, epsilon = epsilon, numiter = numiter))$out.list$state
	}
	else {
		cat("You must enter AIC or BIC for the criteria argument\n")
	}
}

# Use fit third-clustering model function 
#########################################
fit.model <- function (obs, datainfo = NULL, useCloneDists = TRUE, aic = TRUE, bic = FALSE, delta = 1, var.fixed = FALSE, epsilon = 1e-06, numiter = 30000) {
    kb <- datainfo$start 
    if (useCloneDists) {
        dists.pre = kb[2:length(kb)] - kb[1:(length(kb) - 1)] 
        dists = dists.pre/(max(dists.pre)) 
    } else {
        dists <- rep(1, length(kb))
    }
    covars <- as.matrix(dists) 
    obs.ord <- obs[order(kb)] 
    kb.ord <- kb[order(kb)] 
    ind.nonna <- which(!is.na(obs.ord)) 
    data <- obs.ord[ind.nonna] 
    kb <- kb.ord[ind.nonna] 
    numobs <- length(data) 
    if (numobs > 5) { 
        temp3 <- clara(data, 3)
        init.mean.three <- vector(length = 3)
        init.mean.three <- temp3$medoids
        init.var.three <- vector(length = 3)
        if (var.fixed == FALSE) {
            for (i in 1:3) {
                if (length(temp3$data[temp3$clustering == i]) > 1) {
               		init.var.three[i] <- log(sqrt(var(temp3$data[temp3$clustering == i])))
            	} else {
            		init.var.three[i] <- log(0.5)
            	}
            }
        } else {
            init.var.three[1:3] <- log(sqrt(var(data)))
        }
        z3.init <- c(init.mean.three[, 1], init.var.three, -0.7, 
            -0.7, -3.6, -3.6, -3.6, -3.6, -3.6, -3.6, 0)
        z.pre <- run.nelder(numobs, z3.init, data, covars, var.fixed, epsilon, numiter, i)
        if (!is.nan(z.pre$x[1])) {
            z3 <- find.param.three(z.pre, var.fixed)
        } else {
            z3 <- NULL
        }
		z <- z3
		trans.mat <- list()
        for (j in 1:(length(data) - 1)) {
            trans.mat[[j]] <- z$LH.trans + exp(-(covars[j, 1]^(z$rate1)) * prod(covars[j, -1])) * z$RH.trans
        }
        Vit.seg <- Viterbi.three(data, z, trans.mat)
        maxstate.unique <- unique(Vit.seg)
        mean <- rep(0, length(data))
        var <- rep(0, length(data))
        for (m in 1:length(maxstate.unique)) {
            mean[Vit.seg == maxstate.unique[m]] <- mean(data[Vit.seg == 
              maxstate.unique[m]])
            var[Vit.seg == maxstate.unique[m]] <- var(data[Vit.seg == maxstate.unique[m]])
        }
        out <- cbind(matrix(Vit.seg, ncol = 1), matrix(mean, ncol = 1), matrix(var, ncol = 1))
        out.all <- matrix(NA, nrow = length(kb.ord), ncol = 4)
        out.all[ind.nonna, 1:3] <- out
        out.all[, 4] <- obs.ord
        out.all <- as.data.frame(out.all)
        dimnames(out.all)[[2]] <- c("state", "mean", "var", "obs")
        numstates <- length(unique(Vit.seg))
    } else {
        out.all <- matrix(NA, nrow = length(kb.ord), ncol = 4)
        out.all[ind.nonna, 1] <- c(rep(1, numobs))
        out.all[ind.nonna, 2] <- c(rep(mean(obs.ord), numobs))
        out.all[ind.nonna, 3] <- c(rep(var(obs.ord), numobs))
        out.all[ind.nonna, 4] <- obs.ord
        out.all <- as.data.frame(out.all)
        dimnames(out.all)[[2]] <- c("state", "mean", "var", "obs")
        numstates = 1
    }
    list(out.list = out.all, nstates.list = numstates)
}

# Calculate coordinate for filter
#################################
CalculateCoordinate <- function(x, y){
  if (index[x] == T){
    RhsMatrix <- matrix(c(-1*Intercept,-(DATA.filter[x, 9]+DATA.filter[x, 8]/Slope)), nrow=2, ncol=1, byrow=T)
    Result <- Coef.Matrix %*% RhsMatrix
    Result[y,1]
  }
}

# Make Venn diagram
###################
pcentFun <- function(x,y) {
	100 * (x / y)
}
# !!! Before use this function, please perform this instructions: http://stackoverflow.com/a/15315369/4424721
#############################################################################################################
MakeVennDiagram <- function(x, set, v1, v2) {
	compare.by.venn.df <- x[, grep(set, names(x)[8:ncol(x)], perl=T, value=T)]
	compare.by.venn.df <- compare.by.venn.df[complete.cases(compare.by.venn.df), ]
	category <- names(compare.by.venn.df)
	venn.matrix <- vennCounts(compare.by.venn.df)
	areas <- venn.matrix[,ncol(venn.matrix)]
	names(areas) <- apply(venn.matrix[, c(1:(ncol(venn.matrix)-1))], 1, paste, collapse="")
	areasPcent <- round(pcentFun(areas, sum(areas)), digits=2)
	VennSet <- Venn(SetNames = category, Weight = areas)
	VennList <- compute.Venn(VennSet, doWeights = TRUE)
	VennList_pc <- VennList
	for (i in c(1:length(areasPcent))) {
		position <- grep(names(areasPcent)[i], VennList_pc@FaceLabels$Signature)
		VennList_pc@FaceLabels$Signature[position] <- paste(areasPcent[i], "%", sep="")
	}
	gp <- VennThemes(VennList)
	for (i in names(gp$FaceText)) gp$FaceText[[i]]$fontsize <- 12
	for (i in names(gp$SetText)) gp$SetText[[i]]$fontsize <- 15
	pdf(file=file.path(prefixDir, outputBio, paste("Venn Diagram between ", v2, " from ", v1, ".pdf", sep="")), width=8, height=12.76)
	par(mfrow=c(1, 2))
	par(mai=c(0.7, 0.7, 0.7, 0.5))
	if (ncol(venn.matrix) > 4) {
		rowNumberMakevp <<- 1
		plot(VennList)
		rowNumberMakevp <<- 2
		plot(VennList_pc, show = list(FaceText = c("signature"), DarkMatter = TRUE))
		} else {
			rowNumberMakevp <<- 1
			plot(VennList, gpList=gp)
			rowNumberMakevp <<- 2
			plot(VennList_pc, gpList=gp, show = list(FaceText = c("signature"), DarkMatter = TRUE))
		}
	dev.off()
}

# Generate heatmap from custom selection
########################################
MakeHeatmapToPdf <- function(input, item1, item2, vs_name, method) {
	input <- as.matrix(input)
	options(warn=-1)
	pdf(file=file.path(prefixDir, outputHeatmap, paste("Heatmap_", item1, "_vs_", item2, "_on_", paste(vs_name, collapse="_"), "_", method, "_", currentDate, ".pdf", sep="")), width=14, height=14)
	heatmap.2(x=input, col=heatmapColors, breaks=seq(from=-1, to=1, by=0.01), Rowv=T, Colv=T, dendrogram="both", trace="none", cellnote=input, notecol="white", notecex = 1.5, margins=c(7,7), main=paste(method, "'s correlation coefficients and hierarchical clustering on\n", "'", paste(vs_name, collapse=", "), "'", sep=""), cex.lab=1.1, cexRow=0.8, srtRow=90, cexCol=0.8, srtCol=0, lmat=matrix(c(4,2,3,1), ncol=2), lwid=c(0.1, 0.9), lhei=c(0.15, 0.85), key=T, density.info="none")
	options(warn=0)
	dev.off()
}
HeatmapBySelection <- function(dataSet, corrMethod, use.opt) {  
	for (m in corrMethod) {
		log <- list()
		bio.set <- list(TB = tissue.bio.set, PB = protein.bio.set)
		for (element.set in bio.set) {
			for (tHa in element.set) {
				for (tHb in element.set) {
					if (length(element.set) == length(bio.set$TB)) {
						selection.filename <- bio.set$PB
					} else {
						selection.filename <- bio.set$TB
					}
					if (tHa == tHb) {
						log[[tHa]] <- tHa
						ds.selection <- cbind(dataSet[, c(1:7)], dataSet[,grep(paste(tHa, "\\..*",sep=""), names(dataSet), perl=T, value=T)])
						corr.file <- assign(paste(m, ".cor", sep=""), as.data.frame(PearsonAndSpearmanCorrelations(ds.selection, m, use.opt)))
						MakeHeatmapToPdf(input=corr.file, item1=tHa, item2=tHb, vs_name=selection.filename, method=m)
					} else {
						tgh <- assign(name, paste(tHa, tHb, sep="_"))
						log[[tgh]] <- paste(tHa, tHb, sep="_")
						if (length(grep(paste("(",tHa,"_",tHb,")|(",tHb,"_",tHa,")",sep=""), log, perl=T)) < 2) {
							ds.selection1 <- cbind(dataSet[, c(1:7)], dataSet[,grep(paste(tHa, "\\..*",sep=""), names(dataSet), perl=T, value=T)])
							ds.selection2 <- cbind(dataSet[, c(1:7)], dataSet[,grep(paste(tHb, "\\..*",sep=""), names(dataSet), perl=T, value=T)])
							corr.file <- assign(paste(m, ".cor", sep=""), as.data.frame(PearsonAndSpearmanCorrelationsHeatmapMod(ds.selection1, ds.selection2, m, use.opt)))
							MakeHeatmapToPdf(input=corr.file, item1=tHa, item2=tHb, vs_name=selection.filename, method=m)
						} else {
							log[[tgh]] <- NULL
						}
					}
				}
			}
		}
	}
}
# Function for generate multiplot grafs
#######################################
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
  }
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
    }
  }
}

#######################################
################# END #################
#######################################
print("Run script")

# Make samples list file
##########################
MakeSamplesListFile(sourceDir, damIdLocation)

# Load GATC counts in data frame
################################
if (startCol == 0) {
	step01 <- read.delim(alreadyRun, header=T, as.is=T, dec=".")
	startCol <- ncol(step01)
	gatcs <- step01
} else {
	gatcs <- read.delim(gatcFile, header=T, as.is=T, dec=".")
}
if (ludoLabel == F) {
	if (onlyEdge == T) {
		samplesList <- samplesList[grep("edge" ,samplesList$id), ]
	} else {
		modS <- samplesList[grep("edge", samplesList$id), ]
		modS$id <- sub("(.+)(edge)(.+)", paste("\\1", "all", "\\3", sep=""), modS$id, perl=T)
		modS$conditions <- sub("(.+)edge", paste("\\1", "all", sep=""), modS$conditions, perl=T)
	}
} else {
	modS <- samplesList
	modS$id <- sub("(.+)(\\.[0-9]?)$", paste("\\1", "_all", "\\2", sep=""), modS$id, perl=T)
	modS$conditions <- sub("(.+)", paste("\\1", "_all", sep=""), modS$conditions, perl=T)
}
if (needSomeFiles == T) {
	samplesList <- samplesList[useSomeFiles, ]
}
gatcs <- cbind(gatcs, matrix(data=NA, nrow=nrow(gatcs), ncol=nrow(samplesList)))
for (i in 1:nrow(samplesList)){
	colnames(gatcs)[startCol+i] <- samplesList$id[i]
	load(file=samplesList$path[i])
	if (all(gatcs$ID.il == reads2GATC$ID)) gatcs[, startCol + i] <- reads2GATC$count
}
rm(i)
if (ludoLabel == F) {
	if (onlyEdge != T) {
		modG <- gatcs[, c(1:7, grep("edge", names(gatcs)))]
		names(modG)[8:ncol(modG)] <- gsub("(.+)(edge)(.+)", paste("\\1", "all", "\\3", sep=""), names(modG)[8:ncol(modG)], perl=T)
		for (enzyme in names(modG[8:ncol(modG)])) {
			S <- gsub("(.+)(all)(.+)", paste("\\1", "edge", "\\3", sep=""), enzyme, perl=T)
			E <- gsub("(.+)(all)(.+)", paste("\\1", "inner", "\\3", sep=""), enzyme, perl=T)
			modG[[enzyme]] <- gatcs[[S]] + gatcs[[E]]
		}
	gatcs <- modG
	samplesList <- modS
	# rm(modS, modG)
	}
} else {
	names(gatcs)[8:ncol(gatcs)] <- sub("(.+)(\\.[0-9]?)$", "\\1_all\\2", names(gatcs)[8:ncol(gatcs)], perl=T)
	samplesList <- modS
}
currentDate <- format(Sys.time(), "%d-%m-%Y")
load.gatc.df <- paste("DF_Counts_Step_01_Raw_Counts_", currentDate, ".csv", sep="")
WriteIntermediateFiles(source=gatcs, output.file=load.gatc.df)


# Remove unimportant data
#########################
	print("Remove unimportant data")
	DATA <- gatcs[!(gatcs$chr %in% c("U", "M", "Uextra" )), ]
	use.chr.only <- paste("DF_Counts_Step_02_Useful_Chrs_Only_", currentDate, ".csv", sep="")
WriteIntermediateFiles(source=DATA, output.file=use.chr.only)

# Filter DATA
#############
print("Start filtering data")
set <- unique(sub("(.*)\\.([0-9]?)$", "\\1", names(DATA)[8:length(DATA)], perl=T))
DATA.outfilter <- DATA[, 1:7]
for (i in set){
	protein.set <- grep(i, names(DATA))
	print(paste("Filter data from", i, sep=" "))
	if (length(protein.set) == 1) {
		print(paste("You have one replicate from ", i, ". Skipped!", sep=""))
		DATA.outfilter[[grep(i, names(DATA), value=T)]] <- DATA[, protein.set]
	}
	if (length(protein.set) >= 2) {
		if (length(protein.set) > 2) { 
			cor.df <- as.data.frame(matrix(0, nrow=length(protein.set), ncol=length(protein.set)))
			colnames(cor.df) <- protein.set
			rownames(cor.df) <- protein.set
			for (cor1 in protein.set) {
				for (cor2 in protein.set) {
					if (cor2 != cor1) {
					cor.df[grep(cor1, rownames(cor.df)), grep(cor2, colnames(cor.df))]	<- cor(DATA[, cor1], DATA[, cor2], method="pearson", use="pairwise.complete.obs")
					}
				}
			}	
			sample2 <- as.numeric(colnames(cor.df)[max(grep(max(cor.df), cor.df))])
			sample1 <- as.numeric(rownames(cor.df)[grep(max(cor.df), cor.df[, max(grep(max(cor.df), cor.df))])])
			DATA.filter <- DATA[,c(1:7, sample1, sample2)]
		} else {
			DATA.filter <- DATA[,c(1:7, grep(i, names(DATA)))]
		}
		for (name in c(8:9)) DATA.filter[, name][DATA[, name]==0] <- NA
		index <- complete.cases(DATA.filter[,8], DATA.filter[,9])
		Rep.Slope <- lm(DATA.filter[,9] ~ DATA.filter[,8])
		Intercept <- Rep.Slope$coefficients[[1]]
		Slope <- Rep.Slope$coefficients[[2]]
		Coef.Matrix <- solve(matrix(c(Slope,-1,-1/Slope,-1), nrow=2, ncol=2, byrow=T))
		DATA.nrow <- c(1:length(index))
	
		X.Coord <- sapply(DATA.nrow, CalculateCoordinate, y=1)
		X.Coord[sapply(X.Coord, is.null)] <- NA
		DATA.filter$Rep.Slope.X <- unlist(X.Coord)
	
		Y.Coord <- sapply(DATA.nrow, CalculateCoordinate, y=2)
		Y.Coord[sapply(Y.Coord, is.null)] <- NA
		DATA.filter$Rep.Slope.Y <- unlist(Y.Coord)
	
		index.sorted <- order(DATA.filter$Rep.Slope.X, decreasing=F, na.last=NA)
		results.df <- as.data.frame(matrix(nrow=0, ncol=3))
	
		DATA.list.out <- list()
		bmp(filename=file.path(prefixDir, outputCleanStat, paste("Correlations_between_", names(DATA.filter)[8], "_and_", names(DATA.filter)[9], ".bmp", sep="")), width=1600, height=800, units = "px")
		par(mfrow=c(1, 2))
		par(mai=c(1.5, 1.5, 0.5, 0.5))
		par(cex=1.3)
		for (bin.size in c(0, 100, 200, 500, seq(1000, 10000, 1000), 15000, seq(20000, 50000, 10000), sum(!is.na(DATA.filter$Rep.Slope.X)))) {
			if (bin.size == 0) {
				Pearson.Cor <- cor(DATA.filter[, 8], DATA.filter[, 9], method="pearson", use="pairwise.complete.obs")
				Spearman.Cor <- cor(DATA.filter[, 8], DATA.filter[, 9], method="spearman", use="pairwise.complete.obs")
				RemovedGATCs <- nrow(DATA[DATA[,protein.set[1]] == 0 | DATA[,protein.set[2]] == 0, ])
				LeftGATCs <- nrow(DATA) - RemovedGATCs
				results.df <- rbind(results.df, data.frame("Bin.Size"=bin.size, "Number.of.Removed.GATCs"=RemovedGATCs, "Pearson.Cor"=Pearson.Cor, "Number.of.left.GATCs"=LeftGATCs))
				plot(x=DATA.filter[, 8], y=DATA.filter[, 9], cex=0.3, xlab=names(DATA.filter)[8], ylab=names(DATA.filter)[9], las=1, bty="l", pch=".", text(x=min(DATA.filter[, 8], na.rm=T), y=max(DATA.filter[, 9], na.rm=T)*0.75, adj=0, main="Original correlation" , labels=c(paste("pearson = ", round(Pearson.Cor, digits=2), "\n\n", sep=""), paste("spearman = ", round(Spearman.Cor, digits=2), sep=""))))

			} else {
				number.of.bins <- ceiling(sum(!is.na(DATA.filter$Rep.Slope.X))/bin.size)
				for (item in 1:number.of.bins) {
					if (item < number.of.bins) DATA.ordered.subset <- DATA.filter[(index.sorted[(1+bin.size*(item-1)):(bin.size*item)]),]
					if (item == number.of.bins) DATA.ordered.subset <- DATA.filter[(index.sorted[(1+bin.size*(item-1)):(sum(!is.na(DATA.filter$Rep.Slope.X)))]),]
					DATA.ordered.subset$Rep.Ratio <- log2(DATA.ordered.subset[, 8]/DATA.ordered.subset[, 9])
					filter.index <- (DATA.ordered.subset$Rep.Ratio > (mean(DATA.ordered.subset$Rep.Ratio) - 2*sd(DATA.ordered.subset$Rep.Ratio))) & (DATA.ordered.subset$Rep.Ratio < (mean(DATA.ordered.subset$Rep.Ratio) + 2*sd(DATA.ordered.subset$Rep.Ratio)))
					for (name in c(8:9)) DATA.ordered.subset[[names(DATA.ordered.subset)[name]]][filter.index == F] <- NA
					if (item == 1) {
						DATA.filtered <- DATA.ordered.subset
						} else {
							DATA.filtered <- rbind(DATA.filtered, DATA.ordered.subset)
						}
				}
				DATA.merged <- DATA[, c(1:7, grep(names(DATA.filter)[8], names(DATA)), grep(names(DATA.filter)[9], names(DATA)))]
				for (name in paste(names(DATA.filter)[8:9], "filt", sep=".")) DATA.merged[,name] <- NA
				index.merge <- match(DATA.merged$ID[index.sorted], DATA.filtered$ID)
				for (name in 10:11) DATA.merged[index.sorted, name] <- DATA.filtered[index.merge, name-2]
				DATA.merged$note <- "kept"
				DATA.merged$note[(!is.na(DATA.merged[, 8] + DATA.merged[, 9]) + !is.na(DATA.merged[, 10] + DATA.merged[, 11])) == 1] <- "filtered"
				RemovedGATCs <- sum(DATA.merged$note == "filtered")
				Pearson.Cor <- cor(DATA.merged[, 10], DATA.merged[, 11], method="pearson", use="pairwise.complete.obs")
				LeftGATCs <- nrow(DATA) - RemovedGATCs
				results.df <- rbind(results.df, data.frame("Bin.Size"=bin.size, "Number.of.Removed.GATCs"=RemovedGATCs, "Pearson.Cor"=Pearson.Cor, "Number.of.left.GATCs"=LeftGATCs))
				for (name in 10:11) DATA.merged[is.na(DATA.merged[, name]), name] <- 0
				DATA.list.out[[paste("bins", bin.size, sep="_")]] <- DATA.merged[, 10:11]
				names(DATA.list.out[[paste("bins", bin.size, sep="_")]]) <- sub("(.*)\\.filt", "\\1", names(DATA.list.out[[paste("bins_", bin.size, sep="")]]), perl=T)
				write.results.df <- file.path(outputCleanStat, paste("Clean_results_for_", i, ".csv", sep=""))
				WriteIntermediateFiles(source=results.df, output.file=write.results.df)
			}
		}
		results.df.part <- results.df[results.df$Pearson.Cor > (max(results.df$Pearson.Cor) - 0.1*sd(results.df$Pearson.Cor)),]
		Use.Bin.Size <- results.df.part[results.df.part$Number.of.Removed.GATCs == min(results.df.part$Number.of.Removed.GATCs), "Bin.Size"][1]
		AllUsefullData <- sum(DATA[[names(DATA.list.out[[paste("bins", Use.Bin.Size, sep="_")]])[1]]] > 0 & DATA[[names(DATA.list.out[[paste("bins", Use.Bin.Size, sep="_")]])[2]]] > 0)
		UsefullData <- sum(DATA.list.out[[paste("bins", Use.Bin.Size, sep="_")]][, 1] > 0 & DATA.list.out[[paste("bins", Use.Bin.Size, sep="_")]][, 2] > 0)
		statCurrentItem <- data.frame(Sample=i, Original.Correlation=results.df[1, "Pearson.Cor"], Correlation=results.df.part[results.df.part$Number.of.Removed.GATCs == min(results.df.part$Number.of.Removed.GATCs), "Pearson.Cor"], Bin.Size=Use.Bin.Size, Original.NonZero.Value=AllUsefullData, Cleaned.NonZero.Value=UsefullData, Differences=AllUsefullData-UsefullData, Notes=paste("Used replicates: ", names(DATA.filter)[8], ", ", names(DATA.filter)[9], ".", if (length(grep("(.*\\.)3?$", names(DATA.filter)[8:9], perl=T)) != 0) {paste(" And replace replicate number in ", grep("(.*\\.)3?$", names(DATA.filter)[8:9], perl=T, value=T), " from 3 to 2", sep="")}, sep=""))
		plot(x=DATA.list.out[[paste("bins", Use.Bin.Size, sep="_")]][, 1], y=DATA.list.out[[paste("bins", Use.Bin.Size, sep="_")]][, 2], cex=0.3, xlab=names(DATA.list.out[[paste("bins", Use.Bin.Size, sep="_")]])[1], ylab=names(DATA.list.out[[paste("bins", Use.Bin.Size, sep="_")]])[2], las=1, bty="l", pch=".", main=paste("Correlation with bin size", Use.Bin.Size, sep=" "), text(x=min(DATA.list.out[[paste("bins", Use.Bin.Size, sep="_")]][, 1], na.rm=T), y=max(DATA.list.out[[paste("bins", Use.Bin.Size, sep="_")]][, 2], na.rm=T)*0.75, adj=0, labels=c(paste("pearson = ", round(results.df.part[results.df.part$Bin.Size == Use.Bin.Size, "Pearson.Cor"], digits=2), "\n\n", sep=""), paste("spearman = ", round(cor(DATA.list.out[[paste("bins", Use.Bin.Size, sep="_")]][, 1], DATA.list.out[[paste("bins", Use.Bin.Size, sep="_")]][, 2], method="spearman", use="pairwise.complete.obs"), digits=2), sep=""))))
		dev.off()
		if (exists("stat.clean.df") == T) {
			stat.clean.df <- rbind(stat.clean.df, statCurrentItem)
		} else {
			stat.clean.df <- statCurrentItem
		}
		DATA.outfilter <- cbind(DATA.outfilter, DATA.list.out[[paste("bins", Use.Bin.Size, sep="_")]])
	}
}
print("End all filtering data")
write.stat.clean.df <- file.path(outputCleanStat, "Filter_statistics_for_all_samples.csv")
WriteIntermediateFiles(source=stat.clean.df, output.file=write.stat.clean.df)
DATA <- DATA.outfilter

samplesList <- samplesList[samplesList$id %in% names(DATA)[8:length(DATA)], ]
if (length(grep("(.*\\.)3?$", samplesList$id, perl=T)) != 0){
samplesList$id <- sub("(.*\\.)3?$", "\\12", samplesList$id, perl=T)
samplesList$replicate <- sub("3", "2", samplesList$replicate)
names(DATA)[8:length(DATA)] <- sub("(.*\\.)3?$", "\\12", names(DATA)[8:length(DATA)], perl=T)
}
# Combine data into one
#######################
	if (needCombine == T) {
		print("Combine data into one")
		DATA <- cbind(DATA[,1:7], CombineSamples(DATA))
		opt.sum.samples <- paste("DF_Counts_Step_03_Summed_Samples_", currentDate, ".csv", sep="")
		WriteIntermediateFiles(source=DATA, output.file=opt.sum.samples)
	}

# Counts statistics
###################
	print("Counts statistics")
chrs <- unique(DATA$chr)
	DATA.only <- DATA[, 8:ncol(DATA)]
stat <- as.data.frame(matrix(data=NA, nrow=length(chrs), ncol=ncol(DATA.only)+4, byrow=F, dimnames=NULL))
	names(stat) <- c("chr", "GATCs.number", "chr.length.bp", "chr.length.proportion", colnames(DATA.only)[1:(ncol(DATA.only))])
	stat$chr <- chrs
	stat$chr.length.bp <- c(23011544, 21146708, 24543557, 27905053, 1351857, 22422827, 368872, 3288761, 2555491, 2517507, 204112, 347038)
	genome.length <- sum(stat$chr.length.bp)
stat$chr.length.proportion <- round(100 * stat$chr.length.bp / genome.length, digits=2)
	for (j in 1:(ncol(DATA.only))){
		for (i in 1:length(chrs)){
			Data.only.chr <- DATA.only[(DATA$chr == chrs[i]), j]
				if (j == 1) stat$GATCs.number[i] <- length(Data.only.chr)
					stat[i, 4+j] <- sum(Data.only.chr)
					rm(Data.only.chr)
		}
		rm(i)
	}
rm(j)
	statistics.a <- paste("DF_Counts_Step_04_Statistics_A_", currentDate, ".csv", sep="")
WriteIntermediateFiles(source=stat, output.file=statistics.a)
	for (j in 1:(ncol(DATA.only))){
		totalCounts <- sum(stat[, 4+j])
		for (i in 1:length(chrs)){
			stat[i, 4+j] <- round(100 * stat[i, 4+j] / totalCounts, digits=2)
		}
		rm(i, totalCounts)
	}
rm(j)
	statistics.b <- paste("DF_Counts_Step_04_Statistics_B_", currentDate, ".csv", sep="")
WriteIntermediateFiles(source=stat, output.file=statistics.b)

# Add Pseudo counts
###################
DATAs <- list(DATA=DATA)
	if (usePseudoCounts == T) {
		print("Add Pseudo counts")
		for ( i in pseudoCounts) {
			num <- sub("^([0-1]*)(.?)([0-1]*$)", "\\1\\3", i)
			DATA.pseudo <- assign(paste("pseudo", num, sep=""), DATA)
			DATA.pseudo[, 8:ncol(DATA.pseudo)] <- DATA[, 8:ncol(DATA)] + i
			pseudo.filename <- assign(paste("pseudo.fn", num, sep=""), paste("DF_Counts_Step_05_Pseudo_", num, "_Added_", currentDate, ".csv", sep=""))
			DATA.pseudo.strname <- assign(paste("pseudo", num, sep=""), paste("pseudo", num, sep=""))
			WriteIntermediateFiles(source=DATA.pseudo, output.file=pseudo.filename)
			DATAs[[DATA.pseudo.strname]] <- DATA.pseudo
		}
		rm(i)
	}

# Correlation on Counts
#######################
print("Correlation on Counts")
MainCorrelations(dataSet=DATA, corrMethod=corrMethod, suffixCSV="07a_On_Counts", createPDF=F, counts=T)

#################################
# Run many functions from Step_06
#################################

# Declare variables
DATAs.rpm <- DATAs
DATAs.norm <- DATAs
DATAs.norm.ave <- DATAs
####################################
print("Run calculate reads per million")
for (name in names(DATAs.rpm)) {

# Calculation reads per million
###############################
	for (i in 8:(ncol(DATAs.rpm[[name]]))){
		column.sum <- sum(DATAs.rpm[[name]][, i])
			DATAs.rpm[[name]][, i] <- DATAs.rpm[[name]][, i] / column.sum * 10^6
			rm(column.sum)
	}
	rm(i)
	calc.rpm.file <- paste("DF_Counts_Step_06_RPMs_", name, "_", currentDate, ".csv", sep="")
	WriteIntermediateFiles(source=DATAs.rpm[[name]], output.file=calc.rpm.file)

# Correlation on Channels
#########################
		print("Correlation on Channels")
		MainCorrelations(dataSet=DATAs.rpm[[name]], corrMethod=corrMethod, labelHeatmap=labelHeatmap, suffixCSV="07b_On_Channels_A", suffixPDF="07b_On_Channels", corr.on.file=calc.rpm.file, createPDF=T)

# ACF plots for all GATC fragments
##################################
		print("ACF plots for all GATC fragments")
		AcfOnData(dataSet=DATAs.rpm[[name]], labelAcf=labelAcf, method="acf", suffixPDF="08_ACF_Channels", ylab.val="ACF on seq counts", na.data=F)

# Plot boxplots on channels
###########################
		print("Plot boxplots on channels")
		bmp(filename=file.path(prefixDir, paste("DF_Counts_Step_09_Boxplot_Channels_", name, "_", currentDate, ".bmp", sep="")), width=2000, height=1000, units="px")
		par(mar=c(12, 8, 0.5, 0.5))
		boxplot(DATAs.rpm[[name]][, 8:(ncol(DATAs.rpm[[name]]))], names=colnames(DATAs.rpm[[name]])[8:(ncol(DATAs.rpm[[name]]))], las=2, ylab="RPM")
		dev.off()

# DAM Normalization
###################
		print("DAM Normalization")
		DATAs.norm[[name]] <- DATAs.norm[[name]][, -c(8:ncol(DATAs.norm[[name]]))]
		listNorm <- samplesList[1:5]
		listNorm$normalization <- paste(listNorm$tissue, listNorm$conditions, listNorm$replicate, sep="")
		uniqueSamples <- unique(listNorm$normalization)
		for (sample in uniqueSamples) {
			tissue.id <- subset(subset(listNorm, normalization == sample), protein != "DAM")$id
				dam.id <- subset(subset(listNorm, normalization == sample), protein == "DAM")$id
				for (protein in tissue.id) {
					tissue.norm <- paste(protein, ".norm", sep="")
						DATAs.norm[[name]][[tissue.norm]] <- log2(DATAs.rpm[[name]][[protein]] / DATAs.rpm[[name]][[dam.id]])
				}
		}
	for (i in 8:(ncol(DATAs.norm[[name]]))){
		nan.index <- is.nan(DATAs.norm[[name]][, i])
			inf.index <- is.infinite(DATAs.norm[[name]][, i])
			DATAs.norm[[name]][nan.index, i] <- NA
			DATAs.norm[[name]][inf.index, i] <- NA
			rm(nan.index)
			rm(inf.index)
	}
	rm(i)
		dam.norm <- assign(paste("dam.norm", name, sep="."), paste("DF_Counts_Step_10_Dam_norm_", name, "_", currentDate, ".csv", sep=""))
		WriteIntermediateFiles(source=DATAs.norm[[name]], output.file=dam.norm)
} 
rm(name)
############################
# Stop counting from Step_06
############################
	print("Run some functions from normalized data")
#################################
# Run many functions from Step_10
#################################
	for (name in names(DATAs.norm)) {
		print(paste("Start calculate from", name, sep=" "))
# Correlation on Normalized data
################################
			MainCorrelations(dataSet=DATAs.norm[[name]], corrMethod=corrMethod, labelHeatmap=labelHeatmap, use.opt="pairwise.complete.obs", suffixCSV="11_On_Normalized_A", suffixPDF="11_On_Normalized", corr.on.file=dam.norm, createPDF=T)

# Scatter Plots on Normalized data
##################################
			ScatterPlotting3D(dataSet=DATAs.norm[[name]], tag=name)
			print(paste("Stop calculate from", name, sep=" "))
			print(paste("Start calculate from averaged", name, sep=" "))
# Averaging Replicates
######################
			DATAs.norm.ave[[name]] <- DATAs.norm.ave[[name]][, -c(8:ncol(DATAs.norm.ave[[name]]))]
			listNormAve <- listNorm[!(listNorm$protein == "DAM"), ]  # remove row's with DAM
			listNormAve$normalizationAve <- paste(listNormAve$tissue, listNormAve$protein, listNormAve$conditions, sep=".")
			uniqueAveSamples <- unique(listNormAve$normalizationAve)
			for (item in uniqueAveSamples) {
				item.norm.ave <- paste(item, ".norm.ave", sep="")
					vector <- grep(paste(item, "\\.[0-9]", "\\.norm", sep=""), names(DATAs.norm[[name]]))
					weHaveOneReplicate <- length(vector)
					if (weHaveOneReplicate > 1) {
						DATAs.norm.ave[[name]][[item.norm.ave]] <- rowMeans(DATAs.norm[[name]][, c(vector)])
					} else {
						print(paste("This", item, "have one replicate!", sep=" "))
							DATAs.norm.ave[[name]][[item.norm.ave]] <- DATAs.norm[[name]][, c(vector)]
					}	
			}
		dam.norm.ave <- assign(paste("dam.norm.ave", name, sep="."), paste("DF_Counts_Step_12_Averaged_", name, "_", currentDate, ".csv", sep=""))
			WriteIntermediateFiles(source=DATAs.norm.ave[[name]], output.file=dam.norm.ave)

# ACF plots on Averaged
#######################
			print("ACF plots on Averaged")
			AcfOnData(dataSet=DATAs.norm.ave[[name]], labelAcf=labelAcf, method="acf", suffixPDF="14_ACF_Averaged", ylab.val="ACF on rpms", na.data=T)

# Density on Averaged
#######################
			print("Density on Averaged")
			AcfOnData(dataSet=DATAs.norm.ave[[name]], labelAcf=labelAcf, method="density", suffixPDF="15_Density_Averaged", ylab.val="density", na.data=T)

# Create WIG & GFF on Averaged data
###################################
			print("Create WIG & GFF on Averaged data")
			DamIdSeqToWigGff(dataSet=DATAs.norm.ave[[name]])

# If we analyze only one protein that don't need to run Correlation counts and Scatter plots
############################################################################################
			if ((ncol(DATAs.norm.ave[[name]])-7) != 1) {
# Correlations on Averaged NormData
###################################
				print("Correlations on Averaged NormData")
					MainCorrelations(dataSet=DATAs.norm.ave[[name]], corrMethod=corrMethod, labelHeatmap=labelHeatmap, use.opt="pairwise.complete.obs", suffixCSV="13_On_Averaged_A", suffixPDF="13_On_Averaged", corr.on.file=dam.norm.ave, createPDF=T)
			} else {
				print("I found that I calculate only one protein data, because I don't run Correlation counts")
				print(paste("Stop calculate from averaged", name, sep=" "))
			}
	}

###############################
# Start calculate BioHMM data #
###############################

# By sample use only DATA, without pseudo counts

chromosomesVector <- unique(DATAs.norm.ave$DATA$chr)
HMM.data <- list()
DOMAIN.data <- list()
DOMAIN.data.filt <- list()
DOMAIN.data.anti <- list()
for (listItem in 8:length(DATAs.norm.ave$DATA)) {
	DATA.name <- sub("(.*)(_edge|_all)\\.norm\\.ave", "\\1", names(DATAs.norm.ave$DATA)[listItem])
	HMM.data[[DATA.name]] <- list()
	print(paste("Start calculate BioHMM data from", DATA.name, sep=" "))
	for (chr in 1:length(chromosomesVector)) {
		chrCoord <- which(DATAs.norm.ave$DATA$chr == chromosomesVector[chr])
		dataframe.name <- paste("chr", chromosomesVector[chr], sep="_")
		HMM.data[[DATA.name]][[dataframe.name]] <- DATAs.norm.ave$DATA[chrCoord, c(1:7, listItem)]

		names(HMM.data[[DATA.name]][[dataframe.name]])[8] <- "DamID.value"
		HMM.data[[DATA.name]][[dataframe.name]] <- HMM.data[[DATA.name]][[dataframe.name]][!is.na(HMM.data[[DATA.name]][[dataframe.name]]$DamID.value), ]

		print(paste("I'm run BioHMM data from chromosome", chromosomesVector[chr], sep=" "))
		HMM.data[[DATA.name]][[dataframe.name]]$BioHMM.output <- runBioHMM(mval=HMM.data[[DATA.name]][[dataframe.name]]$DamID.value, datainfo=HMM.data[[DATA.name]][[dataframe.name]], useCloneDists=T)

		classifier <- aggregate(x=HMM.data[[DATA.name]][[dataframe.name]]$DamID.value, by=list(HMM.data[[DATA.name]][[dataframe.name]]$BioHMM.output), FUN=mean)
 		
 		if (nrow(classifier) == 3){
			# if "1s" are targets and "2s" are non-targets
			non_value <- (which(classifier == min(classifier$x), arr.ind=T))[1]
			domain_value <- (which(classifier == max(classifier$x), arr.ind=T))[1]
			ambiguous_value <- c(1:3)[!(c(1:3) %in% c(non_value, domain_value))]
			HMM.data[[DATA.name]][[dataframe.name]]$target[(HMM.data[[DATA.name]][[dataframe.name]]$BioHMM.output == non_value)] <- -1
			HMM.data[[DATA.name]][[dataframe.name]]$target[(HMM.data[[DATA.name]][[dataframe.name]]$BioHMM.output == ambiguous_value)] <- 0
			HMM.data[[DATA.name]][[dataframe.name]]$target[(HMM.data[[DATA.name]][[dataframe.name]]$BioHMM.output == domain_value)] <- 1
		} else {
			print(paste(DATA.name, ", ", dataframe.name, " - not enough data!", sep=""))
		}
		if ("target" %in% names(HMM.data[[DATA.name]][[dataframe.name]])){
			tmp_rle <- rle(HMM.data[[DATA.name]][[dataframe.name]]$target)
			tmp_rle$values[tmp_rle$lengths <= 2 & tmp_rle$values == 1] <- 0
			HMM.data[[DATA.name]][[dataframe.name]]$target.filt <- inverse.rle(tmp_rle)

			domains <- FeatureCalls.to.GFF.like(start.coordinate=HMM.data[[DATA.name]][[dataframe.name]]$start, end.coordinate=HMM.data[[DATA.name]][[dataframe.name]]$end, feature.type=HMM.data[[DATA.name]][[dataframe.name]]$target)
			domains <- cbind(chr=HMM.data[[DATA.name]][[dataframe.name]]$chr[1], domains, stringsAsFactors=F)

			antidomains <- domains[domains$value!=1,]
			antidomains$value <- 1
			antidomains <- cbind(antidomains, NA, NA, NA, NA, NA)
			antidomains <- antidomains[,c(1,5,6,2,3,4,7,8,9)]
			names(antidomains) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
			antiprotein.name <- tolower(sub("([a-zA-Z]+)\\.([a-zA-Z0-9]+)\\.([a-zA-Z_0-9]+)", "\\2.anti.domain", DATA.name, perl=T))
			DATA.antidomain.name <- sub("(.*)", "\\1.anti.domains", DATA.name, perl=T)
			antidomains$feature[domains$score!=1] <- antiprotein.name
			antidomains[, c(2,7:9)] <- "."

			domains <- domains[domains$value==1,]
			domains <- cbind(domains, NA, NA, NA, NA, NA)
			domains <- domains[,c(1,5,6,2,3,4,7,8,9)]
			names(domains) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
			protein.name <- tolower(sub("([a-zA-Z]+)\\.([a-zA-Z0-9]+)\\.([a-zA-Z_0-9]+)", "\\2.domain", DATA.name, perl=T))
			DATA.domain.name <- sub("(.*)", "\\1.domains", DATA.name, perl=T)
			domains$feature[domains$score==1] <- protein.name
			domains[, c(2,7:9)] <- "."

			if (chr == 1) {DOMAIN.data[[DATA.domain.name]] <- domains;} else {DOMAIN.data[[DATA.domain.name]] <- rbind(DOMAIN.data[[DATA.domain.name]], domains);}
			if (chr == 1) {DOMAIN.data.anti[[DATA.antidomain.name]] <- antidomains;} else {DOMAIN.data.anti[[DATA.antidomain.name]] <- rbind(DOMAIN.data.anti[[DATA.antidomain.name]], antidomains);}


			domains <- FeatureCalls.to.GFF.like(start.coordinate=HMM.data[[DATA.name]][[dataframe.name]]$start, end.coordinate=HMM.data[[DATA.name]][[dataframe.name]]$end, feature.type=HMM.data[[DATA.name]][[dataframe.name]]$target.filt)
			if (length(domains$value[domains$value==1]) > 0) {
				domains <- cbind(chr=HMM.data[[DATA.name]][[dataframe.name]]$chr[1], domains, stringsAsFactors=F)
				domains <- domains[domains$value==1,]
				domains <- cbind(domains, NA, NA, NA, NA, NA)
				domains <- domains[,c(1,5,6,2,3,4,7,8,9)]
				names(domains) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
				protein.name <- tolower(sub("([a-zA-Z]+)\\.([a-zA-Z0-9]+)\\.([a-zA-Z_0-9]+)", "\\2.domain", DATA.name, perl=T))
				DATA.filtdomain.name <- sub("(.*)", "\\1.domains", DATA.name, perl=T)
				domains$feature[domains$score==1] <- protein.name
				domains[, c(2,7:9)] <- "."
				if ( chr == 1 ) {DOMAIN.data.filt[[DATA.filtdomain.name]] <- domains;} else {DOMAIN.data.filt[[DATA.filtdomain.name]] <- rbind(DOMAIN.data.filt[[DATA.filtdomain.name]], domains);}
			}
		} else {
			print(paste("No data in ", DATA.name, " - ", dataframe.name, ".", sep=""))
		}
	}
	gff.file <- file.path(prefixDir, outputGff, outputDomain, "constitutive", paste(DATA.domain.name, "_constitutive.gff", sep=""))
	gff.file.filt <- file.path(prefixDir, outputGff, outputDomain, "filt", paste(DATA.domain.name, "_filt.gff", sep=""))
	gff.file.anti <- file.path(prefixDir, outputGff, outputDomain, "anti", paste(DATA.domain.name, "_anti.gff", sep=""))

	for (i in c("constitutive", "filt", "anti")) {
		ifelse(!dir.exists(file.path(prefixDir, outputGff, outputDomain, i)), dir.create(file.path(prefixDir, outputGff, outputDomain, i), showWarnings=FALSE), FALSE)
	}

	write.table(DOMAIN.data[[DATA.domain.name]], file=gff.file, sep="\t", row.names=F, col.names=F, quote=F, dec=".", append=F)
	write.table(DOMAIN.data.filt[[DATA.filtdomain.name]], file=gff.file.filt, sep="\t", row.names=F, col.names=F, quote=F, dec=".", append=F)
	write.table(DOMAIN.data.anti[[DATA.antidomain.name]], file=gff.file.anti, sep="\t", row.names=F, col.names=F, quote=F, dec=".", append=F)
}
print("Congratulations!!!")
stop("Do not need to go further", call. = FALSE)
# GFFFormatting(input_df=HMM.data[[DATA.name]][[dataframe.name]], output_df=DOMAIN.data, typeOfFeature=HMM.data[[DATA.name]][[dataframe.name]]$target, data_name=DATA.name, chromosome=chr)

# GFFFormatting <- function(input_df, output_df, typeOfFeature, data_name, chromosome) {
# 	domains <- FeatureCalls.to.GFF.like(start.coordinate=input_df$start, end.coordinate=input_df$end, feature.type=typeOfFeature)
# 	domains <- cbind(chr=input_df$chr[1], domains, stringsAsFactors=F)
# 	domains <- domains[domains$value==1,]
# 	domains <- cbind(domains, NA, NA, NA, NA, NA)
# 	domains <- domains[,c(1,5,6,2,3,4,7,8,9)]
# 	names(domains) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
# 	protein.name <- tolower(sub("([a-zA-Z]+)\\.([a-zA-Z0-9]+)\\.([a-zA-Z_0-9]+)", "\\2.domain", data_name, perl=T))
# 	DATA.domain.name <<- sub("(.*)", "\\1.domains", data_name, perl=T)
# 	domains$feature[domains$score==1] <- protein.name
# 	domains[, c(2,7:9)] <- "."
# 	if (chromosome == 1) {output_df[[DATA.domain.name]] <<- domains;} else {output_df[[DATA.domain.name]] <<- rbind(output_df[[DATA.domain.name]], domains);}
# }

print("Search Biology meaning...")

# Declare functions
###################

MakeReportDomainsSize <- function(inputData, dscr, includeHet=T) {
	domainSizeReport <- as.data.frame(matrix(NA, nrow=0, ncol=5))
	names(domainSizeReport) <- c("Item.name", "Genome.size", "Domain.size", "Ratio", "Selected.chromosomes")
	getChrList <- names(Dmelanogaster)[-c(which(names(Dmelanogaster) == c("chrU", "chrM", "chrUextra")))]
	if (includeHet != T) {
		getChrList <- getChrList[-grep("Het", getChrList)]
		modificator <- "without_Heterochromatin"
	} else {
		modificator <- ""
	}
	chr.lengths <- integer()
	for (i in getChrList) chr.lengths <- append(chr.lengths, length(DNAString(Dmelanogaster[[i]])))
	AllGenome <- sum(chr.lengths)
	getChrList <- sub("chr(.*)", "\\1", getChrList, perl=T)
	for (item in names(inputData)) {
		if (includeHet != T) {
		domainsSize <- sum(inputData[[item]]$end[inputData[[item]]$seqname %in% getChrList]-inputData[[item]]$start[inputData[[item]]$seqname %in% getChrList])
		fifthCol <- paste(unique(inputData[[item]]$seqname[inputData[[item]]$seqname %in% getChrList]), collapse=", ")
		} else {
		domainsSize <- sum(inputData[[item]]$end-inputData[[item]]$start)
		fifthCol <- paste(unique(inputData[[item]]$seqname), collapse=", ")
		}
		partOfDomainFromGenome <- paste(round(pcentFun(domainsSize, AllGenome), digits=3), "%", sep="")
		domainSizeReport[grep(paste("^",item, "$", sep=""), names(inputData)),] <- c(item, AllGenome, domainsSize, partOfDomainFromGenome, fifthCol)
	}
	reportFile <- file.path(prefixDir, outputBio, paste("Report_about", dscr, modificator,"domains_part_from_genome_D.melanogaster.csv", sep="_"))
	write.table(domainSizeReport, file=reportFile, sep=";", row.names=F, col.names=T, quote=F, dec=".", append=F)
	# dfHist <- domainSizeReport[,c(1,4)]
	# dfHist$Ratio <- as.numeric(sub("(.*)%$", "\\1", dfHist$Ratio, perl=T))
	# reportHist <- ggplot(dfHist, aes(x=Item.name, y=Ratio)) + geom_histogram(colour="black", binwidth=0.3)
}



IntersectDomain <- function(inputData, Tset, Pset, dscr, ScoreValue) {
	DOMAIN.data.part <- inputData[grep(paste("(", paste(tissue.bio.set, collapse="|"), ")\\.(", paste(protein.bio.set, collapse="|"), ")\\.(", paste(conditions.bio.set, collapse="|"), ")\\..*", sep="") ,names(inputData), perl=T)]
	COMPARE.data <- list()
	for (i in c(2:5)) {
		if (i == 2) {
			tissue_list <- combn(Tset, i, simplify=F)
		} else {
			tissue_list <- append(tissue_list, combn(Tset, i, simplify=F))
		}
	}
	for (protein in Pset) {
		for (tissueSet in tissue_list) {
			setsSize <- length(tissueSet)
			compare_list <- list()
			chr_compare <- list()
			for (tissue in tissueSet) {
				name <- grep(paste(tissue,"\\.",protein,"\\..*", sep=""),names(DOMAIN.data.part), value=T, perl=T)
				chr_compare[[tissue]] <- unique(DOMAIN.data.part[[name]]$seqname)
			}
			if (setsSize < 2) {
				print("The comparison is not applicable!")
			} else if (setsSize == 2) {
				chr_list <- intersect(chr_compare[[names(chr_compare[1])]], chr_compare[[names(chr_compare[2])]])
			} else if (setsSize == 3) {
				chr_list <- intersect(chr_compare[[names(chr_compare[3])]], intersect(chr_compare[[names(chr_compare[1])]], chr_compare[[names(chr_compare[2])]]))
			} else if (setsSize == 4) {
				chr_list <- intersect(chr_compare[[names(chr_compare[4])]], intersect(chr_compare[[names(chr_compare[3])]], intersect(chr_compare[[names(chr_compare[1])]], chr_compare[[names(chr_compare[2])]])))
			} else if (setsSize == 5){
				chr_list <- intersect(chr_compare[[names(chr_compare[5])]], intersect(chr_compare[[names(chr_compare[4])]], intersect(chr_compare[[names(chr_compare[3])]], intersect(chr_compare[[names(chr_compare[1])]], chr_compare[[names(chr_compare[2])]]))))
			} else {
				print("")
			}
			for (x in chr_list) {
				for (tissue in tissueSet) {
					name <- grep(paste(tissue,"\\.",protein,"\\..*", sep=""),names(DOMAIN.data.part), value=T, perl=T)
					dfIRanges <- DOMAIN.data.part[[name]][DOMAIN.data.part[[name]]$seqname == x, c("start", "end")] 
					objIRanges <- with(dfIRanges, IRanges(start, end))
					compare_list[[name]] <- objIRanges
				}
				if (setsSize < 2) {
				print("The comparison is not applicable!")
				} else if (setsSize == 2) {
					compareResult <- as.data.frame(intersect(compare_list[[names(compare_list[1])]], compare_list[[names(compare_list[2])]]))
				} else if (setsSize == 3) {
					compareResult <- as.data.frame(intersect(compare_list[[names(compare_list[3])]], intersect(compare_list[[names(compare_list[1])]], compare_list[[names(compare_list[2])]])))
				} else if (setsSize == 4) {
					compareResult <- as.data.frame(intersect(compare_list[[names(compare_list[4])]], intersect(compare_list[[names(compare_list[3])]], intersect(compare_list[[names(compare_list[1])]], compare_list[[names(compare_list[2])]]))))
				} else if (setsSize == 5) {
					compareResult <- as.data.frame(intersect(compare_list[[names(compare_list[5])]], intersect(compare_list[[names(compare_list[4])]], intersect(compare_list[[names(compare_list[3])]], intersect(compare_list[[names(compare_list[1])]], compare_list[[names(compare_list[2])]])))))
				}
				if (nrow(compareResult) != 0) {
					compareResult <- compareResult[compareResult$width!=1,]
					compareResult$seqname <- x
					compareResult <- compareResult[, c("seqname", "start", "end")]
					if (exists("compareDf") == T) {
						compareDf <- rbind(compareDf, compareResult)
					} else {
						compareDf <- compareResult
					}
				} else next()
			}
			compareDf <- cbind(compareDf, as.data.frame(matrix(NA, ncol=6, nrow=nrow(compareDf))))
			compareDf <- compareDf[,c(1,5,6,2,3,4,7,8,9)]
			names(compareDf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
			protein.name <- tolower(paste(protein, dscr, "domain", sep="."))
			compareDf$feature <- protein.name
			compareDf$score <- ScoreValue
			compareDf[, c(2,7:9)] <- "."
			compare.domain.name <- paste(sapply(tissueSet, function(tiss) grep(paste(tiss,"\\.",protein,"\\..*", sep=""),names(DOMAIN.data.part), value=T, perl=T)), collapse="_vs_")
			COMPARE.data[[compare.domain.name]] <- compareDf
			gff.file <- file.path(prefixDir, outputGff, outputDomain, dscr, paste(compare.domain.name, dscr, "gff", sep="."))
			ifelse(!dir.exists(file.path(prefixDir, outputGff, outputDomain, dscr)), dir.create(file.path(prefixDir, outputGff, outputDomain, dscr), showWarnings=FALSE), FALSE)
			write.table(COMPARE.data[[compare.domain.name]], file=gff.file, sep="\t", row.names=F, col.names=F, quote=F, dec=".", append=F)
			rm(compareDf)
		}
	}
	return(COMPARE.data)
}

# Heatmap Correlation
kc.list <- list.files("/home/anton/backup/input/KC", full.names=T)
DATA.venn <- DATAs.norm.ave$DATA[, c(1:7)]
DATA.kc <- as.data.frame(matrix(NA, ncol=0, nrow=nrow(DATAs.norm.ave$DATA)))
for (i in kc.list) {
	kc <- read.delim(i, header=T, sep="\t", skip=17, stringsAsFactors=F)
	compare.vector <- match(DATA.venn[,1], kc[,1])
	kc.fullname <- sub("([a-zA-Z0-9]*)-([a-zA-Z0-9]*)_xx.*", paste("\\1", "\\2", "std_all.norm.ave", sep="."), basename(i), perl=T)
	DATA.kc[[kc.fullname]] <- kc$score[compare.vector]
	for (x in c("score", "bound")) {
	kc.name <- sub("([a-zA-Z0-9]*)-([a-zA-Z0-9]*)_xx.*", paste("\\1", "\\2", x, sep="."), basename(i), perl=T)
	DATA.venn[[kc.name]] <- kc[[x]][compare.vector]
	}
}
DATA.middle.part <- cbind(DATAs.norm.ave$DATA[,c(1:7)], DATAs.norm.ave$DATA[,grep(paste("(", paste(tissue.bio.set, collapse="|"), ")\\.(", paste(protein.bio.set, collapse="|"), ")\\.(", paste(conditions.bio.set, collapse="|"), ")_(edge|all).*", sep="") ,names(DATAs.norm.ave$DATA), perl=T, value=T)])

DATA.part <- cbind(DATA.middle.part, DATA.kc)
print("Run heatmap generate. Many NULL message!")
HeatmapBySelection(dataSet=DATA.part, corrMethod=corrMethod, use.opt="pairwise.complete.obs")

# Generate Venn Diagram
for (i in names(DATA.middle.part)[8:ncol(DATA.middle.part)]) {
	name <- sub("(.*)_(all|edge).*", "\\1", i, perl=T)
	name.score <- paste(name, "score", sep=".")
	name.bound <- paste(name, "bound", sep=".")
	DATA.venn[[name.score]] <- DATA.middle.part[[i]]
	for (x in names(HMM.data[[name]])) {
		if ("target" %in% names(HMM.data[[name]][[x]])){
			if (exists("target.df") == F) {
				target.df <- HMM.data[[name]][[x]][, c("ID", "target")]
			} else {
				target.df <- rbind(target.df, HMM.data[[name]][[x]][, c("ID", "target")])
			}
		}
	}
	compare.vector <- match(DATA.venn$ID, target.df$ID)
	DATA.venn[[name.bound]] <- target.df$target[compare.vector]
	DATA.venn[[name.bound]][DATA.venn[[name.bound]] == -1] <- 0
	rm(target.df)
}
print("Make Venn Diagram")
for (t.s in tissue.bio.set)	MakeVennDiagram(x=DATA.venn, set=paste(t.s,".*bound", sep=""), v1=t.s, v2="LAM, HP1, PC")
for (t.s in tissue.bio.set)	MakeVennDiagram(x=DATA.venn, set=paste(t.s,".*\\.(LAM|HP1).*bound", sep=""), v1=t.s, v2="LAM, HP1")
for (p.s in protein.bio.set) MakeVennDiagram(x=DATA.venn, set=paste(".*\\.", p.s, ".*bound", sep=""), v1=p.s, v2="Kc167, BR, FB, NRN, Glia")
for (p.s in protein.bio.set) MakeVennDiagram(x=DATA.venn, set=paste("(Kc167|BR|FB)\\.", p.s, ".*bound", sep=""), v1=p.s, v2="Kc167, BR, FB")

# Count area domain from whole genome
#####################################
KC.DOMAIN <- list()
KC.DOMAIN.ANTI <- list()
for (i in kc.list) {
	kc <- read.delim(i, header=T, sep="\t", skip=17, stringsAsFactors=F)
	name <- sub("([a-zA-Z0-9].+)-([a-zA-Z0-9].+)_xx_.*", "\\1.\\2.m.domains", basename(i) ,perl=T)
	for (y in unique(kc$seqname)) {
		kc_domains <- FeatureCalls.to.GFF.like(start.coordinate=kc$start[kc$seqname == y], end.coordinate=kc$end[kc$seqname == y], feature.type=kc$bound[kc$seqname == y])
		kc_domains <- cbind(chr=y, kc_domains, stringsAsFactors=F)

		kc_anti_domains <- kc_domains[kc_domains$value!=1,]
		kc_anti_domains$value <- 1
		kc_anti_domains <- cbind(kc_anti_domains, NA, NA, NA, NA, NA)
		kc_anti_domains <- kc_anti_domains[,c(1,5,6,2,3,4,7,8,9)]
		names(kc_anti_domains) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
		kc.antidomain.name <- sub("([a-zA-Z0-9].+)-([a-zA-Z0-9].+)_xx_.*", "\\1.\\2.m.anti.domains", basename(i) ,perl=T)
		antiprotein.name <- tolower(kc.antidomain.name)
		kc_anti_domains$feature <- antiprotein.name
		kc_anti_domains[, c(2,7:9)] <- "."
		kc_anti_domains$seqname <- sub("chr(.*)", "\\1", kc_anti_domains$seqname, perl=T)
		KC.DOMAIN.ANTI[[kc.antidomain.name]] <- rbind(KC.DOMAIN.ANTI[[kc.antidomain.name]], kc_anti_domains)

		kc_domains <- kc_domains[kc_domains$value==1,]
		kc_domains <- cbind(kc_domains, NA, NA, NA, NA, NA)
		kc_domains <- kc_domains[,c(1,5,6,2,3,4,7,8,9)]
		names(kc_domains) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
		kc.domain.name <- sub("([a-zA-Z0-9].+)-([a-zA-Z0-9].+)_xx_.*", "\\1.\\2.m.domains", basename(i) ,perl=T)
		protein.name <- tolower(kc.domain.name)
		kc_domains$feature <- protein.name
		kc_domains[, c(2,7:9)] <- "."
		kc_domains$seqname <- sub("chr(.*)", "\\1", kc_domains$seqname, perl=T)
		KC.DOMAIN[[kc.domain.name]] <- rbind(KC.DOMAIN[[kc.domain.name]], kc_domains)
	}
	gff.file <- file.path(prefixDir, outputGff, outputDomain, paste(name, "gff", sep="."))
	gff.file.anti <- file.path(prefixDir, outputGff, outputDomain, paste(name, "anti.gff", sep="."))

	write.table(KC.DOMAIN[[kc.domain.name]], file=gff.file, sep="\t", row.names=F, col.names=F, quote=F, dec=".", append=F)
	write.table(KC.DOMAIN.ANTI[[kc.antidomain.name]], file=gff.file.anti, sep="\t", row.names=F, col.names=F, quote=F, dec=".", append=F)

}
DOMAIN.data <- append(DOMAIN.data, KC.DOMAIN)
DOMAIN.data.filt <- append(DOMAIN.data.filt, KC.DOMAIN)
DOMAIN.data.anti <- append(DOMAIN.data.anti, KC.DOMAIN.ANTI)

MakeReportDomainsSize(DOMAIN.data, "constitutive")
MakeReportDomainsSize(DOMAIN.data, "constitutive", includeHet=F)

MakeReportDomainsSize(DOMAIN.data.filt, "filtered")
MakeReportDomainsSize(DOMAIN.data.filt, "filtered", includeHet=F)

MakeReportDomainsSize(DOMAIN.data.anti, "anti")
MakeReportDomainsSize(DOMAIN.data.anti, "anti", includeHet=F)

# Make boxplot & histogram from BioHMM data
for (i in names(HMM.data)) {
	for (x in names(HMM.data[[i]])){
		setDomain <- unique(HMM.data[[i]][[x]]$BioHMM.output)
		if ("target" %in% names(HMM.data[[i]][[x]])){
			if (exists("vecHist2") == F) {
				tmp <- rle(HMM.data[[i]][[x]]$target)
				vecHist2 <- tmp$length[tmp$values == 1]
			} else {
				tmp <- rle(HMM.data[[i]][[x]]$target)
				vecHist2 <- append(vecHist2, tmp$length[tmp$values == 1])
			}
		}
		if (length(setDomain) < 3) {
			if (length(setDomain) == 2) {
				if (all(sort(setDomain) == c(1,3))) labelsDomain <- c("-1", "1")
				if (all(sort(setDomain) == c(1,2))) labelsDomain <- c("-1", "0")
				if (all(sort(setDomain) == c(2,3))) labelsDomain <- c("0", "1")
				} else {
					if (setDomain == 1) labelsDomain <- "-1"
					if (setDomain == 2) labelsDomain <- "0"
					if (setDomain == 3) labelsDomain <- "1"
				}
		} else {
			labelsDomain <- c("-1", "0", "1")
		}
		if(exists("dfBox") == F) {
			dfBox <- data.frame("DamID.value"=HMM.data[[i]][[x]]$DamID.value, "BioHMM.output"=factor(HMM.data[[i]][[x]]$BioHMM.output, labels=labelsDomain), "chr"=factor(HMM.data[[i]][[x]]$chr, labels=x))
		} else {
			dfBox <- rbind(dfBox, data.frame("DamID.value"=HMM.data[[i]][[x]]$DamID.value, "BioHMM.output"=factor(HMM.data[[i]][[x]]$BioHMM.output, labels=labelsDomain), "chr"=factor(HMM.data[[i]][[x]]$chr, labels=x)))
		}
	}
	plot1 <- ggplot(aes(y = DamID.value, x = chr, fill = BioHMM.output), data = dfBox) + geom_boxplot() + labs(title=paste(i, "BioHMM output", sep="\n"))
	domain.name <- paste(i, "domains", sep=".")
	vecHist <- (DOMAIN.data[[domain.name]]$end - DOMAIN.data[[domain.name]]$start)/1000
	plot2 <- ggplot() + aes(vecHist)
	if (max(vecHist) >= 300) {
		plot2 <- plot2 + geom_histogram(binwidth = 1, colour="black", fill="white") + coord_cartesian(xlim = c(0, 300)) + labs(x=paste("Domain size in kb. Image has been scaled. Max value is equal", max(vecHist), sep=" "))
	} else {
		plot2 <- plot2 + geom_histogram(binwidth = max(vecHist)/300, colour="black", fill="white") + coord_cartesian(xlim = c(0, max(vecHist))) + labs(x="Domain size in kb")
	}
	plot3 <- ggplot() + aes(vecHist2)
	if (max(vecHist2) >= 300) {
		plot3 <- plot3 + geom_histogram(binwidth = 1, colour="black", fill="white") + coord_cartesian(xlim = c(0, 300)) + labs(x=paste("Domain size in GATC fragments. Image has been scaled. Max value is equal", max(vecHist2), sep=" "))
	} else {
		plot3 <- plot3 + geom_histogram(binwidth = max(vecHist2)/300, colour="black", fill="white") + coord_cartesian(xlim = c(0, max(vecHist2))) + labs(x="Domain size in GATC fragments")
	}
	bmp(filename=file.path(prefixDir, outputBio, outputBioBoxplot, paste("Graph_on_BioHMM_Data_from", i, ".bmp", sep="")), width=1920, height=1080, units = "px")
	print(multiplot(plot1, plot2, plot3))
	dev.off()
	rm(dfBox, vecHist2)
}
print("End of bio-analyzing")

# Domain intersection

COMPARE.data <- IntersectDomain(inputData=DOMAIN.data, Tset=tissueIntersectSet, Pset=proteinIntersectSet, dscr="constitutive", ScoreValue=1)
COMPARE.data.filt <- IntersectDomain(inputData=DOMAIN.data.filt, Tset=tissueIntersectSet, Pset=proteinIntersectSet, dscr="filt", ScoreValue=1)
COMPARE.data.antidomains <- IntersectDomain(inputData=DOMAIN.data.anti, Tset=tissueIntersectSet, Pset=proteinIntersectSet, dscr="anti", ScoreValue=-1)

# Write stistics from compared data
MakeReportDomainsSize(COMPARE.data, "compare_constitutive")
MakeReportDomainsSize(COMPARE.data.antidomains, "compare_anti")

#Make histogramm from compared data
for (i in names(COMPARE.data)) {
	vecHist <- (COMPARE.data[[i]]$end - COMPARE.data[[i]]$start)/1000
	plot2 <- ggplot() + aes(vecHist)
	if (max(vecHist) >= 300) {
		plot2 <- plot2 + geom_histogram(binwidth = 1, colour="black", fill="white") + coord_cartesian(xlim = c(0, 300)) + labs(x=paste("Domain size in kb. Image has been scaled. Max value is equal", max(vecHist), sep=" "))
	} else {
		plot2 <- plot2 + geom_histogram(binwidth = max(vecHist)/300, colour="black", fill="white") + coord_cartesian(xlim = c(0, max(vecHist))) + labs(x="Domain size in kb")
	}
	bmp(filename=file.path(prefixDir, outputBio, outputBioBoxplot, paste("Histogram_on_Compared_BioHMM_Data_from", i, ".bmp", sep="")), width=1920, height=1080, units = "px")
	print(plot2)
	dev.off()
}

#######################################################################
############## GATCs overlaps Genes ###################################
#######################################################################
AssignGATCsToGenes <- function(Genes, GATCs){
Results <- as.data.frame(nonzero(regionOverlap(Genes, GATCs)))
colnames(Results) <- c("RowNumber", "ColNumber")
CompareResults <- cbind(GATCs[Results$ColNumber, c("ID","chr","start","end")], data.frame("Genes.ID" = Genes$ID[Results$RowNumber], "TSS" = Genes$start[Results$RowNumber], "strand" = Genes$strand[Results$RowNumber]))
rownames(CompareResults) <- NULL
return(CompareResults)
}

# Intersect with promoters only
gatcG <- read.delim(gatcFile4Genes, header=T, sep="\t", stringsAsFactors=F)
Genes <- read.delim(genesFilePath, header=T, as.is=T, dec=".")
Genes <- Genes[Genes$chr %in% intersect(unique(Genes$chr), unique(gatcG$chr)),]
gatcG <- gatcG[gatcG$chr %in% intersect(unique(Genes$chr), unique(gatcG$chr)),]
Genes[Genes$strand == "+", "end"] <- Genes[Genes$strand == "+", "start"]
Genes[Genes$strand == "-", "start"] <- Genes[Genes$strand == "-", "end"]
GATCvsGenes <- AssignGATCsToGenes(Genes, gatcG)
idx.gatc.ori <- unlist(sapply(GATCvsGenes$ID, function(x) grep(x, gatcG$ID)))
GATCvsGenes.List <- list()

# Intersect with gene body
GenesOrig <- read.delim(genesFilePath, header=T, as.is=T, dec=".")
GATCvsGenesOrig <- AssignGATCsToGenes(GenesOrig, gatcG)
idx.gatc.ori.body <- unlist(sapply(GATCvsGenesOrig$ID, function(x) grep(x, gatcG$ID)))
GATCvsGenesOrig.List <- list()


for (i in tissueExprSet) {
	tissueSelection <- list.files(path=exprDataDir, pattern=paste(i, "_.*", sep=""), full.names=T)
	if (length(tissueSelection) > 1) {
		tissueNames <- file_path_sans_ext(list.files(path=exprDataDir, pattern=paste(i, "_.*", sep="")))
		for (y in tissueSelection) {
			exprTissue <- read.delim(y, header=F, stringsAsFactors=F, sep="\t")
			if (exists("rnaSeq")==F) {
				rnaSeq <- data.frame(exprTissue[, 2])
			} else {
				rnaSeq <- cbind(rnaSeq, exprTissue[, 2])
			}
		}
		names(rnaSeq) <- tissueNames
		rownames(rnaSeq) <- exprTissue[,1]
		if (sum(is.na(rnaSeq)) != 0) {
			rnaSeq[is.na(rnaSeq),] <- 0
	 	}
		conditions <- rep(i, 2)
		exprData <-newCountDataSet(countData=rnaSeq, conditions=conditions)
		exprData <- estimateSizeFactors(exprData)
		exprData <- estimateDispersions(exprData)
		exprResults <- nbinomTest(exprData, i, i)
	} else {
		exprResults <- read.delim(tissueSelection, header=F, stringsAsFactors=F, sep="\t")
		names(exprResults) <- c("id", "baseMean")
	}
	tmpDF <- GATCvsGenes
	tmpDF.Orig <- GATCvsGenesOrig
	if (i == "Kc167") {
		tempDATA <- DATA.part 
	} else {
		tempDATA <- DATAs.norm.ave$DATA
	}
	for (z in protein.bio.set) {
		dataSet <- grep(paste(i, "\\.", z, "\\.(", paste(conditions.bio.set, collapse="|"), ")_(all|edge).*",sep=""), names(tempDATA), value=T, perl=T)
		df <- tempDATA[, c("ID", dataSet)]
		if (i == "Kc167") {
			hmmSet <- grep(paste(i, "\\.", z, "\\.bound", sep=""), names(DATA.venn), value=T, perl=T)
			hmm.df <- DATA.venn[, c(1:7, grep(hmmSet, names(DATA.venn)))]
			hmm.df <- hmm.df[!(is.na(hmm.df[,ncol(hmm.df)])),]
			names(hmm.df)[ncol(hmm.df)] <- "target"
		} else {
			hmmSet <- grep(paste(i, "\\.", z, "\\.(", paste(conditions.bio.set, collapse="|"), ")$",sep=""), names(HMM.data), value=T, perl=T)
			hmm.df <- rbind.fill(HMM.data[[hmmSet]])
		}
		idx.hmm <- which(df$ID %in% hmm.df$ID)
		df <- cbind(df, as.data.frame(matrix(NA, ncol=2, nrow=nrow(df))))
		names(df)[3:4] <- paste(hmmSet, c("target", "target.filt"), sep="_")
		df[,3][idx.hmm] <- hmm.df$target
		if (i != "Kc167") {
			df[,4][idx.hmm] <- hmm.df$target.filt
		}
		tmpDF <- cbind(tmpDF, df[idx.gatc.ori, c(2:4)])
		tmpDF.Orig <- cbind(tmpDF.Orig, df[idx.gatc.ori.body, 2])
		names(tmpDF.Orig)[ncol(tmpDF.Orig)] <- names(df)[2]
	}
	idx.fbgn.ori <- match(tmpDF$Genes.ID, exprResults$id)
	tmpDF$expression <- exprResults$baseMean[idx.fbgn.ori]
	GATCvsGenes.List[[i]] <- tmpDF

	tmpDF.Orig <- tmpDF.Orig[, -c(1:4, 6, 7)]
	tmpDF.Orig <- ddply(tmpDF.Orig,"Genes.ID",numcolwise(mean, na.rm = TRUE))
	for (colItem in c(2:ncol(tmpDF.Orig))) tmpDF.Orig[is.nan(tmpDF.Orig[,colItem]) == T, colItem] <- NA
	idx.fbgn.ori.body <- match(tmpDF.Orig$Genes.ID, exprResults$id)
	tmpDF.Orig$expression <- exprResults$baseMean[idx.fbgn.ori.body]
	GATCvsGenesOrig.List[[i]] <- tmpDF.Orig

	rm(rnaSeq, tmpDF, tmpDF.Orig)
}

# selector <- grep("WID\\.PC\\.(Late|Early)", names(HMM.data), value=T)

# selector <- grep("SG\\.LAM\\.wt", names(HMM.data), value=T)

# selector <- grep("BR\\.LAM\\.m$", names(HMM.data), value=T)
# for (wid in selector) {
# 	widdf <- rbind.fill(HMM.data[[wid]])
# 	dataSet <- grep(paste(wid, "_(all|edge).*", sep=""), names(DATAs.norm.ave$DATA))
# 	gatcdf <- DATAs.norm.ave$DATA[, c(1:7, dataSet)]
# 	index <- which(gatcdf$ID %in% widdf$ID)
# 	gatcdf <- cbind(gatcdf, as.data.frame(matrix(NA, ncol=3, nrow=nrow(gatcdf))))
# 	names(gatcdf)[c(1,9:11)] <- c("GATC.ID", "BioHMM.output", "target", "target.filt")
# 	gatcdf$BioHMM.output[index] <- widdf$BioHMM.output
# 	gatcdf$target[index] <- widdf$target
# 	gatcdf$target.filt[index] <- widdf$target.filt
# 	# widFile <- file.path(prefixDir, outputBio, paste("USA_HMM_data_3state_from_", wid, ".csv", sep=""))
# 	# write.table(gatcdf, file=widFile, sep=";", row.names=F, col.names=T, quote=F, dec=".", append=F)
# }

# Count statistic

Stat_BoxList <- list()
Stat_BoxList_3ds <- list()
for (i in names(GATCvsGenes.List)) {
	for (x in protein.bio.set) {
		target.pos <- grep(paste(".*", x, ".*target$", sep=""), names(GATCvsGenes.List[[i]]), value=T)
		expr.pos <- ncol(GATCvsGenes.List[[i]])
		vectorPlus <- GATCvsGenes.List[[i]][[target.pos]] == 1 & !(is.na(GATCvsGenes.List[[i]][[target.pos]]))
		vectorMinus <- GATCvsGenes.List[[i]][[target.pos]] == -1 & !(is.na(GATCvsGenes.List[[i]][[target.pos]]))
		vectorZero <- GATCvsGenes.List[[i]][[target.pos]] == 0 & !(is.na(GATCvsGenes.List[[i]][[target.pos]]))

		vectorZeroMinus <- (GATCvsGenes.List[[i]][[target.pos]] == -1 | GATCvsGenes.List[[i]][[target.pos]] == 0) & !(is.na(GATCvsGenes.List[[i]][[target.pos]]))

		target_one <- GATCvsGenes.List[[i]][vectorPlus, target.pos]
		target_zero_minus <- rep(0, length(vectorZeroMinus[vectorZeroMinus==T]))

		target_minus <- GATCvsGenes.List[[i]][vectorMinus, target.pos]
		target_zero <-GATCvsGenes.List[[i]][vectorZero, target.pos]

		expr_one <- GATCvsGenes.List[[i]][vectorPlus, expr.pos]
		expr_zero_minus <- GATCvsGenes.List[[i]][vectorZeroMinus, expr.pos]

		expr_minus <- GATCvsGenes.List[[i]][vectorMinus, expr.pos]
		expr_zero <-GATCvsGenes.List[[i]][vectorZero, expr.pos]

		if (exists("dfBox") == F) {
			dfBox <- data.frame("protein" = x, "bind" = c(target_one, target_zero_minus), "expression" = c(expr_one, expr_zero_minus))
			dfBox3st <- data.frame("protein" = x, "bind" = c(target_one, target_zero, target_minus), "expression" = c(expr_one, expr_zero, expr_minus))
		} else {
			dfBox <- rbind(dfBox, data.frame("protein" = x, "bind" = c(target_one, target_zero_minus), "expression" = c(expr_one, expr_zero_minus)))
			dfBox3st <- rbind(dfBox3st, data.frame("protein" = x, "bind" = c(target_one, target_zero, target_minus), "expression" = c(expr_one, expr_zero, expr_minus)))
		}
	}
	boxplot_folder <- file.path(prefixDir, outputBio, outputExpr)
	ifelse(!dir.exists(file.path(prefixDir, outputBio, outputExpr)), dir.create(file.path(prefixDir, outputBio, outputExpr), showWarnings=FALSE), FALSE)
	dfBox$bind <- factor(dfBox$bind, labels=c("not bind", "bind"))
	Stat_BoxList[[i]] <- dfBox
	Stat_BoxList_3ds[[i]] <- dfBox3st
	dfBox$expression <- log2(dfBox$expression + 1)

	if (length(unique(dfBox3st$bind)) < 3) {
		dfBox3st$bind <- factor(dfBox3st$bind, labels=c("ambiguous", "bind"))
	} else {
		dfBox3st$bind <- factor(dfBox3st$bind, labels=c("not bind", "ambiguous", "bind"))
	}
	dfBox3st$expression <- log2(dfBox3st$expression + 1)

	plot_box <- ggplot(aes(y = expression, x = bind, fill=bind), data = dfBox) + geom_boxplot() + labs(title=paste(i, " tissue.", "\nBioHMM output", sep="")) + facet_grid(. ~ protein)

	# plot_box <- ggplot(aes(y = expression, x = bind, fill=bind), data = dfBox) + geom_boxplot() + facet_grid(. ~ protein) + theme(legend.position = "none") + theme(strip.text.x = element_text(size=25))

	# plot_box2 <- plot_box + theme(legend.position = "none") + theme(strip.text.x = element_text(size=25), axis.text=element_text(size=16), axis.title.x=element_blank(), axis.title.y=element_text(size=22)) + ylab("Genes log2 (RNA tag count + 1)")
	# plot_box2 + scale_x_continuous(label=c("TSS not bound", "TSS bound"))


	plot_box3st <- ggplot(aes(y = expression, x = bind, fill=bind), data = dfBox3st) + geom_boxplot() + labs(title=paste(i, " tissue.", "\nBioHMM output. 3 state interpretation.", sep="")) + facet_grid(. ~ protein)
	plot_hist <- ggplot(dfBox, aes(x=expression, fill=bind)) + geom_histogram(colour="black", binwidth=0.3) + scale_y_sqrt() + facet_grid(bind ~ protein)
	pdf(file=file.path(boxplot_folder, paste("Boxplot_expression_data_with_binding_from_", i, ".pdf", sep="")), width=14, height=14)
	print(plot_box)
	dev.off()
	pdf(file=file.path(boxplot_folder, paste("3State_Boxplot_expression_data_with_binding_from_", i, ".pdf", sep="")), width=14, height=14)
	print(plot_box3st)
	dev.off()
	pdf(file=file.path(boxplot_folder, paste("Histogram_expression_data_with_binding_from_", i, ".pdf", sep="")), width=18, height=12)
	print(plot_hist)
	dev.off()
	rm(dfBox, dfBox3st)
}

# Count P-value by wilcox-text
by(Stat_BoxList$BR$expression, INDICES = list(Stat_BoxList$BR$protein, Stat_BoxList$BR$bind), wilcox.test)

for (i in names(Stat_BoxList)) {
	for (y in levels(Stat_BoxList[[i]]$protein)) {
		wilcoxStat <- wilcox.test(Stat_BoxList[[i]][Stat_BoxList[[i]]$protein == y & Stat_BoxList[[i]]$bind == "bind", "expression"], Stat_BoxList[[i]][Stat_BoxList[[i]]$protein == y & Stat_BoxList[[i]]$bind == "not bind", "expression"])
		print(paste(y, " in ", i, ". P-value: ", wilcoxStat$p.value, sep=""))
	}
}

#  for (i in c("FB", "Glia", "NRN")) {
#  	Stat_BoxList_3ds[[i]]$bind <- factor(Stat_BoxList_3ds[[i]]$bind, labels=c("not bind", "ambiguous", "bind"))
#  }
# Stat_BoxList_3ds$Kc167$bind <- factor(Stat_BoxList_3ds$Kc167$bind, labels=c("not bind", "bind"))
for (i in names(Stat_BoxList_3ds)) {
	for (y in levels(Stat_BoxList_3ds[[i]]$protein)) {
		bind_vs_notbind <- wilcox.test(Stat_BoxList_3ds[[i]][Stat_BoxList_3ds[[i]]$protein == y & Stat_BoxList_3ds[[i]]$bind == "bind", "expression"], Stat_BoxList_3ds[[i]][Stat_BoxList_3ds[[i]]$protein == y & Stat_BoxList_3ds[[i]]$bind == "not bind", "expression"])
		bind_vs_ambiguous <- wilcox.test(Stat_BoxList_3ds[[i]][Stat_BoxList_3ds[[i]]$protein == y & Stat_BoxList_3ds[[i]]$bind == "bind", "expression"], Stat_BoxList_3ds[[i]][Stat_BoxList_3ds[[i]]$protein == y & Stat_BoxList_3ds[[i]]$bind == "ambiguous", "expression"])
		notbind_vs_ambiguous <- wilcox.test(Stat_BoxList_3ds[[i]][Stat_BoxList_3ds[[i]]$protein == y & Stat_BoxList_3ds[[i]]$bind == "not bind", "expression"], Stat_BoxList_3ds[[i]][Stat_BoxList_3ds[[i]]$protein == y & Stat_BoxList_3ds[[i]]$bind == "ambiguous", "expression"])
		print(paste("Bind vs not bind. ", y, " in ", i, ". P-value: ", bind_vs_notbind$p.value, sep=""))
		print(paste("Bind vs ambiguous. ", y, " in ", i, ". P-value: ", bind_vs_ambiguous$p.value, sep=""))
		print(paste("Not bind vs ambiguous. ", y, " in ", i, ". P-value: ", notbind_vs_ambiguous$p.value, sep=""))
	}
}



# intersect domain data and combination with expression dataframe
# Genes variable report us expression data, DOMAIN.data and COMPARE.data as well as their "anti" analogues report us coordinates of domains
# chromos <- "2R"
# usePandT <- "BR.LAM.m.anti.domains"
# sampleGenes <- with(Genes[Genes$chr == chromos,], IRanges(start, end))
# sampleDomains <- with(DOMAIN.data.anti[[usePandT]][DOMAIN.data.anti[[usePandT]]$seqname == chromos,], IRanges(start, end))
# Genes_vs_Domains <- as.data.frame(intersect(sampleGenes, sampleDomains))
# shareHKG <- paste(round((nrow(Genes_vs_Domains)/nrow(Genes))*100, digits=2), "%", sep="")
# print(shareHKG)


# chromos <- "2R"
# usePandT <- "BR.LAM.m.domains"
# sampleGenes <- with(Genes[Genes$chr == chromos,], IRanges(start, end))
# sampleDomains <- with(DOMAIN.data[[usePandT]][DOMAIN.data[[usePandT]]$seqname == chromos,], IRanges(start, end))
# Genes_vs_Domains <- as.data.frame(intersect(sampleGenes, sampleDomains))
# shareHKG <- paste(round((nrow(Genes_vs_Domains)/nrow(Genes))*100, digits=2), "%", sep="")
# print(shareHKG)

#Make Scatter plots between gene expression and DamId targets

# MakeScatterPlotsExpressionDamID <- function(DataList, tissueRange, label) {
# 	for (i in tissueRange) {
# 		selectColumns <- grep(paste(i, ".*norm\\.ave", sep="") ,names(DataList[[i]]), perl=T)
# 		ifelse(!dir.exists(file.path(prefixDir, outputBio, outputScttr)), dir.create(file.path(prefixDir, outputBio, outputScttr), showWarnings=FALSE), FALSE)
# 		for (y in selectColumns) {
# 			scatterName <- sub(paste("(",i, "\\.)(LAM|HP1|PC).*", sep=""), "\\1\\2", names(DataList[[i]][y]))
# 			Cor.P <- round(cor(DataList[[i]][,y], DataList[[i]][,ncol(DataList[[i]])], method="pearson", use="pairwise.complete.obs"), digits=2)
# 			Cor.S <- round(cor(DataList[[i]][,y], DataList[[i]][,ncol(DataList[[i]])], method="spearman", use="pairwise.complete.obs"), digits=2)
# 			# plot(x=DataList[[i]][,y], y=log2(DataList[[i]][,ncol(DataList[[i]])]), cex=0.3, xlab=scatterName, ylab="Genes expression value log2(RNA tag count)", text(x=max(DataList[[i]][,y], na.rm=T)*0.85, y=max(log2(DataList[[i]][,ncol(DataList[[i]])]), na.rm=T)*0.8, labels=c(paste("r = ", Cor.P, "\n\n", sep=""), paste("s = ", Cor.S, sep=""))))
# 			scatter_plot_one <- ggplot(DataList[[i]], aes(DataList[[i]][,y], log2(DataList[[i]][,ncol(DataList[[i]])])))+geom_point(alpha=1/10, colour="red", size=3) + xlab(scatterName) + ylab("Genes expression value log2(RNA tag count)") + geom_text(data = data.frame(), size = 4, hjust=0, aes(max(DataList[[i]][,y], na.rm=T)*0.65, max(log2(DataList[[i]][,ncol(DataList[[i]])]), na.rm=T)*0.9, label = c(paste("Pearson.Cor = ", Cor.P, "\n\n", sep=""), paste("Spearman.Cor = ", Cor.S, sep="")))) + theme_bw()
# 			bmp(filename=file.path(prefixDir, outputBio, outputScttr, paste("Scatter_plots_", scatterName, "_vs_expression_values. Type_", label,".bmp", sep="")), width=800, height=800, units = "px")
# 			par(mfrow=c(1, 1))
# 			par(mai=c(1.5, 1.5, 0.7, 0.5))
# 			par(cex=1.5)
# 			print(scatter_plot_one)
# 			dev.off()
# 			rm(scatter_plot_one)
# 		}
# 	}
# }

# MakeScatterPlotsExpressionDamID(DataList=GATCvsGenes.List, tissueRange=c("BR", "FB", "Kc167"), label="by_promotors")
# MakeScatterPlotsExpressionDamID(DataList=GATCvsGenesOrig.List, tissueRange=c("BR", "FB", "Kc167"), label="by_full_gene_body")

for (i in c("BR", "FB", "Kc167")) {
	selectColumns <- grep(paste(i, ".*norm\\.ave", sep="") ,names(GATCvsGenes.List[[i]]), perl=T)
	ifelse(!dir.exists(file.path(prefixDir, outputBio, outputScttr)), dir.create(file.path(prefixDir, outputBio, outputScttr), showWarnings=FALSE), FALSE)
	for (y in selectColumns) {
		scatterName <- sub(paste("(",i, "\\.)(LAM|HP1|PC).*", sep=""), "\\1\\2", names(GATCvsGenes.List[[i]][y]))
		Cor.P <- round(cor(GATCvsGenes.List[[i]][,y], GATCvsGenes.List[[i]][,ncol(GATCvsGenes.List[[i]])], method="pearson", use="pairwise.complete.obs"), digits=2)
		Cor.S <- round(cor(GATCvsGenes.List[[i]][,y], GATCvsGenes.List[[i]][,ncol(GATCvsGenes.List[[i]])], method="spearman", use="pairwise.complete.obs"), digits=2)
		# plot(x=GATCvsGenes.List[[i]][,y], y=log2(GATCvsGenes.List[[i]][,ncol(GATCvsGenes.List[[i]])]), cex=0.3, xlab=scatterName, ylab="Genes expression value log2(RNA tag count)", text(x=max(GATCvsGenes.List[[i]][,y], na.rm=T)*0.85, y=max(log2(GATCvsGenes.List[[i]][,ncol(GATCvsGenes.List[[i]])]), na.rm=T)*0.8, labels=c(paste("r = ", Cor.P, "\n\n", sep=""), paste("s = ", Cor.S, sep=""))))
		scatter_plot_one <- ggplot(GATCvsGenes.List[[i]], aes(GATCvsGenes.List[[i]][,y], log2(GATCvsGenes.List[[i]][,ncol(GATCvsGenes.List[[i]])])))+geom_point(alpha=1/10, colour="red", size=3) + xlab(scatterName) + ylab("Genes expression value log2(RNA tag count)") + geom_text(data = data.frame(), size = 4, hjust=0, aes(max(GATCvsGenes.List[[i]][,y], na.rm=T)*0.65, max(log2(GATCvsGenes.List[[i]][,ncol(GATCvsGenes.List[[i]])]), na.rm=T)*0.9, label =c(paste("Pearson.Cor = ", Cor.P, "\n\n", sep=""), paste("Spearman.Cor = ", Cor.S, sep="")))) + theme_bw()
		bmp(filename=file.path(prefixDir, outputBio, outputScttr, paste("Scatter_plots_", scatterName, "_vs_expression_values. Type_by_promotors.bmp", sep="")), width=800, height=800, units = "px")
		par(mfrow=c(1, 1))
		par(mai=c(1.5, 1.5, 0.7, 0.5))
		par(cex=1.5)
		print(scatter_plot_one)
		rm(Cor.P, Cor.S, scatter_plot_one)
		dev.off()
	}
}

for (i in c("BR", "FB", "Kc167")) {
	selectColumns <- grep(paste(i, ".*norm\\.ave", sep="") ,names(GATCvsGenesOrig.List[[i]]), perl=T)
	ifelse(!dir.exists(file.path(prefixDir, outputBio, outputScttr)), dir.create(file.path(prefixDir, outputBio, outputScttr), showWarnings=FALSE), FALSE)
	for (y in selectColumns) {
		scatterName <- sub(paste("(",i, "\\.)(LAM|HP1|PC).*", sep=""), "\\1\\2", names(GATCvsGenesOrig.List[[i]][y]))
		Cor.P <- round(cor(GATCvsGenesOrig.List[[i]][,y], GATCvsGenesOrig.List[[i]][,ncol(GATCvsGenesOrig.List[[i]])], method="pearson", use="pairwise.complete.obs"), digits=2)
		Cor.S <- round(cor(GATCvsGenesOrig.List[[i]][,y], GATCvsGenesOrig.List[[i]][,ncol(GATCvsGenesOrig.List[[i]])], method="spearman", use="pairwise.complete.obs"), digits=2)
		# plot(x=GATCvsGenesOrig.List[[i]][,y], y=log2(GATCvsGenesOrig.List[[i]][,ncol(GATCvsGenesOrig.List[[i]])]), cex=0.3, xlab=scatterName, ylab="Genes expression value log2(RNA tag count)", text(x=max(GATCvsGenesOrig.List[[i]][,y], na.rm=T)*0.85, y=max(log2(GATCvsGenesOrig.List[[i]][,ncol(GATCvsGenesOrig.List[[i]])]), na.rm=T)*0.8, labels=c(paste("r = ", Cor.P, "\n\n", sep=""), paste("s = ", Cor.S, sep=""))))
		scatter_plot_one <- ggplot(GATCvsGenesOrig.List[[i]], aes(GATCvsGenesOrig.List[[i]][,y], log2(GATCvsGenesOrig.List[[i]][,ncol(GATCvsGenesOrig.List[[i]])])))+geom_point(alpha=1/10, colour="red", size=3) + xlab(scatterName) + ylab("Genes expression value log2(RNA tag count)") + geom_text(data = data.frame(), size = 4, hjust=0, aes(max(GATCvsGenesOrig.List[[i]][,y], na.rm=T)*0.65, max(log2(GATCvsGenesOrig.List[[i]][,ncol(GATCvsGenesOrig.List[[i]])]), na.rm=T)*0.9, label =c(paste("Pearson.Cor = ", Cor.P, "\n\n", sep=""), paste("Spearman.Cor = ", Cor.S, sep="")))) + theme_bw()
		bmp(filename=file.path(prefixDir, outputBio, outputScttr, paste("Scatter_plots_", scatterName, "_vs_expression_values. Type_by_full_gene_body.bmp", sep="")), width=800, height=800, units = "px")
		par(mfrow=c(1, 1))
		par(mai=c(1.5, 1.5, 0.7, 0.5))
		par(cex=1.5)
		print(scatter_plot_one)
		rm(Cor.P, Cor.S, scatter_plot_one)
		dev.off()
	}
}


SplitDfFromColToRow <- function (data) {
	if (ncol(data) %% 2 == 0) {
		ColCounts <- ncol(data) / 2
		dfTempList <- as.list(rnorm(ColCounts))
		names(dfTempList) <- paste("a", 1:length(dfTempList), sep = "")
		for (envDF in c(1:ColCounts)) {
			tempDF <- data[,c(envDF*2-1, envDF*2)]
			tempDF$tissue <- sub("([a-zA-Z0-9]+)\\..*target$", "\\1", names(tempDF)[1], perl=T)
			names(tempDF)[1:2] <- c("bound", "expression")
			dfTempList[[names(dfTempList)[envDF]]] <- tempDF
			rm(tempDF)
		}
		return(rbind.fill(dfTempList))
	}
}
# Find common GATC's by protein
for (protein in proteinIntersectSet){
	for (i in names(GATCvsGenes.List)) {
		vecSel <- grep(paste(i, "\\.", protein, "\\..*target$", sep=""), names(GATCvsGenes.List[[i]]), perl=T, value=T)
		if (exists("vecSelnames") == T) {
			vecSelnames <- c(vecSelnames, vecSel, paste(i, "expression", sep="_"))
		} else {
			vecSelnames <- c(vecSel, paste(i, "expression", sep="_"))	
		}
		if (exists("dfSel") == T) {
			dfSel <- cbind(dfSel, data.frame(vecSel=GATCvsGenes.List[[i]][[vecSel]], "expression"=GATCvsGenes.List[[i]][, ncol(GATCvsGenes.List[[i]])]))
		} else {
			dfSel <- data.frame(vecSel=GATCvsGenes.List[[i]][[vecSel]], "expression"=GATCvsGenes.List[[i]][, ncol(GATCvsGenes.List[[i]])])
		}
	}
	colnames(dfSel) <- vecSelnames
	setOfcolumns <- grep(".*target$", names(dfSel), perl=T)
	for (column in setOfcolumns) {
		if (length(grep("-1", unique(dfSel[, column]), value=T)) == 1) dfSel[dfSel[,column] == -1 & !is.na(dfSel[, column]), column] <- 0 
	}
	for (y in c(2:5)) {
		if (y == 2) {
			columns_list <- combn(setOfcolumns, y, simplify=F)
		} else {
			columns_list <- append(columns_list, combn(setOfcolumns, y, simplify=F))
		}
	}
	for (column_item_list in columns_list) {
		column_item <- unlist(column_item_list)
		if (length(column_item) == 2) {
			dfPlus <- dfSel[dfSel[, column_item[1]] == 1 & dfSel[, column_item[2]] == 1, sort(c(column_item, column_item + 1))]
			dfZero <- dfSel[dfSel[, column_item[1]] == 0 & dfSel[, column_item[2]] == 0, sort(c(column_item, column_item + 1))]
		} else if (length(column_item) == 3) {
			dfPlus <- dfSel[dfSel[, column_item[1]] == 1 & dfSel[, column_item[2]] == 1 & dfSel[, column_item[3]] == 1, sort(c(column_item, column_item + 1))]
			dfZero <- dfSel[dfSel[, column_item[1]] == 0 & dfSel[, column_item[2]] == 0 & dfSel[, column_item[3]] == 0, sort(c(column_item, column_item + 1))]
		} else if (length(column_item) == 4) {
			dfPlus <- dfSel[dfSel[, column_item[1]] == 1 & dfSel[, column_item[2]] == 1 & dfSel[, column_item[3]] == 1 & dfSel[, column_item[4]] == 1, sort(c(column_item, column_item + 1))]
			dfZero <- dfSel[dfSel[, column_item[1]] == 0 & dfSel[, column_item[2]] == 0 & dfSel[, column_item[3]] == 0 & dfSel[, column_item[4]] == 0, sort(c(column_item, column_item + 1))]
		} else if (length(column_item) == 5) {
			dfCompare <- dfSel
		} else {
			print("")
		}
		if (length(column_item) >= 2 & length(column_item) != 5) {
			dfPlus <- na.omit(dfPlus)
			dfZero <- na.omit(dfZero)
			dfCompare <- rbind(dfPlus, dfZero)
		} else {
			dfCompare <- na.omit(dfCompare)
		}
		dfCompare <- SplitDfFromColToRow(dfCompare)
		dfCompare$bound <- factor(dfCompare$bound, labels=c("Not bound", "Bound"))
		dfCompare$expression <- log2(dfCompare$expression + 1)
		dfCompare$tissue <- as.factor(dfCompare$tissue)
		boxplot_compare_folder <- "Boxplot_by_common_domains"
		protein_category_folder <- protein
		ifelse(!dir.exists(file.path(prefixDir, outputBio, outputExpr, boxplot_compare_folder)), dir.create(file.path(prefixDir, outputBio, outputExpr, boxplot_compare_folder), showWarnings=FALSE), FALSE)
		ifelse(!dir.exists(file.path(prefixDir, outputBio, outputExpr, boxplot_compare_folder, protein_category_folder)), dir.create(file.path(prefixDir, outputBio, outputExpr, boxplot_compare_folder, protein_category_folder), showWarnings=FALSE), FALSE)
		plot4 <- ggplot(aes(y=expression, x=bound, fill=bound), data=dfCompare) + geom_boxplot() + labs(title=paste("Gene expression in the", protein, "domains common for tissues", paste(levels(dfCompare$tissue), collapse=", "), sep=" ")) + facet_grid(. ~ tissue)
		pdf(file=file.path(prefixDir, outputBio, outputExpr, boxplot_compare_folder, protein_category_folder, paste(paste("Gene expression in the", protein, "domains common for tissues", paste(levels(dfCompare$tissue), collapse=" vs "), sep=" "), "pdf", sep=".")), width=14, height=14)
		print(plot4)
		dev.off()
	}
	rm(dfSel, vecSelnames, columns_list)
}

# Count housekeeping genes located in various combinations of constitutive «anti»domains
HKGenes <- read.csv(HKGenesPath, stringsAsFactors=F)
HKGenes <- GATCvsGenes[GATCvsGenes$Genes.ID %in% HKGenes$ID,]
HKGenes <- cbind(HKGenes[, c("Genes.ID", "chr", "TSS")], HKGenes[, "TSS"])
names(HKGenes) <- c("id", "chr", "start", "end")
for (comb in names(COMPARE.data.antidomains)) {
	a1 <- COMPARE.data.antidomains[[comb]][, c("seqname", "start", "end")]
	HKRanges <- with(HKGenes, IRanges(start, end))
	combRanges <- with(a1, IRanges(start, end))
	a2 <- as.data.frame(intersect(HKRanges, combRanges))
	HKGenesArea <- as.character(HKGenes[HKGenes$start %in% a2$start, "id"])
	stringName <- gsub("(?!M)[a-z\\._0-9MHT]*?\\.anti\\.domains", "", comb, perl=T)
	if (exists("HKCount") == T & exists("HKGenesFile") == T) {
		HKCount <-rbind(HKCount, data.frame("stringBetween"=stringName, "count"=length(HKGenesArea)))
		HKGenesFile[[stringName]] <- data.frame(Genes=HKGenesArea, stringsAsFactors=F)
	} else {
		HKCount <- data.frame("stringBetween"=stringName, "count"=length(HKGenesArea))
		HKGenesFile <- setNames(list(x=data.frame(Genes=HKGenesArea, stringsAsFactors=F)), stringName)
	}
}
write.table(HKCount, file=file.path(prefixDir, outputBio, "Housekeepeng_Genes_counts.csv"), sep=";", row.names=F, col.names=T, quote=F, dec=".", append=F, eol="\r\n")
save(HKGenesFile, file=file.path(prefixDir, outputBio, "Housekeepeng_Genes_list.RData"))
png(filename=file.path(prefixDir, outputBio, "HouseKeepengGenes statistics.png"), width=800, height=800)
plot(HKCount$count)
dev.off()

rm(HKCount, HKGenesFile)