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

# Declare variables
###################
	prefixDir <- "RUN22-06-2015" # output directory into working directory
	onlyEdge <- F # use only edge reads to counts or not
	workDir <- getwd()	# working directory (WD)
	sourceDir <- "/home/anton/backup/output/RUN22-06-2015" # location your RData files. You can specify the highest folder as it is possible. Searching runs recursively.
	damIdLocation <- "/home/anton/data/DAM/RUN/damid_description.csv" # location your DamID-Description file
	outputGff <- "gff"	# output folder for gff in WD
	outputWig <- "wig"	# output folder for wig in WD
	outputScttr <- "scatter_plots"	# output folder for scatter plots in WD
	outputDomain <- "domains"
	outputCleanStat <- "clean_stat"
	outputBio <- "Bio"
	outputHeatmap <- "Heatmap"
	startCol <- 7	# the number of last column in GATCs file, default "7"
	gatcFile <- paste(workDir, "GATCs_mod.txt", sep="/")	# location you GATCs file
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
	tissue.bio.set <- c("Kc167", "BR", "FB", "NRN", "Glia") # Part of tissies from all
	protein.bio.set <- c("LAM", "HP1", "PC") # Part of protein from all
	conditions.bio.set <- c("m", "m_25mkM4HT", "mf_min")

# Create folders
################
	dir.create(file.path(workDir, prefixDir), showWarnings = FALSE)
	dir.create(file.path(workDir, prefixDir, outputGff), showWarnings = FALSE)
	dir.create(file.path(workDir, prefixDir, outputWig), showWarnings = FALSE)
	dir.create(file.path(workDir, prefixDir, outputScttr), showWarnings = FALSE)
	dir.create(file.path(workDir, prefixDir, outputGff, outputDomain), showWarnings = FALSE)
	dir.create(file.path(workDir, prefixDir, outputCleanStat), showWarnings = FALSE)
	dir.create(file.path(workDir, prefixDir, outputHeatmap), showWarnings = FALSE)

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

# Use fit model function
########################
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
        temp2 <- clara(data, 2)
        init.mean.two <- temp2$medoids
        init.var.two <- vector()
        if (var.fixed == FALSE) {
            for (i in 1:2) {
                if (length(temp2$data[temp2$clustering == i]) > 
                  1) 
                  init.var.two[i] <- log(sqrt(var(temp2$data[temp2$clustering == 
                    i])))
                else init.var.two[i] <- log(0.5)
            }
        } else {
            init.var.two[1:2] <- log(sqrt(var(data)))
        }
        z2.init <- c(init.mean.two[, 1], init.var.two, -1, -3.6, 
            -3.6, 0)
        z.pre <- run.nelder(numobs, z2.init, data, 
            covars, var.fixed, epsilon, numiter, i)
        if (!is.nan(z.pre$x[1])) {
            z2 <- find.param.two(z.pre, var.fixed)
        } else {
            z2 <- NULL
        }
        if (aic) {
            factor <- 2
        } else if (bic) {
            factor <- log(numobs) * delta
        } else {
            stop("No criteria selected")
        }
        z <- z2
        nstates <- 2
        trans.mat <- list()
        for (j in 1:(length(data) - 1)) {
            trans.mat[[j]] <- z$LH.trans + exp(-(covars[j, 
              1]^(z$rate1)) * prod(covars[j, -1])) * z$RH.trans
        }
        Vit.seg <- Viterbi.two(data, 
            z, trans.mat)
        maxstate.unique <- unique(Vit.seg)
        mean <- rep(0, length(data))
        var <- rep(0, length(data))
        for (m in 1:length(maxstate.unique)) {
            mean[Vit.seg == maxstate.unique[m]] <- mean(data[Vit.seg == 
              maxstate.unique[m]])
            var[Vit.seg == maxstate.unique[m]] <- var(data[Vit.seg == 
              maxstate.unique[m]])
        }
        out <- cbind(matrix(Vit.seg, ncol = 1), matrix(mean, 
            ncol = 1), matrix(var, ncol = 1))
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
if (onlyEdge == T) {
	samplesList <- samplesList[grep("edge" ,samplesList$id), ]
} else {
	modS <- samplesList[grep("edge", samplesList$id), ]
	modS$id <- gsub("(.+)(edge)(.+)", paste("\\1", "all", "\\3", sep=""), modS$id, perl=T)
	modS$conditions <- gsub("(.+)edge", paste("\\1", "all", sep=""), modS$conditions, perl=T)
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
 		
 		if (nrow(classifier) == 2){
			# notify if there is a problem with groupping the BioHMM outputs ("1s" and "2s") 
			if (classifier$Group.1[1] !=1){
				# print(paste(DATA.name, ", ", dataframe.name, " - 1st classifier is not 1!"), sep="")
			}
			if (classifier$Group.1[2] !=2){
				# print(paste(DATA.name, ", ", dataframe.name, " - 1st classifier is not 2!"), sep="")
			}
			 
			# if "1s" are targets and "2s" are non-targets
			if ((classifier$x[1] > 0) & (classifier$x[2] < 0)) {
			# print(paste(DATA.name, ", ", dataframe.name, " - '1s' are targets and '2s' are non-targets", sep=""))
			HMM.data[[DATA.name]][[dataframe.name]]$target[(HMM.data[[DATA.name]][[dataframe.name]]$BioHMM.output == 1)] <- 1
			HMM.data[[DATA.name]][[dataframe.name]]$target[(HMM.data[[DATA.name]][[dataframe.name]]$BioHMM.output == 2)] <- 0
			}
			# if "1s" are non-targets and "2s" are targets
			if ((classifier$x[1] < 0) & (classifier$x[2] > 0)) {
			# print(paste(DATA.name, ", ", dataframe.name, " - '1s' are non-targets and '2s' are targets", sep=""))
			HMM.data[[DATA.name]][[dataframe.name]]$target[(HMM.data[[DATA.name]][[dataframe.name]]$BioHMM.output == 2)] <- 1
			HMM.data[[DATA.name]][[dataframe.name]]$target[(HMM.data[[DATA.name]][[dataframe.name]]$BioHMM.output == 1)] <- 0
			}

			# if it is not clear what is what
			if ((classifier$x[1] < 0) & (classifier$x[2] < 0)) {
				# print(paste(DATA.name, ", ", dataframe.name, " - all data less then zero!", sep=""))
			}

			# if it is not clear what is what
			if ((classifier$x[1] > 0) & (classifier$x[2] > 0)) {
				# print(paste(DATA.name, ", ", dataframe.name, " - all data more then zero!", sep=""))
			}
		} else {
			print(paste(DATA.name, ", ", dataframe.name, " - not enough data!", sep=""))
		}
		if ("target" %in% names(HMM.data[[DATA.name]][[dataframe.name]])){
			domains <- FeatureCalls.to.GFF.like(start.coordinate=HMM.data[[DATA.name]][[dataframe.name]]$start, end.coordinate=HMM.data[[DATA.name]][[dataframe.name]]$end, feature.type=HMM.data[[DATA.name]][[dataframe.name]]$target)
			domains <- cbind(chr=HMM.data[[DATA.name]][[dataframe.name]]$chr[1], domains, stringsAsFactors=F)
			domains <- domains[domains$value==1,]
			domains <- cbind(domains, NA, NA, NA, NA, NA)
			domains <- domains[,c(1,5,6,2,3,4,7,8,9)]
			names(domains) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
			protein.name <- tolower(sub("([a-zA-Z]+)\\.([a-zA-Z0-9]+)\\.([a-zA-Z_0-9]+)", "\\2.domain", DATA.name, perl=T))
			DATA.domain.name <- sub("(.*)", "\\1.domains", DATA.name, perl=T)
			domains$feature[domains$score==1] <- protein.name
			domains[, c(2,7:9)] <- "."
			if (chr == 1) {DOMAIN.data[[DATA.domain.name]] <- domains;} else {DOMAIN.data[[DATA.domain.name]] <- rbind(DOMAIN.data[[DATA.domain.name]], domains);}
		} else {
			print(paste("No data in ", DATA.name, " - ", dataframe.name, ".", sep=""))
		}
	}
	gff.file <- file.path(prefixDir, outputGff, outputDomain, paste(DATA.domain.name, "gff", sep="."))
	write.table(DOMAIN.data[[DATA.domain.name]], file=gff.file, sep="\t", row.names=F, col.names=F, quote=F, dec=".", append=F)
}
print("Congratulations!!!")

print("Search Biology meaning...")
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
	compare.vector <- match(DATA.venn$GATC.ID, target.df$ID)
	DATA.venn[[name.bound]] <- target.df$target[compare.vector]
	rm(target.df)
}
print("Make Venn Diagram")
for (t.s in tissue.bio.set)	MakeVennDiagram(x=DATA.venn, set=paste(t.s,".*bound", sep=""), v1=t.s, v2="LAM, HP1, PC")
for (t.s in tissue.bio.set)	MakeVennDiagram(x=DATA.venn, set=paste(t.s,".*\\.(LAM|HP1).*bound", sep=""), v1=t.s, v2="LAM, HP1")
for (p.s in protein.bio.set) MakeVennDiagram(x=DATA.venn, set=paste(".*\\.", p.s, ".*bound", sep=""), v1=p.s, v2="Kc167, BR, FB, NRN, Glia")
for (p.s in protein.bio.set) MakeVennDiagram(x=DATA.venn, set=paste("(Kc167|BR|FB)\\.", p.s, ".*bound", sep=""), v1=p.s, v2="Kc167, BR, FB")

# Count area domain from whole genome
#####################################
domainSizeReport <- as.data.frame(matrix(NA, nrow=0, ncol=5))
names(domainSizeReport) <- c("Item.name", "Genome.size", "Domain.size", "Ratio", "Selected.chromosomes")
for (item in names(DOMAIN.data)) {
	getChrList <- paste("chr", unique(DOMAIN.data[[item]][,1]), sep="")
	domainsSize <- sum(DOMAIN.data[[item]]$end-DOMAIN.data[[item]]$start)
	chr.lengths <- integer()
	for (i in getChrList) chr.lengths <- append(chr.lengths, length(DNAString(Dmelanogaster[[i]])))
	AllGenome <- sum(chr.lengths)
	partOfDomainFromGenome <- paste(round(pcentFun(domainsSize, AllGenome), digits=3), "%", sep="")
	domainSizeReport[grep(item, names(DOMAIN.data)),] <- c(item, AllGenome, domainsSize, partOfDomainFromGenome, paste(unique(DOMAIN.data[[item]][,1]), collapse=", "))
}
reportFile <- file.path(prefixDir, outputBio, "Report_about_domains_part_from_genome_D.melanogaster.csv")
write.table(domainSizeReport, file=reportFile, sep=";", row.names=F, col.names=T, quote=F, dec=".", append=F)