#!/usr/bin/R

########################
# line 220 functions.R #
########################
#
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

######################
# line 103 domains.R #
######################
#
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

########################
# line 96 expression.R #
########################
#
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

#########################
# line 375 expression.R #
#########################
#
# # Count P-Value from protein combiations
# for (i in names(Stat_TissList)) {
# 	for (x in names(Stat_TissList[[i]])) {
# 		adf <- Stat_TissList[[i]][[x]]
# 		wilcoxStat <- wilcox.test(adf[adf$bound == "Bound", "expression"], adf[adf$bound == "Not bound", "expression"])
# 		BoundMedian <- median(adf$expression[adf$bound == "Bound"])
# 		UnboundMedian <- median(adf$expression[adf$bound == "Not bound"])
# 		if (wilcoxStat$p.value > 0.05) {
# 			print(paste(x, " in ", i, ". P-value: ", wilcoxStat$p.value, sep=""))
# 		} else if (UnboundMedian > BoundMedian) {
# 			print(paste("Normal situation in this samples:", x, "in", i, sep=" "))
# 		} else {
# 			print(paste("Binding of the protein results in enhanced expression of genes in this samples:", x, "in", i, sep=" "))
# 		}
# 	}
# }
# # Count P-Value from tissue combiations
# for (i in names(Stat_ProtList)) {
# 	for (x in names(Stat_ProtList[[i]])) {
# 		for (y in levels(Stat_ProtList[[i]][[x]]$tissue)) {
# 			bdf <- Stat_ProtList[[i]][[x]]
# 			wilcoxStat <- wilcox.test(bdf[bdf$tissue == y & bdf$bound == "Bound", "expression"], bdf[bdf$tissue == y & bdf$bound == "Not bound", "expression"])
# 			BoundMedian <- median(bdf$expression[bdf$tissue == y & bdf$bound == "Bound"])
# 			UnboundMedian <- median(bdf$expression[bdf$tissue == y & bdf$bound == "Not bound"])
# 			if (wilcoxStat$p.value > 0.05) {
# 				print(paste(i, " in combinations tissues ", x, ". Expression from ", y, ". P-value: ", wilcoxStat$p.value, sep=""))
# 			} else if (UnboundMedian > BoundMedian) {
# 				# print(paste("Normal situation in this samples:", i, "in combinations tissues", x, "Expression from", y, sep=" "))
# 			} else {
# 				print(paste("Binding of the protein results in enhanced expression of genes in this samples:", i, "in combinations tissues", x, "Expression from", y, sep=" "))
# 			}
# 		}
# 	}
# }
