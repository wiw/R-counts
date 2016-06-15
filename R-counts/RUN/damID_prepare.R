#!/usr/bin/R
print("Run script")

# Make samples list file
##########################
MakeSamplesListFile(sourceDir, damIdLocation)

# Load GATC counts in data frame
################################
# Проверка на предыдущие запуски программы. Если не было, то тогда gatc файл создается заново
if (startCol == 0) {
	step01 <- read.delim(alreadyRun, header=T, as.is=T, dec=".")
	startCol <- ncol(step01)
	gatcs <- step01
} else {
	gatcs <- read.delim(gatcFile, header=T, as.is=T, dec=".")
}
# Из функции MakeSamplesListFile. Если есть edge/inner риды, тогда одна обработка, если нет - другая 
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
# Нужно ли использовать не весь набор данных в файле с образцами
if (needSomeFiles == T) {
	samplesList <- samplesList[useSomeFiles, ]
}
# Объединение координат gatcs фрагментов и данных DamID
gatcs <- cbind(gatcs, matrix(data=NA, nrow=nrow(gatcs), ncol=nrow(samplesList)))
for (i in 1:nrow(samplesList)){
	colnames(gatcs)[startCol+i] <- samplesList$id[i]
	load(file=samplesList$path[i])
	if (all(gatcs$ID.il == reads2GATC$ID)) gatcs[, startCol + i] <- reads2GATC$count
}
rm(i)
# Если присутствуют edge/inner риды, тогда переменные gatcs и samplesList перезаписываются в соответствии с этими условиями
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
# Фильтрация данных DamID. Производится поиск максимальной корреляции между репликами при различном размере выборки. А затем отбрасыываются те данные которые не попали в выборку
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
DATAs.ave <- DATAs
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

# Averaging Replicates before Normalization
###########################################
		DATAs.ave[[name]] <- DATAs.ave[[name]][, -c(8:ncol(DATAs.ave[[name]]))]
		# listAve <- samplesList[!(samplesList$protein == "DAM"), c(1:5)]  # remove row's with DAM
		listAve <- samplesList[, c(1:5)]
		listAve$average <- paste(listAve$tissue, listAve$protein, listAve$conditions, sep=".")
		uniqueAveSamples <- unique(listAve$average)
		for (item in uniqueAveSamples) {
			item.ave <- paste(item, ".ave", sep="")
				vector <- grep(paste(item, "\\.[0-9]", sep=""), names(DATAs[[name]]))
				weHaveOneReplicate <- length(vector)
				if (weHaveOneReplicate > 1) {
					DATAs.ave[[name]][[item.ave]] <- rowMeans(DATAs[[name]][, vector])
				} else {
					print(paste("This", item, "have one replicate!", sep=" "))
						DATAs.ave[[name]][[item.ave]] <- DATAs[[name]][, vector]
				}	
		}
		dam.ave <- assign(paste("dam.ave", name, sep="."), paste("DF_Counts_Step_9.1_Averaged_", name, "_", currentDate, ".csv", sep=""))
		WriteIntermediateFiles(source=DATAs.ave[[name]], output.file=dam.ave)

# DAM Normalization
###################
		print("DAM Normalization")
		DATAs.norm[[name]] <- DATAs.norm[[name]][, -c(8:ncol(DATAs.norm[[name]]))]
		listNorm <- samplesList[,1:5]
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

# Averaging Replicates after Normalization
##########################################
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