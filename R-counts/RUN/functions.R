#!/usr/bin/R
# Make samples list file
# Функция которая переформатирует файл соответствий файлов fastq их коротким названиям в файл для обработки в R
# Как использовать:
# MakeSamplesListFile(SOURCE, DAMID)
# SOURCE - это папка где находятся ваши .RData файлы
# DAMID - путь до файла, который должен иметь следующий формат
# > Data.set        fastq.file
# > [tissue].[protein].[conditions].[replicate]     [file_name].fastq.gz

# Результат работы будет сохранен в переменную samplesList,
# которая имеет следующие заголовки
# > id tissue protein conditions replicate path
# А кроме этого в рабочей директории будет создан файл rdata_description.csv, с содержимым переменной samplesList

# Функция способна корректно обрабатывать .RData файлы содержащие в имени файла следующие строки: edge, inner, _edge, _inner, paired
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

# Calculate coordinate for filter
# Функция для расчетов координат при использовании фильтрации данных
# Как использовать:
# X.Coord <- sapply(your_dataframe, CalculateCoordinate, y=1)
# Y.Coord <- sapply(your_dataframe, CalculateCoordinate, y=2)
# your_dataframe - dataframe с координатами GATC фрагментов, который содержит 7 столбцов от GATCs.txt и не менее 1 столбца с данными DamID
#################################
CalculateCoordinate <- function(x, y){
  if (index[x] == T){
    RhsMatrix <- matrix(c(-1*Intercept,-(DATA.filter[x, 9]+DATA.filter[x, 8]/Slope)), nrow=2, ncol=1, byrow=T)
    Result <- Coef.Matrix %*% RhsMatrix
    Result[y,1]
  }
}

# Write intermediate files function
# Простая функция, которая записывает на диск содержимое переменной
# Как использовать:
# WriteIntermediateFiles(source, output.file)
# source - переменная, которую необходимо записать на диск, лучше использовать dataframe
# output.file - название файла с расширением (обычно это csv или txt), который будет сохранен в подпапку в вашей рабочей директории. Имя подпапки определяется перменной prefixDir
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
# Мини-функции, которые используются в других функциях: MainCorrelations, HeatmapBySelection.
# Необходимы для построения корреляции Pearson&Spearman "всех со всеми"
# Как использовать:
# PearsonAndSpearmanCorrelations(dataSet, use.method, use.opt)
# dataSet - dataframe с координатами GATC фрагментов, который содержит 7 столбцов от GATCs.txt и не менее 1 столбца с данными DamID
# use.method - два возможных варианта - "pearson", "spearman"
# use.opt - смотри справку по параметру use функции cor
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
# PearsonAndSpearmanCorrelationsHeatmapMod(dataSet1, dataSet2, use.method, use.opt)
# dataSet1, dataSet2 - два разных dataframe со свойствами аналогичным dataSet
# остальные переменные аналогичны переменным из функции PearsonAndSpearmanCorrelations

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
# Функция в тестовом режиме, необходима была для объединения произвольных пар столбцов в один. На вход принимает переменную содержащую таблицу для такой обработки. Сейчас функция не применяется.
##########################
CombineSamples <- function(dataFrame) {
	with(dataFrame, data.frame(DAM.1=`DAM-1.FCC4JPEACXX_L6_R1` + `DAM-1.FCC4JPEACXX_L6_R2`, DAM.2=`DAM-2.FCC4JPEACXX_L6_R1` + `DAM-2.FCC4JPEACXX_L6_R2`, LAM.1=`LAM-1.FCC4JPEACXX_L6_R1` + `LAM-1.FCC4JPEACXX_L6_R2`, LAM.2=`LAM-2.FCC4JPEACXX_L6_R1` + `LAM-2.FCC4JPEACXX_L6_R2`))
}

# Main correlations function
# Функция, которая считает корреляцию Спирмана и Пирсона в наборе данных и выводит эту информацию в виде текстового файла и pdf с тепловой картой в папку prefixDir
# Как использовать:
# MainCorrelations(dataSet, corrMethod, labelHeatmap, use.opt="everything", suffixCSV, suffixPDF, corr.on.file, createPDF=T, counts=F)
# dataSet - dataframe с координатами GATC фрагментов, который содержит 7 столбцов от GATCs.txt и не менее 1 столбца с данными DamID
# corrMethod - текстовый вектор, содержащий варианты корреляций. Смотри справку по параметру method функции cor
# labelHeatmap - вектор, содержащий информацию о том нужно ли включать плотность в тепловую карту или нет. "B" - включает информацию о плотности и убирает легенду, любой другой ключ делает всё наоборот
# use.opt - смотри справку по параметру use функции cor
# suffixCSV - Произвольная метка для имени csv файла
# suffixPDF - Произвольная метка для имени pdf файла, работает только если параметр createPDF = TRUE
# corr.on.file - имя csv файла откуда были взяты данные для расчета корреляции
# createPDF - создавать или нет pdf с тепловой картой
# counts - логическая переменная, должна быть TRUE, когда suffixCSV содержит "Counts"
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
# Функция для расчета автокорреляции и плотности распределения GATC фрагментов. Результат работы pdf файл в папке prefixDir
# Как использовать:
# AcfOnData(dataSet, labelAcf, method, suffixPDF, ylab.val, na.data)
# dataSet - dataframe с координатами GATC фрагментов, который содержит 7 столбцов от GATCs.txt и не менее 1 столбца с данными DamID
# labelAcf - метка "A_ALL" которая указывает на использование всего набора данных, любая другая метка будет использовать положительный набор данных по микрочипам
# method - метод расчета, автокорреляция "acf" или плотность "density"
# suffixPDF - Произвольная метка для имени pdf файла
# ylab.val - подпись оси y
# na.data - если TRUE, тогда значения NA будут удалены из dataSet
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
# Функция для представления обработанных данных в формате wig или gff, файлы сохраняются в папке outputWig или outputGff
# Как использовать:
# DamIdSeqToWigGff(dataSet)
# dataSet - dataframe с координатами GATC фрагментов, который содержит 7 столбцов от GATCs.txt и не менее 1 столбца с данными DamID
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
# Функция для создания визуализации распределения данных, bmp файлы сохраняются в папке outputScttr
# Как использовать:
# ScatterPlotting3D(dataSet, tag)
# dataSet - dataframe с координатами GATC фрагментов, который содержит 7 столбцов от GATCs.txt и не менее 1 столбца с данными DamID
# tag - произвольный текстовый вектор длиной 1
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
# Функция, которая формирует домены объединенные по признаку feature.type
# Как использовать:
# you_variable <- FeatureCalls.to.GFF.like(start.coordinate, end.coordinate, feature.type)
# start.coordinate - столбец в таблице с точками начала района
# end.coordinate - столбец в таблице с точками конца района
# feature.type - столбец в таблице по которому надо производить объединение данных
##############################
FeatureCalls.to.GFF.like <- function(start.coordinate, end.coordinate, feature.type) {
   e.ind <- cumsum(rle(feature.type)$lengths)
   s.ind <- c(1, e.ind[1:(length(e.ind)-1)]+1)
   gff <- data.frame(start=start.coordinate[s.ind], end=end.coordinate[e.ind], value=feature.type[s.ind])
   invisible(gff)
}

# Use fit third-clustering model function
# Это модифицированная функция из пакета snapCGH,
# с установленным ограничением на поиск трех возможных
# cостояний по скрытой марковской модели. Используется
# в функции runBioHMM
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

# Run biological Headen Mark Model function
# Это модифицированная функция из пакета snapCGH.
# Как использовать:
# you_variable <- runBioHMM(mval, datainfo)
# mval - столбец из datainfo, содержащий значения для обработки HMM
# datainfo - dataframe с координатами GATC фрагментов, который содержит 7 столбцов от GATCs.txt и не менее 1 столбца с данными DamID
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


# Make Venn diagram
###################
pcentFun <- function(x, y) {
	100 * (x / y)
}

# !!! Before use this function, please perform this instructions: http://stackoverflow.com/a/15315369/4424721
# Функция создающая pdf файл, содержащий Венн диаграмму, зависящую от веса значений и выраженных в абсолютных и относительных значениях
# Как использовать:
# MakeVennDiagram(x, set, v1, v2)
# x - dataframe с координатами GATC фрагментов, который содержит 7 столбцов от GATCs.txt и не менее 1 столбца с данными DamID
# set - регулярное выражение, описывающее имена столбцов из x используемые для анализа
# v1, v2 - текстовый вектор содержащий названия тканей и/или белков, чьи комбинации необходимо проанализировать
# Пример:
# MakeVennDiagram(x=DATA.venn, set=paste(t.s,".*bound", sep=""), v1=t.s, v2="LAM, HP1, PC")
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
# Функция создает тепловую карту из различных сочетаний столбцов в используемой таблице, файлы pdf создаются в папке prefixDir/outputHeatmap 
# Как использовать:
# HeatmapBySelection(dataSet, corrMethod, use.opt)
# dataSet - dataframe с координатами GATC фрагментов, который содержит 7 столбцов от GATCs.txt и не менее 1 столбца с данными DamID
# corrMethod - смотри справку к параметру method функции cor
# use.opt - смотри справку к параметру use функции cor
# Пример:
# HeatmapBySelection(dataSet=DATA.part, corrMethod=corrMethod, use.opt="pairwise.complete.obs")
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
		bio.set <- list(TB = tissue_bio_set, PB = protein_bio_set)
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
# Функция, которая размещает много различных графиков на одной или нескольких страницах одного файла. Используется вместе с пакетом ggplot2 и функцией ggplot
# Как использовать:
# print(multiplot(..., plotlist, cols, layout))
# ... - перечисление переменных ggplot
# plotlist - переменные ggplot в виде списка
# cols - укажите необходимое число колонок в файле
# layout - это матрица, указывающая расположение переменных на странице. Расчитывается автоматически, но можно и указать явно как расположить элементы
# Пример:
# print(multiplot(plotlist=tid2_list1, cols=3, layout=matrix(c(1:3), ncol = 3, nrow = 1, byrow=T)))
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

# Statistical functions about domains size
# Статистическая функция, создает csv файл в папке prefixDir/outputBio с отчетом о суммарном размере доменов белков в различных тканях в долях от полного генома
# Как использовать:
# MakeReportDomainsSize(inputData, dscr, includeHet=T)
# inputData - объект типа list, содержащий таблицы с заголовками:
#  > seqname	source	feature	start	end	score	strand	frame	attribute
# dscr - отметка в имени файла, опционально
# includeHet - если TRUE, тогда для анализа будут использованы Het хромосомы
# Пример:
# MakeReportDomainsSize(DOMAIN.data.filt, "filtered", includeHet=F)
##########################################
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
		partOfDomainFromGenome <- round(pcentFun(domainsSize, AllGenome), digits=3)
		domainSizeReport[grep(paste("^",item, "$", sep=""), names(inputData)),] <- c(item, AllGenome, domainsSize, partOfDomainFromGenome, fifthCol)
	}
	reportFile <- file.path(prefixDir, outputBio, paste("Report_about", dscr, modificator,"domains_part_from_genome_D.melanogaster.csv", sep="_"))
	write.table(domainSizeReport, file=reportFile, sep=";", row.names=F, col.names=T, quote=F, dec=".", append=F)
}


# Intersect domains between them
# Функция, которая находит общие домены и антидомены в различных комбинациях доменов и антидоменов белков и тканей.
# В результате работы функции создаются gff файлы в папке prefixDir/outputGff/outputDomain/dscr
# Как использовать:
# IntersectDomain(inputData, Tset, Pset, dscr, ScoreValue)
# inputData - объект типа list, содержащий таблицы с заголовками:
#  > seqname	source	feature	start	end	score	strand	frame	attribute
# Tset - текстовый вектор с наборами тканей используемых для анализа
# Pset - текстовый вектор с наборами тканей используемых для анализа
# dscr - маркировка набора данных
# ScoreValue - принимает значения -1 или 1, для антидоменов и доменов соответственно
# Пример:
# IntersectDomain(inputData=DOMAIN.data, Tset=tissue_bio_set, Pset=protein_bio_set, dscr="constitutive", ScoreValue=1)
################################
IntersectDomain <- function(inputData, Tset, Pset, dscr, ScoreValue) {
	DOMAIN.data.part <- inputData[grep(paste("(", paste(tissue_bio_set, collapse="|"), ")\\.(", paste(protein_bio_set, collapse="|"), ")\\.(", paste(conditions_bio_set, collapse="|"), ")\\..*", sep="") ,names(inputData), perl=T)]
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

# Split dataframe from column to row
# Функция, которая преобразует таблицу - меняет местами колонки со строками
# Как использовать:
# SplitDfFromColToRow(data)
# data - объект типа dataframe
####################################
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

# Assign GATCs damID values to Genes position
# Функция, которая устанавливает соответствие между координатами GATC фрагментов и положением генов
# Как использовать:
# your_variable <- AssignGATCsToGenes(Genes, GATCs, und=F, tissue.damid="", use.und.domain=F)
# Genes - объект типа dataframe с заголовками 
# > ID chr start end ratio
# GATCs - объект типа dataframe с заголовками 
# > ID chr start end any_GATC_values
# und - когда TRUE для анализа используются данные по недорепликации
# tissue.damid - только когда use.und.domain = FALSE, метка ткани, которую нужно использовать для анализа
# use.und.domain - использование данные по недорепликации в доменах TRUE, или "сырых" данных по недорепликации FALSE - по умолчанию
# Пример:
# und_vs_gatc <- AssignGATCsToGenes(UndData.list[[i]], gatcUnd, und=T, tissue.damid="FB")
#############################################
AssignGATCsToGenes <- function(Genes, GATCs, und=F, tissue.damid="", use.und.domain=F) {
	gr_genes <- makeGRangesFromDataFrame(Genes, keep.extra.columns=T)
	gr_gatcs <- makeGRangesFromDataFrame(GATCs, keep.extra.columns=T)
	Results <- as.data.frame(findOverlaps(gr_genes, gr_gatcs))
	if (und == F) {
		# For expression data
		CompareResults <- cbind(GATCs[Results$subjectHits, c("ID","chr","start","end")], data.frame("Genes.ID" = Genes$ID[Results$queryHits], "TSS" = Genes$start[Results$queryHits], "strand" = Genes$strand[Results$queryHits]))
	} else {
		if (use.und.domain == F) {
			# For underreplication data
			if (nchar(tissue.damid) > 0) gatcNames <- c("ID","chr","start","end",grep(paste(tissue.damid, ".*", sep=""), names(GATCs), perl=T, value=T)) else gatcNames <- c("ID","chr","start","end")
			CompareResults <- cbind(GATCs[Results$subjectHits, gatcNames], Genes[Results$queryHits, "ratio"])
			names(CompareResults)[ncol(CompareResults)] <- "ratio"
			} else {
				# For underreplication domain data
				CompareResults <- GATCs
				CompareResults$fb_und_domain <- 1
				CompareResults$fb_und_domain[Results$subjectHits] <- 0
			}
	}
	CompareResults <- CompareResults[with(CompareResults, order(chr, ID)), ]
	rownames(CompareResults) <- NULL
	return(CompareResults)
}

# Count Occurences :)
# Функция, которая считает количество вхождений набора символов в строке
# Как использовать:
# CountOccurrencesInCharacter(char, s)
# char - искомая строка или символ
# s - текстовая строка в которой производится поиск
#####################
CountOccurrencesInCharacter <- function(char, s) {
    s2 <- gsub(char,"",as.character(s))
    return ((nchar(as.character(s)) - nchar(s2))/nchar(char))
}