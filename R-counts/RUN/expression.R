#!/usr/bin/R
#######################################################################
####################### GATCs overlaps Genes ##########################
#######################################################################

# Intersect with promoters only
# Подготовительные процедуры, поиск пересечений между координатами GATC фрагментов и как промоторами генов так и полным телом гена
###############################
gatcG <- read.delim(gatcFile4Genes, header=T, sep="\t", stringsAsFactors=F)
Genes <- read.delim(genesFilePath, header=T, as.is=T, dec=".")
Genes <- Genes[Genes$chr %in% intersect(unique(Genes$chr), unique(gatcG$chr)),]
gatcG <- gatcG[gatcG$chr %in% intersect(unique(Genes$chr), unique(gatcG$chr)),]
Genes[Genes$strand == "+", "end"] <- Genes[Genes$strand == "+", "start"]
Genes[Genes$strand == "-", "start"] <- Genes[Genes$strand == "-", "end"] # Gene position like as point

# Таблица соответствия между GATC.ID (ID) и Genes.ID, имеет следующую структуру:
#              ID chr start   end    Genes.ID   TSS strand
# 1 r5GATC2L00040  2L  6923  7693 FBgn0031208  7529      +
# 2 r5GATC2L00081  2L 21359 21643 FBgn0002121 21372      -
# 3 r5GATC2L00091  2L 25086 26562 FBgn0031209 25151      -
# Пересечение происходит между протяженными GATC фрагментами и точечными координатами промоторов
GATCvsGenes <- AssignGATCsToGenes(Genes, gatcG)

# Установка соответствия между именем GATC фрагмента в GATCvsGenes и его положением в gatcG
idx.gatc.ori <- unlist(sapply(GATCvsGenes$ID, function(x) grep(x, gatcG$ID)))
GATCvsGenes.List <- list()

# Intersect with gene body
##########################
GenesOrig <- read.delim(genesFilePath, header=T, as.is=T, dec=".") # Gene occupies a certain domain

# Поиск GATC фрагментов которые полностью или частично перекрывыаются с телом гена
# таблица имеет следующую структуру:
#              ID chr start   end    Genes.ID  TSS strand
# 1 r5GATC2L00040  2L  6923  7693 FBgn0031208 7529      +
# 2 r5GATC2L00041  2L  7694  7716 FBgn0031208 7529      +
# 3 r5GATC2L00042  2L  7717  8552 FBgn0031208 7529      +
# 4 r5GATC2L00043  2L  8553  8796 FBgn0031208 7529      +
# 5 r5GATC2L00044  2L  8797  9694 FBgn0031208 7529      +
# 6 r5GATC2L00047  2L  9803 12441 FBgn0002121 9836      -
# Пересечение происходит между протяженными GATC фрагментами и протяжнными координатами тела гена.
# Отдельно указывается сайт инициации транскрипции (TSS)
GATCvsGenesOrig <- AssignGATCsToGenes(GenesOrig, gatcG)

# Установка соответствия между именем GATC фрагмента в GATCvsGenes и его положением в gatcG
idx.gatc.ori.body <- unlist(sapply(GATCvsGenesOrig$ID, function(x) grep(x, gatcG$ID)))
GATCvsGenesOrig.List <- list()

# Формирование списка dfExpression.List, состоящего из объектов типа dataframe,
# в котором указаны нормированные данные по экспрессии генов в соответствующих тканях
#            id   expression
# 1 FBgn0000003     0.000000
# 2 FBgn0000008   518.270644
# 3 FBgn0000014     8.965697
# 4 FBgn0000015     1.159844
#########################
dfExpression.List <- list()
for (i in tissue_bio_set) {
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
		exprData <- newCountDataSet(countData=rnaSeq, conditions=conditions)
		exprData <- estimateSizeFactors(exprData)
		exprData <- estimateDispersions(exprData)
		exprResults <- nbinomTest(exprData, i, i)
	} else {
		exprResults <- read.delim(tissueSelection, header=F, stringsAsFactors=F, sep="\t")
		names(exprResults) <- c("id", "baseMean")
	}
	dfExpression.List[[i]] <- data.frame(id=exprResults$id, expression=exprResults$baseMean)

	# Таблицы GATCvsGenes* записываются в отдельные переменные tmpDF*,
	# поскольку далее эти переменные будет подвергаться дополнительному форматированию
	tmpDF <- GATCvsGenes
	tmpDF.Orig <- GATCvsGenesOrig
	if (i == "Kc167") {
		tempDATA <- DATA.part
	} else {
		tempDATA <- cbind(DATAs.norm.ave$DATA, DATAs.ave$DATA[,c(8:ncol(DATAs.ave$DATA))])
	}

	# 	Конечный результат этого кода форматирование списка, состоящего из объектов типа dataframe такого вида:
	# 	               ID chr start   end    Genes.ID   TSS strand    BR.LAM.m_all.norm.ave BR.LAM.m_all.ave BR.LAM.m_target BR.LAM.m_target.filt    expression
	# 40  r5GATC2L00040  2L  6923  7693 FBgn0031208  7529      +          -3.1618273              9.5              -1                   -1			  23.92050
	# 81  r5GATC2L00081  2L 21359 21643 FBgn0002121 21372      -                  NA              0.0              NA                   NA			5471.38956
	# 90  r5GATC2L00091  2L 25086 26562 FBgn0031209 25151      -          -1.4555819              6.5              -1                   -1			  25.49604
	# 172 r5GATC2L00173  2L 58765 59732 FBgn0051973 59242      -           0.3967371             33.5               0                    0			 697.05598
	# 187 r5GATC2L00188  2L 66423 67274 FBgn0067779 67044      +          -3.6513258             19.0              -1                   -1			1910.33491
	# 204 r5GATC2L00206  2L 71991 72652 FBgn0031213 72388      +          -3.4053898             10.5              -1                   -1			1336.69455
	# Эта таблица содержит информацию о положении GATC фрагментов, какие промоторы генов с ними пересекаются, соответствующие значения DamID (нормированные усредненные, только усредненные),
	# связан ли данный фрагмент с исследуемым белком и данные по экспрессии гена в конкретном районе.
	for (z in protein_bio_set) {
		dataSet <- grep(paste(i, "\\.", z, "\\.(", paste(conditions_bio_set, collapse="|"), ")_(all|edge).*",sep=""), names(tempDATA), value=T, perl=T)
		df <- tempDATA[, c("ID", dataSet)]
		if (i == "Kc167") {
			hmmSet <- grep(paste(i, "\\.", z, "\\.bound", sep=""), names(DATA.venn), value=T, perl=T)
			hmm.df <- DATA.venn[, c(1:7, grep(hmmSet, names(DATA.venn)))]
			hmm.df <- hmm.df[!(is.na(hmm.df[,ncol(hmm.df)])),]
			names(hmm.df)[ncol(hmm.df)] <- "target"
		} else {
			hmmSet <- grep(paste(i, "\\.", z, "\\.(", paste(conditions_bio_set, collapse="|"), ")$",sep=""), names(HMM.data), value=T, perl=T)
			hmm.df <- rbind.fill(HMM.data[[hmmSet]])
		}
		idx.hmm <- which(df$ID %in% hmm.df$ID)
		df <- cbind(df, as.data.frame(matrix(NA, ncol=2, nrow=nrow(df))))
		names(df)[(ncol(df)-1):ncol(df)] <- paste(hmmSet, c("target", "target.filt"), sep="_")
		df[,(ncol(df)-1)][idx.hmm] <- hmm.df$target
		if (i != "Kc167") {
			df[,ncol(df)][idx.hmm] <- hmm.df$target.filt
		}
		tmpDF <- cbind(tmpDF, df[idx.gatc.ori, c(2:ncol(df))])
		tmpDF.Orig <- cbind(tmpDF.Orig, df[idx.gatc.ori.body, 2])
		names(tmpDF.Orig)[ncol(tmpDF.Orig)] <- names(df)[2]
	}
	idx.fbgn.ori <- match(tmpDF$Genes.ID, exprResults$id)
	if (i != "Kc167") {
		damCol <- grep(paste(i, "\\.DAM.*", sep=""), names(tempDATA), value=T)
		tmpDF$DAM <- tempDATA[idx.gatc.ori, damCol]
		names(tmpDF)[which(names(tmpDF) == "DAM")] <- damCol
	}
	tmpDF$expression <- exprResults$baseMean[idx.fbgn.ori]
	GATCvsGenes.List[[i]] <- tmpDF

	tmpDF.Orig <- tmpDF.Orig[, -c(1:4, 6, 7)]

	# Для случая пересечения GATC фрагментов с телом гена команда ниже усредняет значения DamID в колонках по строкам
	tmpDF.Orig <- ddply(tmpDF.Orig, "Genes.ID", numcolwise(mean, na.rm = TRUE)) # Averaging the DamID values in a column by row
	for (colItem in c(2:ncol(tmpDF.Orig))) tmpDF.Orig[is.nan(tmpDF.Orig[,colItem]) == T, colItem] <- NA
	idx.fbgn.ori.body <- match(tmpDF.Orig$Genes.ID, exprResults$id)
	tmpDF.Orig$expression <- exprResults$baseMean[idx.fbgn.ori.body]
	GATCvsGenesOrig.List[[i]] <- tmpDF.Orig

	rm(rnaSeq, tmpDF, tmpDF.Orig)
}

# Count statistic
# На основе списка GATCvsGenes.List, полученного ранее, формируется графическая статистика зависимости экспрессии генов от связывания белка
# Статистика представлена в виде боксплота (2-х и 3-х факторные состояния) и гистограммы (2-х факторные)
# Также в отдельную переменную записывается список Stat_BoxList*, содержащий таблицы вида:
#   protein bind expression
# 1     LAM bind  0.8969963
# 2     LAM bind  0.5171169
# 3     LAM bind 11.3465865
# 4     LAM bind  0.0000000
# Это переменная необходима для рассчета корреляции между различными найденными вариантами
#################
Stat_BoxList <- list()
Stat_BoxList_3ds <- list()
for (i in names(GATCvsGenes.List)) {
	for (x in protein_bio_set) {
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
# Вычисление коэффициента достоверности, результат выводится сразу в терминал
# NB! код старый, в файле testing.R  реализовано расчет корреляции и добавление результатов сразу на график
# В переменной Stat_BoxList хранится статистика для двух представлений
# В переменной Stat_BoxList_3ds хранится статистика для трех представлений
for (i in names(Stat_BoxList)) {
	for (y in levels(Stat_BoxList[[i]]$protein)) {
		wilcoxStat <- wilcox.test(Stat_BoxList[[i]][Stat_BoxList[[i]]$protein == y & Stat_BoxList[[i]]$bind == "bind", "expression"], Stat_BoxList[[i]][Stat_BoxList[[i]]$protein == y & Stat_BoxList[[i]]$bind == "not bind", "expression"])
		print(paste(y, " in ", i, ". P-value: ", wilcoxStat$p.value, sep=""))
	}
}

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

# Make Scatter plots comparing gene expression with the protein binding levels at gene promoters and with the averaged protein binding levels at gene bodies
# Вывод статистической информации по сравнению экспрессии генов с уровнем белкового связывания в виде графика распределения (scatter plot's)
for (i in c("BR", "FB", "Kc167")) {
	selectColumns <- grep(paste(i, ".*norm\\.ave", sep="") ,names(GATCvsGenes.List[[i]]), perl=T)
	ifelse(!dir.exists(file.path(prefixDir, outputBio, outputScttr)), dir.create(file.path(prefixDir, outputBio, outputScttr), showWarnings=FALSE), FALSE)
	for (y in selectColumns) {
		scatterName <- sub(paste("(",i, "\\.)(LAM|HP1|PC).*", sep=""), "\\1\\2", names(GATCvsGenes.List[[i]][y]))
		Cor.P <- round(cor(GATCvsGenes.List[[i]][,y], GATCvsGenes.List[[i]][,ncol(GATCvsGenes.List[[i]])], method="pearson", use="pairwise.complete.obs"), digits=2)
		Cor.S <- round(cor(GATCvsGenes.List[[i]][,y], GATCvsGenes.List[[i]][,ncol(GATCvsGenes.List[[i]])], method="spearman", use="pairwise.complete.obs"), digits=2)
		scatter_plot_one <- ggplot(GATCvsGenes.List[[i]], aes(GATCvsGenes.List[[i]][,y], log2(GATCvsGenes.List[[i]][,ncol(GATCvsGenes.List[[i]])]))) + geom_point(alpha=1/10, colour="red", size=3) + geom_text(data = data.frame(), size = 16, hjust=0, aes(mean(GATCvsGenes.List[[i]][,y], na.rm=T), max(log2(GATCvsGenes.List[[i]][,ncol(GATCvsGenes.List[[i]])]), na.rm=T)*0.9, label =c(paste("P.Cor = ", Cor.P, "\n\n", sep=""), paste("S.Cor = ", Cor.S, sep="")))) + theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=35), axis.text.y=element_text(size=35))
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
		scatter_plot_one <- ggplot(GATCvsGenesOrig.List[[i]], aes(GATCvsGenesOrig.List[[i]][,y], log2(GATCvsGenesOrig.List[[i]][,ncol(GATCvsGenesOrig.List[[i]])]))) + geom_point(alpha=1/10, colour="red", size=3) + geom_text(data = data.frame(), size = 16, hjust=0, aes(mean(GATCvsGenesOrig.List[[i]][,y], na.rm=T), max(log2(GATCvsGenesOrig.List[[i]][,ncol(GATCvsGenesOrig.List[[i]])]), na.rm=T)*0.9, label =c(paste("P.Cor = ", Cor.P, "\n\n", sep=""), paste("S.Cor = ", Cor.S, sep="")))) + theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=35), axis.text.y=element_text(size=35))
		bmp(filename=file.path(prefixDir, outputBio, outputScttr, paste("Scatter_plots_", scatterName, "_vs_expression_values. Type_by_full_gene_body.bmp", sep="")), width=800, height=800, units = "px")
		par(mfrow=c(1, 1))
		par(mai=c(1.5, 1.5, 0.7, 0.5))
		par(cex=1.5)
		print(scatter_plot_one)
		rm(Cor.P, Cor.S, scatter_plot_one)
		dev.off()
	}
}

# Find common GATC's by protein
# Пересечение данных DamID разных тканей одного белка между собой. Варианты комбинаций генерируются автоматически при помощи функции combn.
# Но выборка из этих комбинаций осуществляется в ограниченном режиме: подбор комбинаций от 2-х до 5-ти тканей. Результат представлен графически в виде биржевой диаграммы
###############################
options("scipen"=4, "digits"=4)
Stat_ProtList <- list()
for (protein in protein_bio_set){
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
			dfPlus <- dfSel[dfSel[, column_item[1]] == 1 & dfSel[, column_item[2]] == 1 & dfSel[, column_item[3]] == 1 & dfSel[, column_item[4]] == 1 & dfSel[, column_item[5]] == 1, sort(c(column_item, column_item + 1))]
			dfZero <- dfSel[dfSel[, column_item[1]] == 0 & dfSel[, column_item[2]] == 0 & dfSel[, column_item[3]] == 0 & dfSel[, column_item[4]] == 0 & dfSel[, column_item[5]] == 0, sort(c(column_item, column_item + 1))]
		} else {
			print("")
		}
		dfPlus <- na.omit(dfPlus)
		dfZero <- na.omit(dfZero)
		dfCompare <- rbind(dfPlus, dfZero)
		dfCompare <- SplitDfFromColToRow(dfCompare)
		dfCompare$bound <- factor(dfCompare$bound, labels=c("Not bound", "Bound"))
		dfCompare$expression <- log2(dfCompare$expression + 1)
		dfCompare$tissue <- as.factor(dfCompare$tissue)
		ProtName <- paste(levels(dfCompare$tissue), collapse="_vs_")
		Stat_ProtList[[protein]][[ProtName]] <- dfCompare
		boxplot_compare_folder <- "Boxplot_by_common_tissues_domains"
		protein_category_folder <- protein
		ifelse(!dir.exists(file.path(prefixDir, outputBio, outputExpr, boxplot_compare_folder)), dir.create(file.path(prefixDir, outputBio, outputExpr, boxplot_compare_folder), showWarnings=FALSE), FALSE)
		ifelse(!dir.exists(file.path(prefixDir, outputBio, outputExpr, boxplot_compare_folder, protein_category_folder)), dir.create(file.path(prefixDir, outputBio, outputExpr, boxplot_compare_folder, protein_category_folder), showWarnings=FALSE), FALSE)
		# count cases and p-value
		df.B <- cbind(ddply(dfCompare, .(tissue), function(val) nrow(val[val$bound == "Bound",])), y=(ddply(dfCompare, .(tissue), function(val) (fivenum(val[val$bound == "Bound", "expression"]))[4] + max(val$expression) * 0.05))$V1)
		df.UB <- cbind(ddply(dfCompare, .(tissue), function(val) nrow(val[val$bound == "Not bound",])), y=(ddply(dfCompare, .(tissue), function(val) (fivenum(val[val$bound == "Not bound", "expression"]))[4] + max(val$expression) * 0.05))$V1)
		df.pValue <- ddply(dfCompare, .(tissue), function(val) (wilcox.test(val[val$bound == "Bound", "expression"], val[val$bound == "Not bound", "expression"]))$p.value)
		df.pValue$V1 <- paste0("rho==", format(df.pValue$V1, digits=7, width=11), sep="") 
 		pVy.coord <- max(dfCompare$expression) * 0.95
		plot4 <- ggplot(aes(y=expression, x=bound), data=dfCompare) + geom_boxplot(aes(fill=bound)) + labs(title=paste("Gene expression in the ", protein, " common tissues targets ", paste(levels(dfCompare$tissue), collapse=", "), ". P-value count by wilcoxon test.", sep="")) + facet_grid(. ~ tissue) + geom_text(data=df.pValue, aes(x=1.5, y=pVy.coord, label=V1), parse=TRUE) + geom_text(data=df.UB, aes(x=1.25, y, label=V1), parse=TRUE) + geom_text(data=df.B, aes(x=2.25, y, label=V1), parse=TRUE)
		# end of
		pdf(file=file.path(prefixDir, outputBio, outputExpr, boxplot_compare_folder, protein_category_folder, paste(paste("Gene expression in the", protein, "common tissues targets", paste(levels(dfCompare$tissue), collapse=" vs "), sep=" "), "pdf", sep=".")), width=14, height=14)
		print(plot4)
		dev.off()
	}
	rm(dfSel, vecSelnames, columns_list)
}

# Find common GATC's by tissues
# Пересечение данных DamID от разных белков одной ткани между собой.
###############################
Stat_TissList <- list()
for (i in names(GATCvsGenes.List)) {
	
	# Select column with one tissue and all three proteins and then add expression value from this tissue
	for (protein in protein_bio_set){
		vecSel <- grep(paste(i, "\\.", protein, "\\..*target$", sep=""), names(GATCvsGenes.List[[i]]), perl=T, value=T)
		if (exists("vecSelnames") == T) {
			vecSelnames <- c(vecSelnames, vecSel)
		} else {
			vecSelnames <- vecSel
		}
		if (exists("dfSel") == T) {
			dfSel <- cbind(dfSel, data.frame(vecSel=GATCvsGenes.List[[i]][[vecSel]]))
		} else {
			dfSel <- data.frame(vecSel=GATCvsGenes.List[[i]][[vecSel]])
		}
	}
	dfSel <- cbind(dfSel, "expression" = GATCvsGenes.List[[i]][, ncol(GATCvsGenes.List[[i]])], "chr" = GATCvsGenes.List[[i]][, which(names(GATCvsGenes.List$BR) == "chr")])
	colnames(dfSel)[1:(ncol(dfSel)-2)] <- vecSelnames
	
	# Filter bound data by replace "-1" to "0" and remove "NA"
	setOfcolumns <- grep(".*target$", names(dfSel), perl=T)
	for (column in setOfcolumns) {
		if (length(grep("-1", unique(dfSel[, column]), value=T)) == 1) dfSel[dfSel[,column] == -1 & !is.na(dfSel[, column]), column] <- 0 
	}
	
	# Create a variant combinations of proteins between them
	for (y in c(1:3)) {
		if (y == 1) {
			columns_list <- combn(setOfcolumns, y, simplify=F)
		} else {
			columns_list <- append(columns_list, combn(setOfcolumns, y, simplify=F))
		}
	}

	# Comparison of the binding between the various combinations of protein and make report
	for (column_item_list in columns_list) {
		column_item <- unlist(column_item_list)
		if (length(column_item) == 1) {
			dfPlus <- dfSel[dfSel[, column_item[1]] == 1, c(column_item, ncol(dfSel)-1, ncol(dfSel))]
			dfZero <- dfSel[dfSel[, column_item[1]] == 0, c(column_item, ncol(dfSel)-1, ncol(dfSel))]
		} else if (length(column_item) == 2) {
			dfPlus <- dfSel[dfSel[, column_item[1]] == 1 & dfSel[, column_item[2]] == 1, c(column_item, ncol(dfSel)-1, ncol(dfSel))]
			dfZero <- dfSel[dfSel[, column_item[1]] == 0 & dfSel[, column_item[2]] == 0, c(column_item, ncol(dfSel)-1, ncol(dfSel))]
		} else if (length(column_item) == 3) {
			dfPlus <- dfSel[dfSel[, column_item[1]] == 1 & dfSel[, column_item[2]] == 1 & dfSel[, column_item[3]] == 1, c(column_item, ncol(dfSel)-1, ncol(dfSel))]
			dfZero <- dfSel[dfSel[, column_item[1]] == 0 & dfSel[, column_item[2]] == 0 & dfSel[, column_item[3]] == 0, c(column_item, ncol(dfSel)-1, ncol(dfSel))]
		} else {
			print("")
		}
		dfPlus <- na.omit(dfPlus)
		dfZero <- na.omit(dfZero)
		dfCompare <- rbind(dfPlus, dfZero)
		plotTitle <- sub(".*(LAM|HP1|PC).*", "\\1", names(dfCompare)[1:(ncol(dfCompare)-2)], perl=T)
		dfCompare <- dfCompare[, c(1, ncol(dfCompare)-1, ncol(dfCompare))]
		colnames(dfCompare) <- c("bound", "expression", "chr")
		dfCompare$bound <- factor(dfCompare$bound, labels=c("Not bound", "Bound"))
		dfCompare$expression <- log2(dfCompare$expression + 1)
		TissName <- paste(plotTitle, collapse="_vs_")
		Stat_TissList[[i]][[TissName]] <- dfCompare
	}
	boxplot_compare_folder <- "Boxplot_by_common_proteins_domains"
	ifelse(!dir.exists(file.path(prefixDir, outputBio, outputExpr, boxplot_compare_folder)), dir.create(file.path(prefixDir, outputBio, outputExpr, boxplot_compare_folder), showWarnings=FALSE), FALSE)

	# count cases and p-value
	Plot_TissList <- list()
	for (item in names(Stat_TissList[[i]])) {
		df5 <- Stat_TissList[[i]][[item]]
		B <- nrow(df5[df5$bound == "Bound",])
		UB <- nrow(df5[df5$bound == "Not bound",])
		x.coord <- c(1.25, 2.25, 1.5)
		y.coord <- c((fivenum(df5[df5$bound == "Not bound", "expression"]))[4] + max(df5$expression) * 0.05, (fivenum(df5[df5$bound == "Bound", "expression"]))[4] + max(df5$expression) * 0.05, max(df5$expression) * 0.95)
		pValue <- (wilcox.test(df5[df5$bound == "Bound", "expression"], df5[df5$bound == "Not bound", "expression"]))$p.value
		# end of
		
		Plot_TissList[[item]] <- ggplot(aes(y=expression, x=bound, fill=bound), data=df5) + geom_boxplot() + labs(title=item) + annotate("text", x = x.coord, y = y.coord, label = c(UB, B, paste("Wilcoxon test\nP-value: ", format(pValue, digits=7, width=11), sep="")), size=4) + theme(legend.position="none")
	}
	pdf(file=file.path(prefixDir, outputBio, outputExpr, boxplot_compare_folder, paste(paste("Gene expression in the", i, "common protein targets", paste(plotTitle, collapse=" vs "), sep=" "), "pdf", sep=".")), width=14, height=14)
	print(multiplot(plotlist=Plot_TissList, cols=3, layout=matrix(seq(1, 3 * ceiling(7/3)), ncol = 3, nrow = ceiling(7/3), byrow=T)))
	dev.off()
	rm(dfSel, vecSelnames, columns_list)
}
options(scipen = 999)

# Count housekeeping genes located in various combinations of constitutive «anti»domains
# Вычисление количества генов "домашнего хозяйства" в различных комбинациях тканей в пределах одного белка
########################################################################################
if (exists("HKCount") == T | exists("HKGenesFile") == T | exists("OrdCount") == T | exists("OrGenesFile") == T) rm(HKCount, HKGenesFile, OrdCount, OrGenesFile)
HKGenes <- read.csv(HKGenesPath, stringsAsFactors=F)
HKGenes <- GATCvsGenes[GATCvsGenes$Genes.ID %in% HKGenes$ID,]
HKGenes <- cbind(HKGenes[, c("Genes.ID", "chr", "TSS")], HKGenes[, "TSS"])
names(HKGenes) <- c("id", "chr", "start", "end")
idx.genes <- match(HKGenes$id, Genes$ID)
OrdGenes <- Genes[c(1:nrow(Genes))[-idx.genes],c(1:4)]
names(OrdGenes) <- names(HKGenes)
for (comb in names(COMPARE.data.antidomains)) {
	a1 <- COMPARE.data.antidomains[[comb]][, c("seqname", "start", "end")]
	combRanges <- with(a1, IRanges(start, end))

	HKRanges <- with(HKGenes, IRanges(start, end))
	a2 <- as.data.frame(intersect(HKRanges, combRanges))
	HKGenesArea <- as.character(HKGenes[HKGenes$start %in% a2$start, "id"])

	OrdRanges <- with(OrdGenes, IRanges(start, end))
	a3 <- as.data.frame(intersect(OrdRanges, combRanges))
	OrdGenesArea <- as.character(OrdGenes[OrdGenes$start %in% a3$start, "id"])

	stringName <- gsub("(?!M)[a-z\\._0-9MHT]*?\\.anti\\.domains", "", comb, perl=T)

	if (exists("HKCount") == T & exists("HKGenesFile") == T | exists("OrdCount") == T & exists("OrGenesFile") == T) {
		HKCount <-rbind(HKCount, data.frame("stringBetween"=stringName, "count"=length(HKGenesArea)))
		HKGenesFile[[stringName]] <- data.frame(Genes=HKGenesArea, stringsAsFactors=F)

		OrdCount <-rbind(OrdCount, data.frame("stringBetween"=stringName, "count"=length(OrdGenesArea)))
		OrdGenesFile[[stringName]] <- data.frame(Genes=OrdGenesArea, stringsAsFactors=F)

	} else {
		HKCount <- data.frame("stringBetween"=stringName, "count"=length(HKGenesArea))
		HKGenesFile <- setNames(list(x=data.frame(Genes=HKGenesArea, stringsAsFactors=F)), stringName)

		OrdCount <- data.frame("stringBetween"=stringName, "count"=length(OrdGenesArea))
		OrdGenesFile <- setNames(list(x=data.frame(Genes=OrdGenesArea, stringsAsFactors=F)), stringName)

	}
}

# "Интелектуальная" сортировка списка комбинаций тканей. Сначала идут парные комбинации, затем тройные и в конце комбинации из пяти тканей
if (exists("df1") == T | exists("df2") == T | exists("df3") == T | exists("df3") == T) rm(df1, df2, df3, df4)
for (i in c(1:nrow(HKCount))) {
	OccurrencesCount <- CountOccurrencesInCharacter("_vs_", HKCount$stringBetween[i])
	if (OccurrencesCount == 1) {
		if (exists("df1") == T) { df1 <- rbind(df1, HKCount[i, ]) } else { df1 <- HKCount[i, ] }
	} else if (OccurrencesCount == 2) {
		if (exists("df2") == T) { df2 <- rbind(df2, HKCount[i, ]) } else { df2 <- HKCount[i, ] }
	} else if (OccurrencesCount == 3) {
		if (exists("df3") == T) { df3 <- rbind(df3, HKCount[i, ]) } else { df3 <- HKCount[i, ] }
	} else if (OccurrencesCount == 4) {
		if (exists("df4") == T) { df4 <- rbind(df4, HKCount[i, ]) } else { df4 <- HKCount[i, ] }
	} else {
		print("Undefined error!")
	}
}
HKCount <- rbind (df1, df2, df3, df4)

if (exists("df1") == T | exists("df2") == T | exists("df3") == T | exists("df3") == T) rm(df1, df2, df3, df4)
for (i in c(1:nrow(OrdCount))) {
	OccurrencesCount <- CountOccurrencesInCharacter("_vs_", OrdCount$stringBetween[i])
	if (OccurrencesCount == 1) {
		if (exists("df1") == T) { df1 <- rbind(df1, OrdCount[i, ]) } else { df1 <- OrdCount[i, ] }
	} else if (OccurrencesCount == 2) {
		if (exists("df2") == T) { df2 <- rbind(df2, OrdCount[i, ]) } else { df2 <- OrdCount[i, ] }
	} else if (OccurrencesCount == 3) {
		if (exists("df3") == T) { df3 <- rbind(df3, OrdCount[i, ]) } else { df3 <- OrdCount[i, ] }
	} else if (OccurrencesCount == 4) {
		if (exists("df4") == T) { df4 <- rbind(df4, OrdCount[i, ]) } else { df4 <- OrdCount[i, ] }
	} else {
		print("Undefined error!")
	}
}
OrdCount <- rbind (df1, df2, df3, df4)
rm(df1, df2, df3, df4)

# Визуализация полученных данных в виде точечного графика
HKCount$stringBetween <- factor(HKCount$stringBetween, levels=unique(as.character(HKCount$stringBetween)))
write.table(HKCount, file=file.path(prefixDir, outputBio, "Housekeepeng_Genes_counts.csv"), sep=";", row.names=F, col.names=T, quote=F, dec=".", append=F, eol="\r\n")
save(HKGenesFile, file=file.path(prefixDir, outputBio, "Housekeepeng_Genes_list.RData"))
HKplot <- ggplot(HKCount, aes(stringBetween, count)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
png(filename=file.path(prefixDir, outputBio, "Housekeepeng_Genes_plot.png"), width=1600, height=800)
print(HKplot)
dev.off()

OrdCount$stringBetween <- factor(OrdCount$stringBetween, levels=unique(as.character(OrdCount$stringBetween)))
write.table(OrdCount, file=file.path(prefixDir, outputBio, "Ordinary_Genes_counts.csv"), sep=";", row.names=F, col.names=T, quote=F, dec=".", append=F, eol="\r\n")
save(OrdGenesFile, file=file.path(prefixDir, outputBio, "Ordinary_Genes_list.RData"))
Ordplot <- ggplot(OrdCount, aes(stringBetween, count)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
png(filename=file.path(prefixDir, outputBio, "Ordinary_Genes_plot.png"), width=1600, height=800)
print(Ordplot)
dev.off()

# Make Venn Diagram by promotors
################################
idx.gatc.venn <- unlist(sapply(GATCvsGenes$ID, function(x) grep(x, DATA.venn$ID)))

DATA.tss.venn <- DATA.venn[idx.gatc.venn, ]
for (t.s in tissue_bio_set)	MakeVennDiagram(x=DATA.tss.venn, set=paste(t.s,".*bound", sep=""), v1=t.s, v2="LAM, HP1, PC")
for (t.s in tissue_bio_set)	MakeVennDiagram(x=DATA.tss.venn, set=paste(t.s,".*\\.(LAM|HP1).*bound", sep=""), v1=t.s, v2="LAM, HP1")
for (p.s in protein_bio_set) MakeVennDiagram(x=DATA.tss.venn, set=paste(".*\\.", p.s, ".*bound", sep=""), v1=p.s, v2="Kc167, BR, FB, NRN, Glia")
for (p.s in protein_bio_set) MakeVennDiagram(x=DATA.tss.venn, set=paste("(Kc167|BR|FB)\\.", p.s, ".*bound", sep=""), v1=p.s, v2="Kc167, BR, FB")
for (p.s in protein_bio_set) MakeVennDiagram(x=DATA.tss.venn, set=paste("(BR|NRN|Glia)\\.", p.s, ".*bound", sep=""), v1=p.s, v2="BR, NRN, Glia")

# Make scatter plots from expression values
##########################################
for (i in c("BR", "FB")) {
	for (y in c("NRN", "Glia")){
		idx.expr.id <- match(dfExpression.List[[y]]$id, dfExpression.List[[i]]$id)
		e1 <- dfExpression.List[[i]]$expression[idx.expr.id]
		e2 <- dfExpression.List[[y]]$expression
		df.e0 <- data.frame(val1=e1, val2=e2)
		colnames(df.e0) <- c(i, y)
		Cor.P <- round(cor(df.e0[[i]], df.e0[[y]], method="pearson", use="pairwise.complete.obs"), digits=2)
		Cor.S <- round(cor(df.e0[[i]], df.e0[[y]], method="spearman", use="pairwise.complete.obs"), digits=2)
		scatter_plot_three <- ggplot(df.e0, aes(df.e0[[i]], df.e0[[y]]))+geom_point(alpha=1/10, colour="red", size=3) + xlab(paste(i, "expression", sep=" ")) + ylab(paste(y, "expression", sep=" ")) + geom_text(data = data.frame(), size = 4, hjust=0, aes(max(df.e0[[i]], na.rm=T)*0.65, max(df.e0[[y]], na.rm=T)*0.9, label =c(paste("Pearson.Cor = ", Cor.P, "\n\n", sep=""), paste("Spearman.Cor = ", Cor.S, sep="")))) + theme_bw()
		bmp(filename=file.path(prefixDir, outputBio, paste("Scatter_plots_of_expression_in_", i, "_and_", y, ".bmp", sep="")), width=800, height=800, units = "px")
		par(mfrow=c(1, 1))
		par(mai=c(1.5, 1.5, 0.7, 0.5))
		par(cex=1.5)
		print(scatter_plot_three)
		rm(Cor.P, Cor.S, scatter_plot_three)
		dev.off()
	}
}