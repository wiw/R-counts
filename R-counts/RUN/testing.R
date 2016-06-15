#!/usr/bin/R
# Набор функций (4 штук) для собственного расчета данных по экспрессии генов в нейрональных и глиальных тканях

# FBgn counter
# Формирование таблицы expressionData, в которой установлено  соответствие между Flybase ID и Transcript.ID, а также значения экспрессии этих генов в различных тканях
##############
fbgnData <- "/home/anton/Doc/GoogleDrive/Work/SeqCalculate/FlybaseIDs/fbgn_annotation_ID_fb_2014_02.tsv"
expressionData <- "/home/anton/backup/input/ExpressionData/DeSalvo\ 2014\ Front\ Neurosci\ (Suppl)/Table1.csv"

fbgnDF <- read.delim(fbgnData, stringsAsFactors=F)
exprDF <- read.csv2(expressionData, stringsAsFactors=F)

exprDF <- exprDF[grep("^CG[0-9]{2,7}$", exprDF$Transcript.ID, perl=T),]
exprDF$Transcript.ID <- as.factor(exprDF$Transcript.ID)
exprDF <- aggregate(exprDF[,c(5:ncol(exprDF))], by=list(exprDF$Transcript.ID), FUN=mean)
names(exprDF)[1] <- "Transcript.ID"

b1 <- grep("\\\\", fbgnDF$gene_symbol, perl=T, value=T)
fbgnDF <- fbgnDF[!(fbgnDF$gene_symbol %in% b1), ]
fbgnDF <- fbgnDF[grep("^CG[0-9]{2,7}$", fbgnDF$annotation_ID, perl=T),]

idx.fbgn <- unlist(sapply(exprDF$Transcript.ID, function(x) grep(paste("^", x, "$", sep=""), fbgnDF$annotation_ID, perl=T)))
fbgnDF <- fbgnDF[idx.fbgn, ]

exprDF <- exprDF[exprDF$Transcript.ID %in% fbgnDF$annotation_ID, ]
expressionData <- cbind(fbgnDF$primary_Fbgn, exprDF)
expressionData$Brain <- rowMeans(expressionData[, grep("Brain*", names(expressionData))])
expressionData$Elav <- rowMeans(expressionData[, grep("Elav*", names(expressionData))])
expressionData$Repo <- rowMeans(expressionData[, grep("Repo*", names(expressionData))])
expressionData$SG <- rowMeans(expressionData[, grep("SG*", names(expressionData))])
expressionData <- expressionData[, -c(3:22)]
names(expressionData)[1] <- "Flybase.ID"

# Compare DeSalvo Data & Shevelyov data & CustomGenerated Data by DeSalvo Data
##############################################################################
for (i in c("BR", "NRN", "Glia")) {
	if (i == "BR") {
		compareVector <- c("Brain")
	} else if (i == "NRN") {
		compareVector <- c("Brain", "Elav")
	} else {
		compareVector <- c("Brain", "Repo")
	}
	for (y in compareVector) {
		idx.desalvo <- match(expressionData$Flybase.ID, dfExpression.List[[i]]$id)
		x1 <- dfExpression.List[[i]]$expression[idx.desalvo]
		x2 <- expressionData[[y]]
		df.x0 <- data.frame(val1=x1, val2=x2)
		colnames(df.x0) <- c(i, y)
		Cor.P <- round(cor(df.x0[[i]], df.x0[[y]], method="pearson", use="pairwise.complete.obs"), digits=2)
		Cor.S <- round(cor(df.x0[[i]], df.x0[[y]], method="spearman", use="pairwise.complete.obs"), digits=2)
		scatterName <- paste(i, "vs", y, sep="_")
		if (i == "BR") {
			scatterTitle <- ".Seq"
		} else {
			scatterTitle <- ".Shevelyov"
		}
		scatter_plot_four <- ggplot(df.x0, aes(df.x0[[i]], df.x0[[y]])) + geom_point(alpha=1/10, colour="red", size=3) + labs(title=paste("Compare expression data in", paste(i, scatterTitle, sep=""), "versus", paste(y, ".DeSalvo", sep=""), sep=" ")) + xlab(paste(i, "expression", sep=" ")) + ylab(paste(y, "expression", sep=" ")) + geom_text(data = data.frame(), size = 4, hjust=0, aes(max(df.x0[[i]], na.rm=T)*0.65, max(df.x0[[y]], na.rm=T)*0.9, label =c(paste("Pearson.Cor = ", Cor.P, "\n\n", sep=""), paste("Spearman.Cor = ", Cor.S, sep="")))) + theme_bw()
		bmp(filename=file.path(prefixDir, outputBio, paste("Compare_expression_data_in", paste(i, scatterTitle, sep=""), "versus", paste(y, ".DeSalvo.bmp", sep=""), sep="_")), width=800, height=800, units = "px")
		par(mfrow=c(1, 1))
		par(mai=c(1.5, 1.5, 0.7, 0.5))
		par(cex=1.5)
		print(scatter_plot_four)
		dev.off()
	}
}

# Compare DamID data in BR vs Glia&NRN from LAM&PC&HP1
######################################################
for (i in protein_bio_set) {
	for (y in c("NRN", "Glia")) {
		brainName <- grep(paste("BR\\.", i, "*", sep=""), names(DATAs.norm.ave$DATA), perl=T, value=T)
		versusName <- grep(paste(y, "\\.", i, "*", sep=""), names(DATAs.norm.ave$DATA), perl=T, value=T)
		d1 <- DATAs.norm.ave$DATA[[brainName]]
		d2 <- DATAs.norm.ave$DATA[[versusName]]
		df.d0 <- data.frame(val1=d1, val2=d2)
		colnames(df.d0) <- c("BR", y)
		Cor.P <- round(cor(df.d0$BR, df.d0[[y]], method="pearson", use="pairwise.complete.obs"), digits=2)
		Cor.S <- round(cor(df.d0$BR, df.d0[[y]], method="spearman", use="pairwise.complete.obs"), digits=2)
		scatter_plot_five <- ggplot(df.d0, aes(df.d0$BR, df.d0[[y]])) + geom_point(alpha=1/10, colour="red", size=3) + labs(title=paste("DamID values from", i, sep=" ")) + xlab("BR DamID value") + ylab(paste(y, "DamID value", sep=" ")) + geom_text(data = data.frame(), size = 4, hjust=0, aes(min(df.d0$BR, na.rm=T)*1.5, max(df.d0[[y]], na.rm=T) * 0.8, label =c(paste("Pearson.Cor = ", Cor.P, "\n\n", sep=""), paste("Spearman.Cor = ", Cor.S, sep="")))) + theme_bw()
		bmp(filename=file.path(prefixDir, outputBio, paste("DamID_Value_in_BR_vs", y, paste("from_", i, ".bmp", sep=""), sep="_")), width=800, height=800, units = "px")
		par(mfrow=c(1, 1))
		par(mai=c(1.5, 1.5, 0.7, 0.5))
		par(cex=1.5)
		print(scatter_plot_five)
		dev.off()
		rm(df.d0, Cor.P, Cor.S)
	}
}

# Compare Brain.DeSlavo expression Data vs Elav, Repo and BR.seq expression data
################################################################################
OutputFolder <- "CompareExprData_btw_Brain.DeSlavo_vs_Elav&Repo&BR"
ifelse(!dir.exists(file.path(prefixDir, outputBio, outputScttr, OutputFolder)), dir.create(file.path(prefixDir, outputBio, outputScttr, OutputFolder), showWarnings=FALSE), FALSE)
for (i in c("Elav", "Repo", "BR")) {
	if (i == "BR") {
		idx.br <- match(expressionData$Flybase.ID, dfExpression.List[[i]]$id)
		z2 <- dfExpression.List[[i]]$expression[idx.br]
	} else {
		z2 <- expressionData[[i]]
	}
	z1 <- expressionData$Brain
	df.z0 <- data.frame(Val1=z1, Val2=z2)
	colnames(df.z0) <- c("Brain", i)
	Cor.P <- round(cor(df.z0$Brain, df.z0[[i]], method="pearson", use="pairwise.complete.obs"), digits=2)
	Cor.S <- round(cor(df.z0$Brain, df.z0[[i]], method="spearman", use="pairwise.complete.obs"), digits=2)
	scatter_plot_six <- ggplot(df.z0, aes(df.z0$Brain, df.z0[[i]])) + geom_point(alpha=1/10, colour="red", size=3) + labs(title=paste("Compare expression data of Brain.DeSalvo vs", i, sep=" ")) + xlab("Brain DeSalvo expression data") + ylab(paste(i, "expression value", sep=" ")) + geom_text(data = data.frame(), size = 4, hjust=0, aes(min(df.z0$Brain, na.rm = T)*1.5, max(df.z0[[i]], na.rm=T)*0.8, label =c(paste("Pearson.Cor = ", Cor.P, "\n\n", sep=""), paste("Spearman.Cor = ", Cor.S, sep="")))) + theme_bw()
	bmp(filename=file.path(prefixDir, outputBio, outputScttr, OutputFolder, paste("Expression_Value_in_Brain.DeSalvo_vs_", i, ".bmp", sep="")), width=800, height=800, units = "px")
	par(mfrow=c(1, 1))
	par(mai=c(1.5, 1.5, 0.7, 0.5))
	par(cex=1.5)
	print(scatter_plot_six)
	dev.off()
}
################################################################################
################################################################################

# Make bar by ReportDomainSize Data
# Три функции, которые формируют гистограмму о долях доменов (и антидоменах) белков LAM, HP1, PC в единичных тканях, их комбинациях
###################################

# Доли доменов в одиночных тканях
fileOne <- read.csv(file.path(workDir, prefixDir, outputBio, "Report_about_constitutive_without_Heterochromatin_domains_part_from_genome_D.melanogaster.csv"), sep=";")
fileOne <- fileOne[order(fileOne$Item.name),]
outputBar <- "Report domains size"
ifelse(!dir.exists(file.path(prefixDir, outputBio, outputBar)), dir.create(file.path(prefixDir, outputBio, outputBar), showWarnings=FALSE), FALSE)

Bar_List <- list()
for (i in protein_bio_set) {
	r1_names <- grep(paste(".*", i, ".*", sep=""), fileOne$Item.name, perl=T, value=T)
	r1_tnames <- sub("([a-zA-Z0-9]*)\\..*$", "\\1", r1_names, perl=T)
	r1_num <- grep(paste(".*", i, ".*", sep=""), fileOne$Item.name)
	r1_ratio <- round(fileOne[r1_num, "Ratio"], digits=2)
	x.coord <- c(1:5)
	y.coord <- r1_ratio + 2
	Bar_List[[i]] <- ggplot(fileOne[r1_num,], aes(Item.name, Ratio, fill=Item.name)) + geom_bar(stat="identity") + coord_cartesian(ylim = c(0,100)) + labs(title=paste(i, "domains size report", sep=" ")) + xlab("Tissues") + ylab("Domains size ratio") + scale_x_discrete(breaks=r1_names, labels=r1_tnames) + theme_bw() + theme(plot.title=element_text(size=26), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20), axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), legend.position = "none") + annotate("text", x = x.coord, y = y.coord, label = r1_ratio, size=5)
}
pdf(file=file.path(prefixDir, outputBio, outputBar, paste("Bar_report", paste(protein_bio_set, collapse=", "), "constitutive_domains_size_from_different_tissues.pdf", sep="_")), width=14, height=14)
print(multiplot(plotlist=Bar_List, cols=3, layout=matrix(c(1:4), ncol = 2, nrow = 2, byrow=T)))
dev.off()

# Доли доменов в комбинациях тканей
fileTwo <- read.csv(file.path(workDir, prefixDir, outputBio, "Report_about_compare_constitutive_without_Heterochromatin_domains_part_from_genome_D.melanogaster.csv"), sep=";")
fileTwo <- fileTwo[order(fileTwo$Item.name),]
Bar_List2 <- list()
for (tset in list(A=tissue_bio_set, B=c("BR", "FB", "Kc167"), C=c("BR", "Glia", "NRN"), D=c("FB", "Glia", "NRN", "Kc167"), E=c("FB", "Glia", "Kc167"), F=c("BR", "FB", "Glia", "NRN"))) {
	if (exists("r2_df") == T) rm(r2_df)
	for (i in protein_bio_set) {
		sample <- paste("^", paste(paste(tset, "\\.", i, "\\.[a-zA-Z0-9_]*\\.domains", sep=""), collapse="_vs_"), "$", sep="")
		if (exists("r2_df") == T) {
				r2_df <- rbind(r2_df, fileTwo[grep(sample, fileTwo$Item.name, perl=T),])
			} else {
				r2_df <- fileTwo[grep(sample, fileTwo$Item.name, perl=T),]
			}
	}
	r2_df <- r2_df[order(r2_df$Item.name), ]
	r2_names <- r2_df$Item.name
	r2_pnames <- sub(".*\\.(LAM|HP1|PC).*", "\\1", r2_df$Item.name, perl=T)
	r2_ratio <- round(r2_df[, "Ratio"], digits=2)
	x.coord <- c(1:nrow(r2_df))
	y.coord <- r2_ratio + 2
	barlist_item_name <- paste(tset, collapse="_")
	Bar_List2[[barlist_item_name]] <- ggplot(r2_df, aes(Item.name, Ratio, fill=Item.name)) + geom_bar(stat="identity") + coord_cartesian(ylim = c(0,100)) + labs(title=paste("Conservative", paste(tset, collapse=", "), "domains size report", sep=" ")) + xlab(paste(tset, collapse=", ")) + ylab("Conservative domains size ratio") + scale_x_discrete(breaks=r2_names, labels=r2_pnames) + theme_bw() + theme(plot.title=element_text(size=15), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20), axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), legend.position = "none") + annotate("text", x = x.coord, y = y.coord, label = r2_ratio, size=5)
	rm(r2_df)
}
pdf(file=file.path(prefixDir, outputBio, outputBar, paste("Bar_report", paste(protein_bio_set, collapse=", "), "conservative_domains_size_from_", paste(tset, collapse=","), ".pdf", sep="_")), width=14, height=14)
print(multiplot(plotlist=Bar_List2, cols=3, layout=matrix(c(1:6), ncol = 3, nrow = 2, byrow=T)))
dev.off()

# Доли антидоменов в различных комбинациях тканей
fileThree <- read.csv(file.path(workDir, prefixDir, outputBio, "Report_about_compare_anti_without_Heterochromatin_domains_part_from_genome_D.melanogaster.csv"), sep=";")
fileThree <- fileThree[order(fileThree$Item.name),]
Bar_List3 <- list()
for (tset in list(A=tissue_bio_set, B=c("BR", "FB", "Kc167"), C=c("BR", "Glia", "NRN"), D=c("FB", "Glia", "NRN", "Kc167"), E=c("FB", "Glia", "Kc167"), F=c("BR", "FB", "Glia", "NRN"))) {
	if (exists("r3_df") == T) rm(r3_df)
	for (i in protein_bio_set) {
		sample <- paste("^", paste(paste(tset, "\\.", i, "\\.[a-zA-Z0-9_]*\\.anti.domains", sep=""), collapse="_vs_"), "$", sep="")
		if (exists("r3_df") == T) {
				r3_df <- rbind(r3_df, fileThree[grep(sample, fileThree$Item.name, perl=T),])
			} else {
				r3_df <- fileThree[grep(sample, fileThree$Item.name, perl=T),]
			}
	}
	r3_df <- r3_df[order(r3_df$Item.name), ]
	r3_names <- r3_df$Item.name
	r3_pnames <- sub(".*\\.(LAM|HP1|PC).*", "\\1", r3_df$Item.name, perl=T)
	r3_ratio <- round(r3_df[, "Ratio"], digits=2)
	x.coord <- c(1:nrow(r3_df))
	y.coord <- r3_ratio + 2
	barlist_item_name <- paste(tset, collapse="_")
	Bar_List3[[barlist_item_name]] <- ggplot(r3_df, aes(Item.name, Ratio, fill=Item.name)) + geom_bar(stat="identity") + coord_cartesian(ylim = c(0,100)) + labs(title=paste("Conservative", paste(tset, collapse=", "), "anti domains size report", sep=" ")) + xlab(paste(tset, collapse=", ")) + ylab("Conservative antidomains size ratio") + scale_x_discrete(breaks=r3_names, labels=r3_pnames) + theme_bw() + theme(plot.title=element_text(size=15), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20), axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), legend.position = "none") + annotate("text", x = x.coord, y = y.coord, label = r3_ratio, size=5)
	rm(r3_df)
}
pdf(file=file.path(prefixDir, outputBio, outputBar, paste("Bar_report", paste(protein_bio_set, collapse=", "), "conservative_anti_domains_size_from_", paste(tset, collapse=","), ".pdf", sep="_")), width=14, height=14)
print(multiplot(plotlist=Bar_List3, cols=3, layout=matrix(c(1:6), ncol = 3, nrow = 2, byrow=T)))
dev.off()

# Make ScatterPlots between Averaged BR/FB DamID data and BR/FB expression dataframe
####################################################################################
for (i in c("BR", "FB")) {
	s1 <- GATCvsGenes.List[[i]]
	protein_ext_set <- c(protein_bio_set, "DAM")
	for (y in protein_ext_set) {
		ave.s <- s1[, grep(paste(i, "\\.", y, "\\.(", paste(conditions_bio_set, collapse="|"), ")_(all|edge)\\.ave$", sep=""), names(s1), perl=T)]
		expr.s <- s1[, ncol(s1)]
		ave.s <- log2(ave.s + 1)
		expr.s <- log2(expr.s + 1)
		df.s0 <- data.frame(val1=ave.s, val2=expr.s)
		colnames(df.s0) <- c("ave", "expr")
		Cor.P <- round(cor(df.s0$ave, df.s0$expr, method="pearson", use="pairwise.complete.obs"), digits=2)
		Cor.S <- round(cor(df.s0$ave, df.s0$expr, method="spearman", use="pairwise.complete.obs"), digits=2)
		scatter_plot_seven <- ggplot(df.s0, aes(ave, expr))+geom_point(alpha=1/10, colour="red", size=3, na.rm=T) + xlab(paste(i, y, "DamID_averaged_data: log2(value +1)", sep="_")) + ylab(paste(i, y, "expression_data: log2(value +1)", sep="_")) + geom_text(data = data.frame(), size = 10, hjust=0, aes(max(df.s0$ave, na.rm=T)*0.5, max(df.s0$expr, na.rm=T)*0.9, label =c(paste("Pearson.Cor = ", Cor.P, "\n\n", sep=""), paste("Spearman.Cor = ", Cor.S, sep="")))) + theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=22), axis.text.y=element_text(size=22), legend.position = "none")
		bmp(filename=file.path(prefixDir, outputBio, paste("Scatter_plots_DamID_averaged_data_vs_expression_values_in_", i, ".", y, ".bmp", sep="")), width=800, height=800, units = "px")
		par(mfrow=c(1, 1))
		par(mai=c(1.5, 1.5, 0.7, 0.5))
		par(cex=1.5)
		print(scatter_plot_seven)
		dev.off()
	}
}

# tid:1 Gene expression from protein binding by tissues, proteins and chromosomes
local({
taskID <- "tid_1"
ifelse(!dir.exists(file.path(prefixDir, outputBio, outputExpr, taskID)), dir.create(file.path(prefixDir, outputBio, outputExpr, taskID), showWarnings=FALSE), FALSE)
for (i in c("BR", "FB")) {
	tid1_list <-list()
	for (y in names(Stat_TissList[[i]])) {
		df6 <- Stat_TissList[[i]][[y]]
		df6 <- df6[df6$chr %in% c("X", "2L", "2R", "3L", "3R", "4"), ]
		df6 <- df6[order(df6$chr), ]
		df6$chr <- factor(df6$chr, levels=unique(df6$chr))
		tid1_item_name <- paste(i, y, sep="_")
		z1 <- na.replace(as.vector(by(df6, df6[, c(1, 3)], nrow)), 0) # Count bound/unbound event for every chromosomes
		z2 <- paste("pV:", format(as.vector(by(df6, df6$chr, function(x) if (length(unique(x[, "bound"])) < 2) return("Not enough data") else (wilcox.test(expression ~ bound, x))$p.value)), digits = 3, width = 8), sep=" ") # Count p.value between bound/unbound expr.Values by chromosomes
		chr_stat <- as.character(rbind(matrix(z1, nrow=2, byrow=F), matrix(z2, nrow=1, byrow=F))) # Combine z1 an z2 like as "z1.1 z1.2 z2.1 z1.3 z1.4 z2.2" & etc...
		x.coord <- as.vector(matrix(c(seq(0.70, by=1, length.out=6), seq(1.10, by=1, length.out=6), seq(1.05, by=1, length.out=6)), nrow=3, byrow=T))
		y.coord <- na.replace(as.vector(rbind(matrix(as.vector(by(df6$expression, df6[, c(1,3)], function(x) (fivenum(x))[4])) + max(df6$expression) * 0.05, nrow=2, byrow=F), matrix(rep(max(df6$expression) * 0.95, 6), nrow=1, byrow=F))), 0)
		tid1_list[[tid1_item_name]] <- ggplot(data=df6, aes(x=chr, y=expression, fill=bound)) + geom_boxplot() + labs(title=tid1_item_name) + annotate("text", x = x.coord, y = y.coord, label = chr_stat, size = 3, na.rm=T) + theme_bw() + theme(axis.title.x=element_blank())
	}
	pdf(file=file.path(prefixDir, outputBio, outputExpr, taskID, paste(taskID, " Gene expression vs binding by chromosomes in ", i, ".pdf", sep="")), width=14, height=14)
	print(multiplot(plotlist=tid1_list, cols=1, layout=matrix(c(1:7), ncol = 1, nrow = 7, byrow=T)))
	dev.off()
}
})

# tid:2 Compare underreplication data of larva fat bodies with our DamID data
options(scipen = 999)
taskID <- "tid_2"
ifelse(!dir.exists(file.path(prefixDir, outputBio, outputExpr, taskID)), dir.create(file.path(prefixDir, outputBio, outputExpr, taskID), showWarnings=FALSE), FALSE)

gatcUnd <- cbind(gatcG[with(gatcG, !(chr %in% "U")), c(1:4)], DATAs.norm.ave$DATA[, grep("FB.*", names(DATAs.norm.ave$DATA))])
gatcUnd <- gatcUnd[with(gatcUnd, order(chr)), ]
gatcUnd$chr <- factor(gatcUnd$chr, levels = unique(gatcUnd$chr))
tid2_fx0 <- function(GPL, NDF, NORM) { 
# GPL - data with chromosome coordinates and DMEL ID's
# NDF - DMEL ID's and microarray x.y position
# NORM - normalized undereplication data with value and microarray x.y position
	gpl <- read.delim(gzfile(GPL), skip=9, stringsAsFactors=F)
	ndf <- read.delim(gzfile(NDF), stringsAsFactors=F)
	norm <- read.delim(gzfile(NORM), skip=1, stringsAsFactors=F)
	norm <- norm[with(norm, order(Y, X)), ]
	gpl <- gpl[with(gpl, which(FEATURE %in% ndf$PROBE_ID)), ]
	out <- data.frame("ID" = ndf$PROBE_ID, "chr" = gpl$CHROMOSOME, "start" = (gpl$RANGE_START + gpl$RANGE_END) / 2, "end" = (gpl$RANGE_START + gpl$RANGE_END) / 2, "ratio" = norm$NormalizedLog2Ratio, stringsAsFactors = F)
	out$chr <- sub("chr(.*)", "\\1", out$chr, perl=T)
	out <- out[with(out, !(chr %in% "U")), ]
	out <- out[with(out, order(chr)), ]
	out$chr <- factor(out$chr, levels=unique(out$chr))
	return(out)
}
UndData.list <- list("data1M"=tid2_fx0("/home/anton/backup/input/URD/GPL11023-20784.txt.gz", "/home/anton/backup/input/URD/GPL11023_252404810019_pseudo_ndf.txt.gz", "/home/anton/backup/input/URD/GSM697657_MA2C_Orr-Weaver_Fatbody_1M_normalized.txt.gz"), "data400K"=tid2_fx0("/home/anton/backup/input/URD/GPL11290-35146.txt.gz", "/home/anton/backup/input/URD/GPL11290_252532110008_pseudo_ndf.txt.gz", "/home/anton/backup/input/URD/GSM697658_MA2C_Orr-Weaver_Fatbody_400K_normalized.txt.gz"))
UndMeanData.list <- list()
for (i in names(UndData.list)) {
	und_vs_gatc <- AssignGATCsToGenes(UndData.list[[i]], gatcUnd, und=T, tissue.damid="FB")
	UndMeanData.list[[i]] <- ddply(und_vs_gatc, "ID", numcolwise(mean, na.rm = TRUE))
}
Und_vs_GatcDF <- with(UndMeanData.list, rbind(data1M, data400K))
Und_vs_GatcDF <- Und_vs_GatcDF[with(Und_vs_GatcDF, order(ID)), ]
Und_vs_GatcDF <- Und_vs_GatcDF[(duplicated(Und_vs_GatcDF$ID) | duplicated(Und_vs_GatcDF$ID, fromLast=T)) == T, ]
Und_vs_GatcDF <- ddply(Und_vs_GatcDF, "ID", numcolwise(mean, na.rm = TRUE))
Und_vs_GatcDF <- Und_vs_GatcDF[rowSums(is.na(Und_vs_GatcDF[, c(4:6)]))!=3, ]
Und_vs_GatcDF[is.na(Und_vs_GatcDF)] <- NA
tid2_list0 <- list()
tid2_fx1 <- function(data) {
	tempList <- list()
	for (i in protein_bio_set) {
		set <- grep(paste(i, "|ratio", sep=""), names(data))
		Cor.P <- paste("Pearson.Cor = ", round(cor(data[, set[1]], data[, set[2]], method="pearson", use="pairwise.complete.obs"), digits=2), sep = "")
		Cor.S <- paste("Spearman.Cor = ", round(cor(data[, set[1]], data[, set[2]], method="spearman", use="pairwise.complete.obs"), digits=2), sep = "")
		allValues <- paste("Values = ", nrow(data[!is.na(data[, set[1]]),]), sep = "")
		x.coord <- rep(min(data[, set[1]], na.rm = T) * 0.5, 3)
		y.coord <- seq(min(data[, set[2]], na.rm = T) * 0.7, by=-0.7, length.out=3)
		tempList[[i]] <- ggplot(data, aes(data[, set[1]], data[, set[2]])) + geom_point(alpha = 1/10, colour = "#65A6D1", size = 3, na.rm = T) + xlab(paste(i, "Fat bodies DamID data", sep=" ")) + ylab("Undereplication data in fat bodies") + annotate("text", x = x.coord, y = y.coord, label = c(Cor.P, Cor.S, allValues), size = 6, na.rm=T) + theme_bw() + theme(axis.title.x = element_text(size = 26), axis.title.y = element_text(size = 26), axis.text.x = element_text(size = 22), axis.text.y = element_text(size = 22))
	}
	return(tempList)
}
tid2_list0 <- tid2_fx1(Und_vs_GatcDF)
bmp(filename=file.path(prefixDir, outputBio, outputExpr, taskID, "Scatter_plots_DamID_data_vs_undereplication_values_in_Fat_bodies.bmp"), width=800*3, height=800, units = "px")
print(multiplot(plotlist=tid2_list0, cols=3, layout=matrix(c(1:3), ncol = 3, nrow = 1, byrow=T)))
dev.off()

fb_und_domain <- read.delim("/home/anton/backup/input/URD/fb_und_domain.txt")
gatcUnd_mod <- gatcUnd
gatcUnd_mod[, c(3, 4)] <- (gatcUnd_mod$start + gatcUnd_mod$end) / 2
UndDomain_vs_GatcDF <- AssignGATCsToGenes(fb_und_domain, gatcUnd_mod, und=T, use.und.domain=T)
UndDomain_vs_GatcDF$fb_und_domain <- factor(UndDomain_vs_GatcDF$fb_und_domain)
tid2_list1 <- list()
tid2_fx2 <- function(data) {
	tempList <- list()
	for (i in protein_bio_set) {
		set <- grep(paste(i, "|fb_und_domain", sep=""), names(data))
		Rep <- data[!is.na(data[, set[1]]) & data[,set[2]] == 1, ]
		UndRep <- data[!is.na(data[, set[1]]) & data[,set[2]] == 0, ]
		und_stat <- c(nrow(Rep), nrow(UndRep), paste("P.value: ", (wilcox.test(Rep[, set[1]], UndRep[, set[1]]))$p.value, sep=""))
		x.coord <- c(2.20, 1.20, 1.5)
		y.coord <- c(fivenum(Rep[, set[1]])[4] + max(data[, set[1]], na.rm=T) * 0.05, fivenum(UndRep[, set[1]])[4] + max(data[, set[1]], na.rm=T) * 0.05, max(data[, set[1]], na.rm=T) * 0.95)
		tempList[[i]] <- ggplot(data=data, aes(x=data[, set[2]], y=data[, set[1]], fill=data[, set[2]])) + geom_boxplot(na.rm=T) + xlab("Undereplication domains in fat bodies") + ylab(paste(i, "Fat bodies DamID data", sep=" ")) + annotate("text", x = x.coord, y = y.coord, label = und_stat, size = 5, na.rm=T) + theme_bw() + theme(axis.title.x=element_blank()) + theme(axis.title.x = element_text(size = 26), axis.title.y = element_text(size = 26), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), legend.position = "none") + scale_x_discrete(breaks=c("0", "1"), labels=c("UndRep", "Rep"))
	}
	return(tempList)
}
options("scipen"=4, "digits"=4)
tid2_list1 <- tid2_fx2(UndDomain_vs_GatcDF)
bmp(filename=file.path(prefixDir, outputBio, outputExpr, taskID, "Boxplot_DamID_data_vs_undereplication_domains_in_Fat_bodies.bmp"), width=800*3, height=800, units = "px")
print(multiplot(plotlist=tid2_list1, cols=3, layout=matrix(c(1:3), ncol = 3, nrow = 1, byrow=T)))
dev.off()	
options(scipen = 999)

# # tid:3 Determind pericentromeric region and make DamID targets vs gene expression Boxplot_DamID_data_vs_undereplication_domains_in_Fat_bodies
# taskID <- "tid_3"
# ifelse(!dir.exists(file.path(prefixDir, outputBio, outputExpr, taskID)), dir.create(file.path(prefixDir, outputBio, outputExpr, taskID), showWarnings=FALSE), FALSE)

# names(GATCvsGenes.List)
# # [1] "BR"    "FB"    "Glia"  "NRN"   "Kc167"

# str(GATCvsGenes.List$BR)
# # 'data.frame':   14832 obs. of  20 variables:
# #  $ ID                   : chr  "r5GATC2L00040" "r5GATC2L00081" "r5GATC2L00091" "r5GATC2L00173" ...
# #  $ chr                  : chr  "2L" "2L" "2L" "2L" ...
# #  $ start                : int  6923 21359 25086 58765 66423 71991 75585 86317 94101 102023 ...
# #  $ end                  : int  7693 21643 26562 59732 67274 72652 76557 87735 95132 102687 ...
# #  $ Genes.ID             : Factor w/ 14832 levels "FBgn0000003",..: 4322 380 4323 11924 13261 4324 4325 486 4326 4327 ...
# #  $ TSS                  : int  7529 21372 25151 59242 67044 72388 76446 87382 94752 102382 ...
# #  $ strand               : Factor w/ 2 levels "-","+": 2 1 1 1 2 2 2 1 2 2 ...
# #  $ BR.LAM.m_all.norm.ave: num  -3.162 NA -1.456 0.397 -3.651 ...
# #  $ BR.LAM.m_all.ave     : num  9.5 0 6.5 33.5 19 10.5 33 147 5.5 2.5 ...
# #  $ BR.LAM.m_target      : num  -1 NA -1 0 -1 -1 -1 -1 -1 NA ...
# #  $ BR.LAM.m_target.filt : num  -1 NA -1 0 -1 -1 -1 -1 -1 NA ...
# #  $ BR.HP1.m_all.norm.ave: num  1.2302 NA 1.8135 1.8952 -0.0889 ...
# #  $ BR.HP1.m_all.ave     : num  196 0 59 91 216 ...
# #  $ BR.HP1.m_target      : num  0 NA 1 0 0 0 0 0 -1 -1 ...
# #  $ BR.HP1.m_target.filt : num  0 NA 1 0 0 0 0 0 -1 -1 ...
# #  $ BR.PC.m_all.norm.ave : num  NA NA NA -1.16 -4.36 ...
# #  $ BR.PC.m_all.ave      : num  0 0 0 23 54.5 24 74 112 7.5 0 ...
# #  $ BR.PC.m_target       : num  NA NA NA 0 -1 -1 -1 0 -1 NA ...
# #  $ BR.PC.m_target.filt  : num  NA NA NA 0 -1 -1 -1 0 -1 NA ...
# #  $ expression           : num  23.9 5471.4 25.5 697.1 1910.3 ...

# Статистика распределения размера доменов белков в одиночных тканях...
for (i in protein_bio_set) {
	stat_list <- list()
	if (exists("stat_df") == T) rm(stat_df)
	for (y in tissue_bio_set) {
		names_select <- grep(paste0(y, "\\.", i, "\\..*domains"), names(DOMAIN.data), perl=T, value=T)
		temp_df <- DOMAIN.data[[names_select]][DOMAIN.data[[names_select]]$seqname %in% c("X", "2L", "2R", "3L", "3R", "4"), ]
		temp_df <- temp_df$end-temp_df$start
		stat_list[[names_select]] <- log10(as.data.frame(temp_df))
		if (exists("stat_df") == T) {
			stat_df <- rbind(stat_df, data.frame("name" = sub("([a-zA-Z0-9]*)(\\..*)", "\\1", names_select, perl=T), "count" = length(temp_df), "mean" = mean(temp_df), "median" = median(temp_df), stringsAsFactors=F))
		} else {
			stat_df <- data.frame("name" = sub("([a-zA-Z0-9]*)(\\..*)", "\\1", names_select, perl=T), "count" = length(temp_df), "mean" = mean(temp_df), "median" = median(temp_df), stringsAsFactors=F)
		}
	}
	# x.coord <- c(rep(40000, 5), rep(60000, 5))
	# y.coord <- rep(seq(from=40000, by=-2000, length.out=5), 2)
	# legend.text <- c(stat_df$name, stat_df$count)
	stat.color <- c("black", "orange", "skyblue" , "lightgreen", "indianred")
	tempdd <- ggplot(NULL, aes(x=temp_df)) + geom_density(colour="black", data=stat_list[[1]]) + geom_density(colour="orange", data=stat_list[[2]]) + geom_density(colour="skyblue", data=stat_list[[3]]) + geom_density(colour="lightgreen", data=stat_list[[4]]) + geom_density(colour="indianred", data=stat_list[[5]]) + geom_vline(aes(xintercept=log10(stat_df$mean)), colour=stat.color) + geom_vline(aes(xintercept=log10(stat_df$median)), colour=stat.color, linetype="longdash") + theme_bw() + theme(axis.title.x=element_blank(), axis.title.y = element_text(size = 26), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), legend.position = "none")
	pdf(file=file.path(prefixDir, outputBio, outputBar, paste0(i, "_statistics in_", paste(tissue_bio_set, collapse="_"), ".pdf")), width=14, height=7)
	print(tempdd)
	dev.off()
	write.csv(stat_df, file=file.path(prefixDir, outputBio, outputBar, paste0(i, "_table_stat_in_", paste(tissue_bio_set, collapse="_"), ".csv")), row.names=F)
}


# ... и в комбинациях тканей
# отдельно формируется текстовый файл содержащий информацию о количестве найденных доменов, среднем и медианном размере этих доменов
for (i in list(A=tissue_bio_set, B=c("BR", "FB", "Kc167"), C=c("BR", "Glia", "NRN"), D=c("FB", "Glia", "NRN", "Kc167"), E=c("FB", "Glia", "Kc167"), F=c("BR", "FB", "Glia", "NRN"), G="BR", H="FB", I="Glia", J="NRN", K="Kc167")) {
	stat_list2 <- list()
	if (exists("stat_df") == T) rm(stat_df)
	for (y in protein_bio_set) {
		if (length(i) > 1) {
			names_select <- grep(paste0("^", paste(paste0(i, "\\.", y, "\\.[a-zA-Z0-9_]*\\.domains"), collapse="_vs_"), "$"), names(COMPARE.data), value=T, perl=T) 
			temp_df <- COMPARE.data[[names_select]][COMPARE.data[[names_select]]$seqname %in% c("X", "2L", "2R", "3L", "3R", "4"), ]
		} else {
			names_select <- grep(paste0(i, "\\.", y, "\\..*domains"), names(DOMAIN.data), perl=T, value=T)
			temp_df <- DOMAIN.data[[names_select]][DOMAIN.data[[names_select]]$seqname %in% c("X", "2L", "2R", "3L", "3R", "4"), ]
		}
		temp_df <- temp_df$end-temp_df$start
		stat_list2[[names_select]] <- log10(as.data.frame(temp_df))
		if (exists("stat_df") == T) {
			stat_df <- rbind(stat_df, data.frame("name" = y, "count" = length(temp_df), "mean" = mean(temp_df), "median" = median(temp_df), stringsAsFactors=F))
		} else {
			stat_df <- data.frame("name" = y, "count" = length(temp_df), "mean" = mean(temp_df), "median" = median(temp_df), stringsAsFactors=F)
		}
	}
	stat.color <- c("black", "orange", "skyblue")
	tempdd <- ggplot(NULL, aes(x=temp_df)) + geom_density(colour="black", data=stat_list2[[1]]) + geom_density(colour="orange", data=stat_list2[[2]]) + geom_density(colour="skyblue", data=stat_list2[[3]]) + geom_vline(aes(xintercept=log10(stat_df$mean)), colour=stat.color) + geom_vline(aes(xintercept=log10(stat_df$median)), colour=stat.color, linetype="longdash") + theme_bw() + theme(axis.title.x=element_blank(), axis.title.y = element_text(size = 26), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), legend.position = "none")
	pdf(file=file.path(prefixDir, outputBio, outputBar, paste0(paste(protein_bio_set, collapse="_"), "_statistics in_", paste(i, collapse="_"), ".pdf")), width=14, height=7)
	print(tempdd)
	dev.off()
	names(stat_df)[1] <- paste(i, collapse="_")
	write.csv(stat_df, file=file.path(prefixDir, outputBio, outputBar, paste0(paste(protein_bio_set, collapse="_"), "_table_stat_in_", paste(i, collapse="_"), ".csv")), row.names=F)
}

# Вычисление количества промоторов разного типа (связаны, несвязаны, NA, сомнительны) в данных по каждому белку каждой ткани
for (i in names(GATCvsGenes.List)) {
	for (prot in protein_bio_set) {
		if (exists("x0") == T) rm(x0)
		for (y in c(-1, 0, 1)) {
			prot_name <- grep(paste0("^", i, "\\.", prot, "\\..*target$"), names(GATCvsGenes.List[[i]]), value=T)
			if (exists("x0") == T) {
				x0 <- c(x0, length(GATCvsGenes.List[[i]][[prot_name]][GATCvsGenes.List[[i]][[prot_name]] == y & !is.na(GATCvsGenes.List[[i]][[prot_name]])]))
			} else {
				x0 <- length(GATCvsGenes.List[[i]][[prot_name]][GATCvsGenes.List[[i]][[prot_name]] == y & !is.na(GATCvsGenes.List[[i]][[prot_name]])])
			}
		}
		x0 <- c(x0, nrow(GATCvsGenes.List[[i]])-sum(x0), nrow(GATCvsGenes.List[[i]]))
		if (exists("df0") == T) {
			df0 <- cbind(df0, data.frame("v1" = x0))
		} else {
			df0 <- data.frame("v1" = x0)
		}
		colnames(df0)[ncol(df0)] <- paste0(i, ".", prot)
	}
}
rownames(df0) <- c("TSS unbound_-1", "TSS amb._0", "TSS bound_1", "TSS NA", "All TSS's")
write.csv(df0, file=file.path(prefixDir, outputBio, "TSS_promotors_stat.csv"), row.names=T)
