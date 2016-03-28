#!/usr/bin/R
# FBgn counter
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

# Compare DAMId data in BR vs Glia&NRN from LAM&PC&HP1
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

# Make bar by ReportDomainSize Data
###################################
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

fileTwo <- read.csv(file.path(workDir, prefixDir, outputBio, "Report_about_compare_constitutive_without_Heterochromatin_domains_part_from_genome_D.melanogaster.csv"), sep=";")
fileTwo <- fileTwo[order(fileTwo$Item.name),]
Bar_List2 <- list()
for (tset in list(A=tissue_bio_set, B=c("BR", "FB", "Kc167"), C=c("BR", "Glia", "NRN"))) {
	if (exists("r2_df") == T) rm(r2_df)
	for (i in protein_bio_set) {
		sample <- paste(paste(paste(tset, "\\.", i, "\\.[a-zA-Z0-9_]*\\.domains", sep=""), collapse="_vs_"), "$", sep="")
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
print(multiplot(plotlist=Bar_List2, cols=3, layout=matrix(c(1:4), ncol = 2, nrow = 2, byrow=T)))
dev.off()

# Make ScatterPlots between Averaged BR/FB DamID data and BR/FB expression dataframe
####################################################################################
for (i in c("BR", "FB")) {
	s1 <- GATCvsGenes.List[[i]]
	for (y in protein_bio_set) {
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