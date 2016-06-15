#!/usr/bin/R
print("Search Biology meaning...")

# Heatmap Correlation
# Считает тепловую карту и корреляции между данными DamID разных наборов. Отличие от предыдущих тепловых карт в том, что добавлены данные по Kc167 клеткам
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
DATA.middle.part <- cbind(DATAs.norm.ave$DATA[,c(1:7)], DATAs.norm.ave$DATA[,grep(paste("(", paste(tissue_bio_set, collapse="|"), ")\\.(", paste(protein_bio_set, collapse="|"), ")\\.(", paste(conditions_bio_set, collapse="|"), ")_(edge|all).*", sep="") ,names(DATAs.norm.ave$DATA), perl=T, value=T)])

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
for (t.s in tissue_bio_set)	MakeVennDiagram(x=DATA.venn, set=paste(t.s,".*bound", sep=""), v1=t.s, v2="LAM, HP1, PC")
for (t.s in tissue_bio_set)	MakeVennDiagram(x=DATA.venn, set=paste(t.s,".*\\.(LAM|HP1).*bound", sep=""), v1=t.s, v2="LAM, HP1")
for (p.s in protein_bio_set) MakeVennDiagram(x=DATA.venn, set=paste(".*\\.", p.s, ".*bound", sep=""), v1=p.s, v2="Kc167, BR, FB, NRN, Glia")
for (p.s in protein_bio_set) MakeVennDiagram(x=DATA.venn, set=paste("(Kc167|BR|FB)\\.", p.s, ".*bound", sep=""), v1=p.s, v2="Kc167, BR, FB")
for (p.s in protein_bio_set) MakeVennDiagram(x=DATA.venn, set=paste("(BR|NRN|Glia)\\.", p.s, ".*bound", sep=""), v1=p.s, v2="BR, NRN, Glia")

# Count area domain from whole genome
# Данные по клеткам Kc уже содержат интерпретацию BioHMM. Поэтому данный код только формирует и форматирует домены и записывает всё в файл
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


# Генерация текстовой статистики по долям доменов от общего размера генома выраженная в п.н.
MakeReportDomainsSize(DOMAIN.data, "constitutive")
MakeReportDomainsSize(DOMAIN.data, "constitutive", includeHet=F)

MakeReportDomainsSize(DOMAIN.data.filt, "filtered")
MakeReportDomainsSize(DOMAIN.data.filt, "filtered", includeHet=F)

MakeReportDomainsSize(DOMAIN.data.anti, "anti")
MakeReportDomainsSize(DOMAIN.data.anti, "anti", includeHet=F)

# Make boxplot & histogram from BioHMM data
# Код генерирует боксплот и гистограммы для каждого сочетания ткань.белок отдельно, которые содержат следующие данные:
# боксплот plot1 - отражает данные DamID по трем состояниям и по хромосомам
# гистограмма plot2 - распределение размеров доменов в т.п.н., ограниченное по оси x до 300
# гистограмма plot3 - распределение размеров доменов в GATC фрагментах, область ограничения по оси x формируется динамически
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

COMPARE.data <- IntersectDomain(inputData=DOMAIN.data, Tset=tissue_bio_set, Pset=protein_bio_set, dscr="constitutive", ScoreValue=1)
COMPARE.data.filt <- IntersectDomain(inputData=DOMAIN.data.filt, Tset=tissue_bio_set, Pset=protein_bio_set, dscr="filt", ScoreValue=1)
COMPARE.data.antidomains <- IntersectDomain(inputData=DOMAIN.data.anti, Tset=tissue_bio_set, Pset=protein_bio_set, dscr="anti", ScoreValue=-1)

# Write statistics from compared data
MakeReportDomainsSize(COMPARE.data, "compare_constitutive")
MakeReportDomainsSize(COMPARE.data, "compare_constitutive", includeHet=F)
MakeReportDomainsSize(COMPARE.data.antidomains, "compare_anti")
MakeReportDomainsSize(COMPARE.data.antidomains, "compare_anti", includeHet=F)


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