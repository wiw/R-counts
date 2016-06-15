#!/usr/bin/R
###############################
# Start calculate BioHMM data #
###############################

# By sample use only DATA, without pseudo counts
# Интерпретация нормированных усредненных данных DamID в области связанные и несвязанные с исследуемым белком.
# Сырые данные записываются в список HMM.data, интерпретированные домены записываются в список DOMAIN.data, 
# фильтрованные домены размером более 3 GATC фрагментов записываются в список DOMAIN.data.filt, 
# антидомены записываются в список DOMAIN.data.anti

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

		# Интерпретация нормированных усредненных данных DamID
		print(paste("I'm run BioHMM data from chromosome", chromosomesVector[chr], sep=" "))
		HMM.data[[DATA.name]][[dataframe.name]]$BioHMM.output <- runBioHMM(mval=HMM.data[[DATA.name]][[dataframe.name]]$DamID.value, datainfo=HMM.data[[DATA.name]][[dataframe.name]], useCloneDists=T)
		
		# Классификатор - определяет число возможных вариантов интерпретации (сейчас три: связано, сомнительно, несвязано) и среднее значение DamID для каждого варианта.
		# Соответственно, минимальное значения классификатора определяется как -1, максимальное - как 1, среднее значение - 0
		classifier <- aggregate(x=HMM.data[[DATA.name]][[dataframe.name]]$DamID.value, by=list(HMM.data[[DATA.name]][[dataframe.name]]$BioHMM.output), FUN=mean)
 		
 		# Поиск, какое значение HMM соответствует минимальному или максимальному значению DamID из классификатора
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
			# Фильтруем данные по размеру домена, не менее 3-х GATC фрагментов
			tmp_rle <- rle(HMM.data[[DATA.name]][[dataframe.name]]$target)
			tmp_rle$values[tmp_rle$lengths <= 2 & tmp_rle$values == 1] <- 0
			HMM.data[[DATA.name]][[dataframe.name]]$target.filt <- inverse.rle(tmp_rle)

			# Формирую домены, форматирую и записываю всё в GFF файл
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