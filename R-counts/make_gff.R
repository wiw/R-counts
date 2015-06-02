rm(list=ls())

library(Biostrings)
library(BSgenome.Dmelanogaster.UCSC.dm3)

for (i in 1:length(seqnames(Dmelanogaster))){
   GATC.temp <- matchPattern("GATC", DNAString(Dmelanogaster[[i]]), max.mismatch=0, fixed=T)
   GATC <- as.data.frame(matrix(data=NA, nrow=1+length(GATC.temp), ncol=9, byrow=F, dimnames=NULL))
   names(GATC) <- c("chr", "src", "exon", "start", "end", "v1", "v2", "v3", "ID")

   GATC$chr <- sub("chr", "", seqnames(Dmelanogaster)[i])
   GATC$src <- "src"
   GATC$exon <- "exon"

   GATC$start[1] <- 1
   GATC$start[2:(1+length(GATC.temp))] <- GATC.temp@ranges@start

   GATC$end[1:(length(GATC.temp))] <- 3+GATC.temp@ranges@start
   GATC$end[1+length(GATC.temp)] <- length(DNAString(Dmelanogaster[[i]]))
   
   GATC[, c(6:8)] <- "."
	num <- seq(1, nrow(GATC[GATC$chr == sub("chr", "", seqnames(Dmelanogaster)[i]), ]))
   GATC$ID <- paste("ID=gene:DmelGATCr5", seqnames(Dmelanogaster)[i], formatC(num, width=5, flag="0"), sep="")

   rm(GATC.temp)

   if (i == 1) GATCs <- GATC
   if (i > 1)  GATCs <- rbind(GATCs, GATC)
}
options(scipen=10)
write.table(GATCs, "DmelGATCfragments-r5_LP120507.gff", quote=F, sep="\t", row.names=F, col.names=F)
