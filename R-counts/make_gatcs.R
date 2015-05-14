rm(list=ls())

library(Biostrings)
library(BSgenome.Dmelanogaster.UCSC.dm3)

for (i in 1:length(seqnames(Dmelanogaster))){
   GATC.temp <- matchPattern("GATC", DNAString(Dmelanogaster[[i]]), max.mismatch=0, fixed=T)
   GATC <- as.data.frame(matrix(data=NA, nrow=1+length(GATC.temp), ncol=7, byrow=F, dimnames=NULL))
   names(GATC) <- c("ID", "chr", "start", "end", "width", "presence.ma", "ID.il")

   GATC$chr <- sub("chr", "", seqnames(Dmelanogaster)[i])

   GATC$start[1] <- 1
   GATC$start[2:(1+length(GATC.temp))] <- GATC.temp@ranges@start

   GATC$end[1:(length(GATC.temp))] <- 3+GATC.temp@ranges@start
   GATC$end[1+length(GATC.temp)] <- length(DNAString(Dmelanogaster[[i]]))
   
   num <- seq(1, nrow(GATC[GATC$chr == sub("chr", "", seqnames(Dmelanogaster)[i]), ]))

   GATC$ID <- paste("r5GATC", seqnames(Dmelanogaster)[i], formatC(num, width=5, flag="0"), sep="")
   rm(GATC.temp)

   if (i == 1) GATCs <- GATC
   if (i > 1)  GATCs <- rbind(GATCs, GATC)
}
GATCs$width <- GATCs$end - GATCs$start + 1

GATCs.temp <- GATCs[GATCs$width != 8,]
pr.ma <- GATCs[GATCs$width == 8, ]

gatc <- read.delim("GATCs.txt", header=T, stringsAsFactors=F, sep="\t")

GATCs.temp$presence.ma <- gatc$presence.ma
GATCs.temp$ID.il <- gatc$ID.il 

GATCs.temp <- rbind(GATCs.temp, pr.ma)

GATCs.temp[is.na(GATCs.temp$presence.ma) == T, 6] <- 0
GATCs <- GATCs.temp[order(GATCs.temp$ID), ]
write.table(GATCs, "GATCs_mod.txt", quote=F, sep="\t", row.names=F, col.names=T)
