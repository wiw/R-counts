edgeDAM1R1 <- read.delim("one_dam1r1_new.sam", header=F, stringsAsFactors=F)
edgeDAM1R1 <- edgeDAM1R1[, c(1,10,11)]
names(edgeDAM1R1) <- c("id", "seq", "qual")
edgeDAM1R1Ori <- read.delim("one_dam1r1_old.sam", header=F, stringsAsFactors=F)
edgeDAM1R1Ori <- edgeDAM1R1Ori[, c(1, 10, 11)]
names(edgeDAM1R1Ori) <- c("id", "seq", "qual")
diffONR1 <- edgeDAM1R1[!(edgeDAM1R1$id %in% edgeDAM1R1Ori$id), 1]
writeLines(diffONR1, con="diff_edge_one_mapped_DAM1R1New_DAM1R1Old.txt", sep="\n")

edgeDAM1R2 <- read.delim("one_dam1r2_new.sam", header=F, stringsAsFactors=F)
edgeDAM1R2 <- edgeDAM1R2[, c(1, 10, 11)]
names(edgeDAM1R2) <- c("id", "seq", "qual")
edgeDAM1R2Ori <- read.delim("one_dam1r2_old.sam", header=F, stringsAsFactors=F)
edgeDAM1R2Ori <- edgeDAM1R2Ori[, c(1, 10, 11)]
names(edgeDAM1R2Ori) <- c("id", "seq", "qual")
diffONR2 <- edgeDAM1R2[!(edgeDAM1R2$id %in% edgeDAM1R2Ori$id), 1]
writeLines(diffONR2, con="diff_edge_one_mapped_DAM1R2New_DAM1R2Old.txt", sep="\n")
