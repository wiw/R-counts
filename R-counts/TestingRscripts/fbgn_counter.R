#!/usr/bin/R
########################################################################
# FBgn counter
########################################################################
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
save(expressionData, file="/home/anton/data/R-script/R-counts/expressionData.R")