#!/usr/bin/R

# Load libraries and set options
library(gplots)
library(ggplot2)
library(snapCGH)
library(cluster)
library(limma)
library(Vennerable)
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(DESeq)
library(Ringo)
library(tools)
library(plyr)
library(gridExtra)
options(scipen = 999)

# Declare variables
###################
prefixDir <- paste("RUN", format(Sys.time(), "%d-%m-%Y"), "ludo", sep="_") # output directory into working directory
onlyEdge <- F # use only edge reads to counts or not
workDir <- "/home/anton/data/R-script/R-counts"	# working directory (WD)
sourceDir <- "/home/anton/backup/ludo-output" # location your RData files. You can specify the highest folder as it is possible. Searching runs recursively.
damIdLocation <- "/home/anton/data/DAM/RUN/damid_description.csv" # location your DamID-Description file
genesFilePath <- "/home/anton/backup/2013-07-09_Assigning_GATCs_to_Genes_(FBgn_IDs)/Drosophila_melanogaster.BDGP5.25.64_Genes_AP130708.txt"
HKGenesPath <- "/home/anton/backup/2013-07-09_Assigning_GATCs_to_Genes_(FBgn_IDs)/houskeeping_genes.csv"
exprDataDir <- "/home/anton/backup/input/ExpressionData"
outputGff <- "gff"	# output folder for gff in WD
outputWig <- "wig"	# output folder for wig in WD
outputScttr <- "scatter_plots"	# output folder for scatter plots in WD
outputDomain <- "domains"
outputCleanStat <- "clean_stat"
outputBio <- "Bio"
outputHeatmap <- "Heatmap"
outputBioBoxplot <- "Boxplot"
outputExpr <- "Expression"
startCol <- 7	# the number of last column in GATCs file, default "7"
gatcFile <- paste(workDir, "GATCs_mod.txt", sep="/")	# location you GATCs file
gatcFile4Genes <- paste(workDir, "GATCs.txt", sep="/")
needCombine <- F	# are you need to combine some columns into one (T or F)? we recommend set the "F"; if you select "T" - edit the "combine" vector on string #85 into this code
usePseudoCounts <- F	# are you need to add pseudo counts into source data (T or F)? Default "T"
pseudoCounts <- c(0.01)		# the vector of pseudo counts, default "c(1)"
corrMethod <- c("pearson", "spearman")	# the vector of type correlation method - please don't change this setting
heatmapColors <- greenred(200)	# check you color into heatmap
labelHeatmap <- c("B", "C")	# the vector of heatmap label - please don't change this setting
labelAcf <- c("A_ALL", "B_MA")	# the vector of acf label - please don't change this setting
writeTemp <- T		# use this option if you need to get the intermediate files, default "T"
needSomeFiles <- F	# if need calculate only some files from source list, default "F"
useSomeFiles <- c(9:16)	# the region of files which need calculate
tissue_bio_set <- c("BR", "FB", "Glia", "NRN", "Kc167") # Part of tissies from all
protein_bio_set <- c("LAM", "HP1", "PC") # Part of protein from all
conditions_bio_set <- c("m", "m_25mkM4HT", "mf_min", "std")

# Create folders
################
dir.create(file.path(workDir, prefixDir), showWarnings = FALSE)
dir.create(file.path(workDir, prefixDir, outputGff), showWarnings = FALSE)
dir.create(file.path(workDir, prefixDir, outputWig), showWarnings = FALSE)
dir.create(file.path(workDir, prefixDir, outputScttr), showWarnings = FALSE)
dir.create(file.path(workDir, prefixDir, outputGff, outputDomain), showWarnings = FALSE)
dir.create(file.path(workDir, prefixDir, outputCleanStat), showWarnings = FALSE)
dir.create(file.path(workDir, prefixDir, outputHeatmap), showWarnings = FALSE)
dir.create(file.path(workDir, prefixDir, outputBio), showWarnings = FALSE)
dir.create(file.path(workDir, prefixDir, outputBio, outputBioBoxplot), showWarnings = FALSE)