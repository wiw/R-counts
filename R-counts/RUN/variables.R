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
# output directory into working directory
prefixDir <- paste("RUN", format(Sys.time(), "%d-%m-%Y"), "ludo", sep="_") 
# use only edge reads to counts or not
onlyEdge <- F
# working directory (WD)
workDir <- "/home/anton/data/R-script/R-counts"	
# location your RData files. You can specify the highest folder as it is possible. Searching runs recursively.
sourceDir <- "/home/anton/backup/ludo-output" 
# location your DamID-Description file
damIdLocation <- "/home/anton/data/DAM/RUN/damid_description.csv" 
# location to your Genes file
genesFilePath <- "/home/anton/backup/2013-07-09_Assigning_GATCs_to_Genes_(FBgn_IDs)/Drosophila_melanogaster.BDGP5.25.64_Genes_AP130708.txt"
# location to your housekeeping genes file
HKGenesPath <- "/home/anton/backup/2013-07-09_Assigning_GATCs_to_Genes_(FBgn_IDs)/houskeeping_genes.csv"
# folder with your expression data files
exprDataDir <- "/home/anton/backup/input/ExpressionData"
# output folders in work directory
outputGff <- "gff"
outputWig <- "wig"
outputScttr <- "scatter_plots"
outputDomain <- "domains"
outputCleanStat <- "clean_stat"
outputBio <- "Bio"
outputHeatmap <- "Heatmap"
outputBioBoxplot <- "Boxplot"
outputExpr <- "Expression"
# the number of last column in GATCs file, default "7"
startCol <- 7
# location you GATCs file
gatcFile <- paste(workDir, "GATCs_mod.txt", sep="/")
gatcFile4Genes <- paste(workDir, "GATCs.txt", sep="/")
# are you need to combine some columns into one (T or F)? we recommend set the "F"; if you select "T" - edit the "combine" vector on string #85 into this code
needCombine <- F
# are you need to add pseudo counts into source data (T or F)? Default "T"
usePseudoCounts <- F
# the vector of pseudo counts, default "c(1)"
pseudoCounts <- c(0.01)
# the vector of type correlation method - please don't change this setting
corrMethod <- c("pearson", "spearman")	
# check you color into heatmap
heatmapColors <- greenred(200)	
# the vector of heatmap label - please don't change this setting
labelHeatmap <- c("B", "C")	
# the vector of acf label - please don't change this setting
labelAcf <- c("A_ALL", "B_MA")	
# use this option if you need to get the intermediate files, default "T"
writeTemp <- T		
# if need calculate only some files from source list, default "F"
needSomeFiles <- F	
# the region of files which need calculate
useSomeFiles <- c(9:16)	
# Part of tissies from all
tissue_bio_set <- c("BR", "FB", "Glia", "NRN", "Kc167") 
# Part of protein from all
protein_bio_set <- c("LAM", "HP1", "PC") 
# Part of conditions from all
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