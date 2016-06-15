#!/usr/bin/R
# The great script to run
# remove.packages("BiocInstaller")
# source("http://bioconductor.org/biocLite.R")
# biocLite()
rm(list=ls())
wds <- "/home/anton/data/R-script/R-counts/RUN"

# Load libraries, variables and create folders
source(file.path(wds, "variables.R"))

# Whether a script run earlier?
###############################
setwd(file.path(workDir, prefixDir))
	alreadyRun <- list.files(prefixDir, "^.*Step_01.*")
	if (length(alreadyRun) == 1) {
		startCol <- 0
	}
setwd(workDir)

# Load functions
source(file.path(wds, "functions.R"))

# Run intial prepare of DamID data
source(file.path(wds, "damID_prepare.R"))

# Run bioHMM function and domain interpretations
source(file.path(wds, "domains.R"))

# Run biological search
source(file.path(wds, "meaning.R"))

# Run calculate on expression data
source(file.path(wds, "expression.R"))

# load garbadge code
# source(file.path(wds, "gccode.R"))

# Run testing code
source(file.path(wds, "testing.R"))