# Hamid Bolouri, Oct 13th, 2020

library(SingleCellExperiment)
library(openCyto)
library(flowClust)
library(data.table)
library(flowCore)
library(flowWorkspace)
library(ggcyto)
library(CytoML)
library(flowDensity)
library(tidyverse)
library(mclust)
library(MASS)
library(IDPmisc)
library(Hmisc)

library(BiocParallel)
MulticoreParam(workers = 16) 


#################################################################################
#							Clean-up auto-gates									#
#################################################################################

ffiles <- dir(pattern = ".fcs$")
fset <- read.flowSet(ffiles, column.pattern = "Width", invert.pattern = TRUE)

i = 0
for (F in 1:length(fset)) {
	 i = i + 1
	 rownames(exprs(fset[[F]])) <- 
		paste0("s", i, ".c", 1:nrow(exprs(fset[[F]])))
} # Label cells uniquely to aid downstream tracking. 

# Wokspace:
wspFile <- "mySamples.wsp"  
ws <- open_flowjo_xml(wspFile)

gs <- flowjo_to_gatingset(ws, name = "Healthy Donor Whole blood samples")
sampleIDs <- sampleNames(gs)

allGates <- gs_get_pop_paths(gs, path = "auto")

gatingRes = list()
for (i in 1:length(gs)) {
	mtrx = NULL
	gh = gs[[i]]
	for (node in allGates) {
		indices <- gh_pop_get_indices(gh, node)
		mtrx <- rbind(mtrx, indices)
	}
	rownames(mtrx) <- allGates
	gatingRes[[i]] <- mtrx
}
names(gatingRes) <- sampleIDs


#-------------------------------------------------------------------------------#
#								Convenience annotations							#
#-------------------------------------------------------------------------------#

# Color Ramp for smoothScatter plots:
myRamp = colorRampPalette(c("dodgerblue", "skyblue", "white",
						    "orange", "red", "darkred"))

# Distinct colors for multi-sample plot
colors30 <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", 
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", 
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D", 
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999", 
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000", 
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

# Probe data:
setwd("/mnt/c/Users/hbolouri/_Cache/COVID19/RawData")
panel <- readRDS("probeAnnotationPanel.RDS")

# For convenience, assign antigen names to probe IDs:
subPanel <- panel[panel$marker_class != "none", ]
subPanel$antigen <- gsub("-", ".", subPanel$antigen)

for (i in 1:nrow(panel)) {
	assign(gsub("-", "_", as.character(panel$antigen[i])), 
	as.character(panel$fcs_colname[i]))
}

davidPops <- c("CD11bp.CD16p.neutrophils",
			   "CD4","CD8","CD19p","NK","monocytes",
			   "DC","eosinophil.CCR3.Siglec8","basophils")


#-------------------------------------------------------------------------------#
#							Run clean-up gates									#
#-------------------------------------------------------------------------------#

# Assignments for validation set:
slpFset = fset
gatingResSLP = gatingRes

cleanupGates = list()
for (i in 1:length(slpFset)) {
	logExp <- apply(exprs(slpFset[[i]]), 2, 
					function(y) y <- log10(y + 1))

	blankGate <- rep(FALSE, nrow(logExp))
	names(blankGate) <- rownames(logExp)
	Singlets = Live = CD45P = blankGate


	# Gate 1 - Cells:
	resMclust <- Mclust(logExp[ , EQbeads], G = 2)

	sizes <- table(resMclust$"classification")
	selClust <- ifelse(sizes[1] > sizes[2], 1, 2)

	Cells <- resMclust$"classification" == selClust
	logExp <- logExp[Cells, ]


	# Gate 2 - Singlets:
	D <- density(logExp[ , DNA1])
	indx <- which.max(D$y)
	peakDNA1 <- D$x[indx]
	sdDNA1 <- sd(logExp[ (logExp[ , DNA1] > 0.9 * peakDNA1 & 
						  logExp[ , DNA1] < 1.1 * peakDNA1), DNA1 ])

	loDNA1 <- peakDNA1 - 3.291 * sdDNA1 # 99.9% CI
	hiDNA1 <- peakDNA1 + 2.807 * sdDNA1 # 99.5% CI

	inGate <- logExp[ , DNA1] > loDNA1 &
			  logExp[ , DNA1] < hiDNA1

	D <- density(logExp[inGate, Event_length])
	M <- peaks(D)$x[which.max(peaks(D)$y)]
	nearby <- logExp[ , Event_length] > 0.8 * M &
			  logExp[ , Event_length] < 1.2 * M
	SD <- sd(logExp[nearby, Event_length])

	threshEL <- M + 2.807 * SD

	gate2 <- names(which( 
			 logExp[ , Event_length] < threshEL &
			 logExp[ , DNA1] > loDNA1 &
			 logExp[ , DNA1] < hiDNA1 
			 ))

	Singlets[gate2] <- TRUE
	logExp <- logExp[gate2, ]


	# Gate 3 - Live:

	# David uses Siglec8-hi eosinphils to set hiCIS
	q99 <- quantile((logExp[ , Siglec8]), prob=0.99)

	# Method 1: 
	D <- density(logExp[logExp[ , Siglec8] > q99 , Cisplatin_198])
	
	peakCIS <- max(peaks(D$x, D$y))

	# Standard deviation of points immediately around the peak:
	sdCIS <- sd(logExp[ (logExp[ , Cisplatin_198] > 0.9 * peakCIS & 
						 logExp[ , Cisplatin_198] < 1.1 * peakCIS), 
						 Cisplatin_198 ])

	ciCIS <- peakCIS + 2.576 * sdCIS # 95% CI

	# Method 2:
	topS8 <- logExp[logExp[ , Siglec8] > q99 , Cisplatin_198]
	mClustRes <- Mclust(topS8, G = 2)
	clust2 <- names(mClustRes$classification)[mClustRes$classification == 2]

	# Combine (see 'Gate 3 diagnostic plot', below):
	hiCIS <- max(ciCIS, max(logExp[clust2, Cisplatin_198]))

	gate3 <- names(which( logExp[ , Cisplatin_198] < hiCIS ))

	Live[gate3] <- TRUE
	logExp <- logExp[gate3, ]


	# Gate 4 - CD45+:
	# Using kmeans for CD16 because of multiple pops in some samples 
	kmClustCD16 <- kmeans(logExp[ , CD16], centers=2, nstart = 20)
	
	loCD16 <- ifelse(mean(logExp[kmClustCD16$cluster == 1, CD16]) <
					 mean(logExp[kmClustCD16$cluster == 2, CD16]), 1, 2)

	mClustCD45 <- Mclust(logExp[kmClustCD16$cluster == loCD16, CD45], G = 2)

	threshCD16 <- max(logExp[(kmClustCD16$cluster == loCD16), CD16])

	D <- density(logExp[(kmClustCD16$cluster == loCD16), CD45], from = 0)
	troughCD45 <- min(peaks(D$x, (max(D$y) - D$y))$x)

	gate4 <- names(which( logExp[ , CD16] > threshCD16 |
						  logExp[ , CD45] > troughCD45 ))

	CD45P[gate4] <- TRUE

	cleanupGates[[i]] <- rbind(Cells, Singlets, Live, CD45P)

	print(paste0(i, "   ", Sys.time()))
}
names(cleanupGates) <- names(gatingResSLP)


saveRDS(cleanupGates, file="cleanupGates.RDS")


pops <- rownames(cleanupGates[[1]])
plts <- pltPreds(gatingResSLP, cleanupGates, pops)

library(gridExtra)

do.call(grid.arrange, plts)


plts <- pltStats(gatingResSLP, cleanupGates, pops)

plts


