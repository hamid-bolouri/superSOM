# Hamid Bolouri, this version 19 Oct 2020

library(SingleCellExperiment)
library(flowCore)
library(flowWorkspace)
library(ggcyto)
library(CytoML)
library(FlowSOM)
library(CATALYST)
library(mclust)
library(MASS)
library(IDPmisc)
library(batchelor)
library(scMerge)
library(scater)
library(cowplot)
library(readxl)
library(reshape2)
library(tidyverse)
library(Hmisc)


library(BiocParallel)
MulticoreParam(workers = 16)


#-------------------------------------------------------------------------------#
#									Functions									#
#-------------------------------------------------------------------------------#

setwd("/mnt/c/Users/hbolouri/")
source("superSOM_funcs.R")

####							Function mkSCE								 ####

mkSCE <- function(gatingRes, fset, autoCD45p) {
	coFactor = 5; sce = list()
	IDs <- names(keyword(fset, "SAMPLE ID"))

	indx <- match(names(gatingRes), gsub(".fcs", "", IDs))
	fset <- fset[indx[!is.na(indx)]]
	gatingRes <- gatingRes[which(!is.na(indx))]

	for (i in 1:length(fset)) {

		print(i)
		#### NB. Uses PREDICTED CD45p cells, some with no FlowJo labels.	 ####
		selCD45p <- names(which(autoCD45p))

		shared <- intersect(selCD45p, 
							rownames(fset[[i]]@exprs))

		colnames(gatingRes[[i]]) <- rownames(fset[[i]]@exprs)
		tmpMat <- t(fset[[i]]@exprs[shared, ])

		iSel <- panel$marker_class == "type" | panel$marker_class == "state"
		selProbes <- as.character(panel$fcs_colname[iSel])
		expMat <- asinh(tmpMat[selProbes, ] / coFactor)

		flowJoGates <- t(gatingRes[[i]][ , shared])

		typesDF <- DataFrame(flowJoGates)
		colnames(typesDF) <- colnames(flowJoGates)
		rownames(typesDF) <- colnames(expMat)

		sce[[i]] <- SingleCellExperiment(assays = list(counts = expMat))
		colData(sce[[i]]) <- typesDF

	}
	names(sce) <- names(gatingRes)
	return(sce)
}


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

for (i in 1:nrow(subPanel)) {
	assign(as.character(subPanel$antigen[i]), 
	as.character(subPanel$fcs_colname[i]))
}

davidPops <- c("CD11bp.CD16p.neutrophils",
			   "CD4","CD8","CD19p","NK","monocytes",
			   "DC","eosinophil.CCR3.Siglec8","basophils")


#-------------------------------------------------------------------------------#
#								Input Query and Ref								#
#-------------------------------------------------------------------------------#

# Load reference data
sceT <- readRDS("trainingSetSCE.RDS")
# sceT is a 'SingleCellExperiment' object where rows are probes and
# each column is a previously labeled cell, with labels stored in
# the metaData of sceT.


# Load auto-gated CD45+ labels from the 'cleanUp script':
cleanupGates <- readRDS("cleanupGates.RDS")

autoCD45p <- do.call(cbind, cleanupGates)["CD45P", ]
names(autoCD45p) <- gsub("^s", "q", names(autoCD45p))

#### 					Make query SCE										 ####

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


sceAllQ <- mkSCE(gatingRes, fset, autoCD45p) 

nQuery <- length(sceAllQ)


#### 					Z-normalize 										 ####

# Z-transform:
sceAllQ <- lapply(sceAllQ, function(x) {
	assay(x, "exprs") <- t(scale(t(assay(x, "counts"))))
	return(x)
})


shared <- intersect(colnames(colData(sceAllQ[[1]])), 
					colnames(colData(sceT)))

colData(sceT) <- colData(sceT)[ , shared]


#-------------------------------------------------------------------------------#
#									Run											#
#-------------------------------------------------------------------------------#

allLabels = list()
for (iQ in 1:nQuery) { 
	lbls <- list()

	sceQ <- sceAllQ[[iQ]]
	colData(sceQ) <- colData(sceQ)[ , shared]

	# Remove Query cells from sceT:
	sceT2 <- sceT[ , setdiff(colnames(sceT), colnames(sceQ))]

	sce <- sce_cbind(list(sceT2, sceQ), exprs = "exprs", 
		   colData_names = colnames(colData(sceT2)))

	# Global clustering of major pops
	
	somX <- 35
	somY <- 35

	selProbes <- rownames(sce)
	
	print(paste0(Sys.time(), "   ", "pred MP.g Clustering started"))
	sce <- CATALYST::cluster(sce, features = selProbes, seed = 1234, 
				   xdim = somX, ydim = somY, maxK = 3)
	print(paste0(Sys.time(), "   ", "pred MP.g Clustering completed"))
	# sce metadata has SOM cluster centroids & colData have cluster_ids

	labelsFJ <- colData(sce)

	# Mask cell-type label info in Query Cells:
	qIndx <- grepl("^q", colnames(sce))
	annotCols <- !grepl("cluster_id|batch", colnames(labelsFJ))
	labelsFJ[qIndx, annotCols] <- FALSE 

	labsByClust <- split(labelsFJ, sce$cluster_id)

	labFracs <- lapply(labsByClust, function(x) {
						x <- x[ , 1:(ncol(x) - 2)]
						apply(x, 2, function(y) {
							fracGd <- sum(y) / length(y)
						})
	})

	# Fractions of all flowJo labels per SOM cluster:
	cellFracsTbl <- do.call(rbind, labFracs)

	# preMajGrps == Boolean matrix of clusters with > 80% of a major pop:
	predMajGrps <- apply(cellFracsTbl[ , davidPops], 2, function(y) y > 0.6)


	# Compensate for high miss rate of DC's:
	empties <- which(rowSums(predMajGrps) == 0) 
	fracEmpties <- cellFracsTbl[empties, colnames(predMajGrps)]

	likelyDC <- names(which(fracEmpties[ , "DC"] > 0.5)) # by def, majority

	predMajGrps[likelyDC, "DC"] <- TRUE 


	# Map predMajGrps back to a per-cell matrix comparable to colData(sce)
	predPerCell <- matrix(FALSE, ncol = ncol(predMajGrps), nrow = ncol(sce))
	colnames(predPerCell) <- colnames(predMajGrps)

	tCD14.DC = 0.75 # Determined by inspection of density of CD14 in DCs

	# Label Q cells in high-purity clusters with FloJo labels
	for (pop in colnames(predMajGrps)) {
		clusterNos <- which(predMajGrps[ , pop])
		
		iCells <- match(sce$cluster_id, clusterNos)
		iCells <- which(!is.na(iCells))
		
		# Filter-out high-CD14 cells mis-classified as DCs
		if (pop == "DC") {
			hiCD14.notDC <- which(assay(sce, "exprs")["Gd160Di", ] > tCD14.DC)
			iCells <- setdiff(iCells, hiCD14.notDC)
		}
		
		predPerCell[iCells, pop]  <- TRUE
	}

	print(paste0("Completed predMP.g  ", Sys.time()))

	# predPerCell contains predictions for all cells in sce.
	# But only those for the query sample are labeled here. 
	lbls[[1]] <- predPerCell[qIndx, ]

	# Feature-selected clustering of major pops
	lbls[[2]] <- predMP.fs(sce)

	# Grouped clustering of sub pops
	lbls[[3]] <- predSP.g(sce, lbls[[1]])

	# feature-based clustering of sub pops
	lbls[[4]] <- predSP.fs(sce, lbls[[1]])

	# Cluster CD19 sub-pops within Regions Of Interest
	lbls[[5]] <- predSP.ROI(sce, lbls[[1]])


	colnames(lbls[[1]]) <- paste0("pred.", colnames(lbls[[1]]))
	colnames(lbls[[2]]) <- paste0("pred0.", colnames(lbls[[2]]))
	colnames(lbls[[3]]) <- paste0("pred1.", colnames(lbls[[3]]))
	colnames(lbls[[4]]) <- paste0("pred2.", colnames(lbls[[4]]))
	colnames(lbls[[5]]) <- paste0("pred3.", colnames(lbls[[5]]))

	lblTbl <-  as.data.frame(do.call(cbind, lbls))


	# Post-processing refinement of sub-T pop predictions
	refinePreds.subT(sceQ) 		# updates lblTbl

	# Post-processing refinement of sub-B pop predictions	
	refinePreds.subB(sceQ) 		# updates lblTbl

	# Post-processing refinement of eosinphil & basophil  predictions
	refinePreds.EosBaso(sceQ)	# updates lblTbl

	# Post-processing refinement of pDC  predictions
	refinePreds.pDC(sce)

	# Merge predictions. Label query cells
	labelQueryCells()
	
	allLabels[[iQ]] <- lblTbl
}
names(allLabels) <- names(sceAllQ)
# LOO done!


# Save results
saveRDS(allLabels, file="allLabels_superSOM_autoClean.RDS")


# plot correlations
allPops <- colnames(allLabels[[1]])
pops2Plt <- gsub("lbl.", "", allPops[grep("^lbl", allPops)])

plts <- pltPreds(sceAllQ, allLabels, pops2Plt)


library(gridExtra)
pdf("LOO_superSOM_autoClean.pdf", 
	width=11, height = 11)
do.call(grid.arrange, plts[c(19, 12:13, 21, 18, 22, 15:17)])
do.call(grid.arrange, plts[c(9:11, 25, 20, 14, 1:4, 23:24, 5:8)])
dev.off()


plt <- pltStats(sceAllQ, allLabels, pops2Plt)

pdf("SLP19_LOO_precisionRecall_autoClean_20Oct20.pdf", 
	height = 7, width = 10)
print(plt)
dev.off()

plt <- pltStats(subSceAllQ, subAllLabels, pops2Plt)

pdf("LOO_precisionRecall_autoClean.pdf", 
	height = 7, width = 10)
print(plt)
dev.off()
























