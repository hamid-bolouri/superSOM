# Hamid Bolouri, 6 Dec 2020


#-------------------------------------------------------------------------------#
#				Convenience function to de-clutter code 						#
#-------------------------------------------------------------------------------#

getPredictions <- function(subSCE, popLbls, purity, pops2ClustL, iGrp) {

	labsByClust <- split(popLbls, subSCE$cluster_id)

	labFracs <- lapply(labsByClust, function(x) {
					   x <- as_tibble(x)
					   x <- x[ , !grepl("^pred|cluster", colnames(x))]
					   
					   if (nrow(x) > 0) {
						  fracGd <- colSums(x) / nrow(x)
						} else fracGd <- rep(0, ncol(x))
					   
					   return(fracGd)
					   })

	# Fractions of all flowJo labels per SOM cluster:
	cellFracsTbl <- do.call(rbind, labFracs)

	# predGrps == Boolean matrix of clusters with > x% of a pop:
	predGrps <- apply(cellFracsTbl, 2, function(y) y > purity)

	colnames(predGrps) <- pops2ClustL[[iGrp]]$qCells

	# Map predGrps back to a per-cell matrix comparable to colData(subSCE)
	predPerCell <- matrix(FALSE, ncol = ncol(predGrps), nrow = ncol(subSCE))
	colnames(predPerCell) <- colnames(predGrps)
	rownames(predPerCell) <- colnames(subSCE)
	
	for (pop in colnames(predGrps)) {
		clusterNos <- as.character(which(predGrps[ , pop]))
		clustersInSCE <- as.character(subSCE$cluster_id)
		
		iCells <- clustersInSCE %in% clusterNos
		
		predPerCell[iCells, pop]  <- TRUE
	}
	return(predPerCell)
}


#-------------------------------------------------------------------------------#
#				Featured-selected prediction of major populations				#
#-------------------------------------------------------------------------------#

predMP.fs <- function(sce) {

	# Slow, so labeling only selected cell pops based on SLP data
	pops2Clust0 <- list(
		CD4.8 = list(
			qCells = c("CD4", "CD8"),
			probes = c("CCR3", "CD56", "CD16", "CD3", "CD19", "CD4", "CD8a"),
			somXY = 30
			),
		B = list(
			qCells = "CD19p",
			probes = c("CCR3", "CD56", "CD16", "CD3", "CD19"),
			somXY = 30
			),
		EOS = list(
			qCells = "eosinophil.CCR3.Siglec8",
			probes = c("CCR3", "CD123", "FceRI", "Siglec8"),
			somXY = 30
			),
		Baso = list(
			qCells = "basophils",
			probes = c("CCR3", "CD123", "FceRI"),
			somXY = 30
			),
		Neut = list(
			qCells = "CD11bp.CD16p.neutrophils",
			probes = c("CCR3", "CD56", "CD16", "CD66b", "CD11b"),
			somXY = 30
			),
		MonoDC = list(
			qCells = c("monocytes", "DC"),
			probes = c("CCR3", "CD56", "CD16", "CD3", "CD19", "CD14", "HLA.DR"),
			somXY = 30
			),
		NK = list(
			qCells = "NK",
			probes = c("CCR3", "CD56", "CD16", "CD3", "CD19", "CD14", "HLA.DR", "CD66b"),
			somXY = 30
			)
	)

	# Co-cluster Ref and Q

	purity = 0.6 # min fraction of same-type cells in a SOM cluster

	cells2Lbl <- unlist(lapply(pops2Clust0, function(x) x$qCells))
	pred0Tbl <- matrix(FALSE, ncol = length(cells2Lbl), 
					  nrow = ncol(sce))
	colnames(pred0Tbl) <- cells2Lbl
	rownames(pred0Tbl) <- colnames(sce)

	for (iGrp in 1:length(pops2Clust0)) {
		
		iProbes <- match(pops2Clust0[[iGrp]]$probes, subPanel$antigen)
		selProbes <- subPanel$fcs_colname[iProbes]
		
		subSCE <- sce[selProbes, ] # All CD45+ cells
		
		popLbls <- tibble(as.data.frame(colData(subSCE)[ , pops2Clust0[[iGrp]]$qCells]))
		
		# Remove cell-type labels of Query cells:
		rownames(popLbls) <- colnames(sce)
		popLbls[colnames(sceQ), ] <- FALSE
		
		# Run co-clustering SOM
		print(paste0(iGrp, " MP.fs started.   ", Sys.time()))	
		subSCE <- CATALYST::cluster(subSCE, features = selProbes, seed = 1234,
				  xdim = pops2Clust0[[iGrp]]$somXY, 
				  ydim = pops2Clust0[[iGrp]]$somXY, maxK = 3)
		
		print(paste0(iGrp, " MP.fs done.   ", Sys.time()))	

		predRes <- getPredictions(subSCE, popLbls, purity, pops2Clust0, iGrp)		
		
		pred0Tbl[colnames(sceQ), colnames(predRes)] <- predRes[colnames(sceQ), ]
	}
	return(pred0Tbl[colnames(sceQ), ])
} 


#-------------------------------------------------------------------------------#
#					SOM grouped-clustering of sub-populations					#
#-------------------------------------------------------------------------------#

predSP.g <- function(sce, lbls1) {

	# define sub-pops and associated informative features 
	pops2Clust1 <- list(
		CD4.1 = list(
			qCells = c("CD4.naive", "CD4.CM", "CD4.EM", "CD4.TEMRA"),
			parentPop = "CD4",
			probes = c("CD45RA", "CCR7"),
			somXY = 30
			),
		CD4.2 = list(
			qCells = "Treg",
			parentPop = "CD4",
			probes = c("CD127", "CD25"),
			somXY = 30
			),
		CD4.3 = list(
			qCells = "CXCR5p",
			probes = c("CD45RA", "CCR7", "CXCR5", "PD.1"),
			parentPop = "CD4",
			somXY = 30
			),
		CD8 = list(
			qCells = c("CD8.naive", "CD8.CM", "CD8.EM", "CD8.TEMRA"),
			parentPop = "CD8",
			probes = c("CCR7", "CD45RA"),
			somXY = 25
			),
		B = list(
			qCells = c("naive.B.cells", "pre.sw.memory.B.cells", 
					   "post.sw.memory.B.cells"),
			parentPop = "CD19p",
			probes = c("CD27", "IgD"),
			somXY = 20
			),
		plas = list(
			qCells = c("plasmablasts"),
			parentPop = "CD19p",
			probes = c("CD20", "CD27", "CD38"),
			somXY = 20
			),
		subDC = list(
			qCells = c("cDC", "pDC"),
			parentPop = "DC",
			probes = c("CD123", "CD11c", "HLA.DR", "CD14"),
			somXY = 20
			)
	)

	# Co-cluster Ref and Q. Label Q.

	purity = 0.6 # min fraction of same-type cells in a SOM cluster

	cells2Lbl <- unlist(lapply(pops2Clust1, function(x) x$qCells))
	predTbl <- matrix(FALSE, ncol = length(cells2Lbl), 
					  nrow = ncol(sce))
	colnames(predTbl) <- cells2Lbl
	rownames(predTbl) <- colnames(sce)

	for (iGrp in 1:length(pops2Clust1)) {
		
		# Select parent cells from FlowJo for Ref and from predictions for Q:
		parentID <- pops2Clust1[[iGrp]]$parentPop	
		
		# parentCells == TRUE if cell is of type parentCell:
		parentCells <- colData(sce)[ , parentID]
		names(parentCells) <- colnames(sce)
		
		# Remove flowJo labels from query parentCells:
		parentCells[colnames(sceQ)] <- FALSE 
		
		# Mark Query parentCells with predicted parentCell values:
		parentCells[names(which(lbls1[ , parentID]))] <- TRUE
		
		iProbes <- match(pops2Clust1[[iGrp]]$probes, subPanel$antigen)
		subProbes <- subPanel$fcs_colname[iProbes]
		
		subSCE <- sce[subProbes, parentCells]
		
		popLbls <- tibble(as.data.frame(colData(subSCE)[ , pops2Clust1[[iGrp]]$qCells]))
		rownames(popLbls) <- colnames(subSCE)
		
		# Remove cell-type labels of Query cells:
		popLbls[intersect(colnames(sceQ), colnames(subSCE)), ] <- FALSE
		
		# Run co-clustering SOM
		subSCE <- CATALYST::cluster(subSCE, features = subProbes, seed = 1234,
				  xdim = pops2Clust1[[iGrp]]$somXY, 
				  ydim = pops2Clust1[[iGrp]]$somXY, maxK = 3)
		
		predRes <- getPredictions(subSCE, popLbls, purity, pops2Clust1, iGrp)
		
		# Map result indices from subSCE back to sce: 
		qCellsInRes <- intersect(colnames(sceQ), rownames(predRes))
		predTbl[qCellsInRes, colnames(predRes)] <- predRes[qCellsInRes, ]
		
		print(paste0(iGrp, " Pred1.   ", Sys.time()))
	}

	return(predTbl[colnames(sceQ), ])
}


#-------------------------------------------------------------------------------#
#			Feature-selection per cell-type SOM sub-type clustering				#
#-------------------------------------------------------------------------------#

predSP.fs <- function(sce, lbls1) {

	####					Using PER CELL features. 						 ####

	pops2Clust2 <- list(
		CD8.EMRA = list( 
			qCells = "CD8.TEMRA",
			parentPop = "CD8",
			probes = c("CCR7", "CD45RA"),
			somXY = 25
			),
		plas = list(
			qCells = c("plasmablasts"),
			parentPop = "CD19p",
			probes = c("CD20", "CD27", "CD38"),
			somXY = 25
			),
		cDC = list(
			qCells = "cDC",
			parentPop = "DC",
			probes = c("CD11c", "CD123"),
			somXY = 25
			),
		pDC = list(
			qCells = "pDC",
			parentPop = "DC",
			probes = c("CD11c", "CD123"),
			somXY = 25
			)
	)

	####				Run SOM on defined cells & probes					 ####

	purity = 0.6 # min fraction of same-type cells in a SOM cluster

	cells2Lbl <- unlist(lapply(pops2Clust2, function(x) x$qCells))
	predTbl2 <- matrix(FALSE, ncol = length(cells2Lbl), 
					  nrow = ncol(sce))
	colnames(predTbl2) <- cells2Lbl
	rownames(predTbl2) <- colnames(sce)


	for (iGrp in 1:length(pops2Clust2)) {
	
		# Select parent cells from FlowJo for Ref and from predictions for Q:
		parentID <- pops2Clust2[[iGrp]]$parentPop	
		
		# parentCells == TRUE if cell is of type parentCell:
		parentCells <- colData(sce)[ , parentID]
		names(parentCells) <- colnames(sce)
		
		# Remove flowJo labels from query parentCells:
		parentCells[colnames(sceQ)] <- FALSE 
		
		# Mark Query parentCells with predicted parentCell values:
		parentCells[names(which(lbls1[ , parentID]))] <- TRUE
		
		iProbes <- match(pops2Clust2[[iGrp]]$probes, subPanel$antigen)
		subProbes <- subPanel$fcs_colname[iProbes]
		
		subSCE <- sce[subProbes, parentCells]
		
		popLbls <- tibble(as.data.frame(colData(subSCE)[ , pops2Clust2[[iGrp]]$qCells]))
		rownames(popLbls) <- colnames(subSCE)
		
		# Remove cell-type labels of Query cells:
		popLbls[intersect(colnames(sceQ), colnames(subSCE)), ] <- FALSE
		
		# Run co-clustering SOM
		subSCE <- CATALYST::cluster(subSCE, features = subProbes, seed = 1234,
				  xdim = pops2Clust2[[iGrp]]$somXY, 
				  ydim = pops2Clust2[[iGrp]]$somXY, maxK = 3)
		
		predRes <- getPredictions(subSCE, popLbls, purity, pops2Clust2, iGrp)
		
		# Map result indices from subSCE back to sce: 
		qCellsInRes <- intersect(colnames(sceQ), rownames(predRes))
		predTbl2[qCellsInRes, colnames(predRes)] <- predRes[qCellsInRes, ]
		
		print(paste0(iGrp, " Pred0.   ", Sys.time()))
	}
	return(predTbl2[colnames(sceQ), ])
}


#-------------------------------------------------------------------------------#
#					ROI & feature-selected SOM sub-type clustering				#
#-------------------------------------------------------------------------------#

# Use training set to define a Region Of Interest (plausible range) per probe.
predSP.ROI <- function(sce, lbls1) {

	####		Define per cell features. Per cell ROIs are added below.	 ####
	
	pops2Clust3 <- list(
		B.naive = list(
			qCells = "naive.B.cells", 
			parentPop = "CD19p",
			probes = c("CD27", "IgD"),
			somXY = 25
			),
		B.pre = list(
			qCells = "pre.sw.memory.B.cells",
			parentPop = "CD19p",
			probes = c("CD27", "IgD"),
			somXY = 25
			),
		B.post = list(
			qCells = "post.sw.memory.B.cells",
			parentPop = "CD19p",
			probes = c("CD27", "IgD"),
			somXY = 25
			)
	)

	#### 			Define ROIs within each of the above pops				 ####

	mat = assay(sce, "exprs")[c(IgD, CD27), ]

	# exclude Query from ROI definition:
	refNaiveB <- sce$naive.B.cells[setdiff(colnames(sce), colnames(sceQ))]

	tIgD <- min(mat[IgD, refNaiveB])
	tCD27 <- max(mat[CD27, refNaiveB])

	ROI.Bnv <- names(which( (mat[IgD, ] > 0.9 * tIgD) &
							(mat[CD27, ] < 1.1 * tCD27) ))

	ROI.Bpre <- names(which( (mat[IgD, ] > 0.9 * tIgD) &
							 (mat[CD27, ] > 0.9 * tCD27) ))

	ROI.Bpost <- names(which( (mat[IgD, ] < 1.1 * tIgD) &
							 (mat[CD27, ] > 0.9 * tCD27) ))

	pops2Clust3[["B.naive"]]$ROI <- ROI.Bnv
	pops2Clust3[["B.pre"]]$ROI <- ROI.Bpre
	pops2Clust3[["B.post"]]$ROI <- ROI.Bpost

	####					Run SOM on defined cells & ROI					 ####

	purity = 0.6 # min fraction of same-type cells in a SOM cluster

	cells2Lbl <- unlist(lapply(pops2Clust3, function(x) x$qCells))
	predTbl3 <- matrix(FALSE, ncol = length(cells2Lbl), 
					  nrow = ncol(sce))
	colnames(predTbl3) <- cells2Lbl
	rownames(predTbl3) <- colnames(sce)


	for (iGrp in 1:length(pops2Clust3)) {
		
		# Select parent cells from FlowJo for Ref and from predictions for Q:
		parentID <- pops2Clust3[[iGrp]]$parentPop	
		
		# parentCells == TRUE if cell is of type parentCell:
		parentCells <- colData(sce)[ , parentID]
		names(parentCells) <- colnames(sce)
		
		# Remove flowJo labels from query parentCells:
		parentCells[colnames(sceQ)] <- FALSE 
		
		# Mark Query parentCells with predicted parentCell values:
		parentCells[names(which(lbls1[ , parentID]))] <- TRUE
		
		iProbes <- match(pops2Clust3[[iGrp]]$probes, subPanel$antigen)
		subProbes <- subPanel$fcs_colname[iProbes]
		
		subSCE <- sce[subProbes, parentCells]
		
		popLbls <- tibble(as.data.frame(colData(subSCE)[ , pops2Clust3[[iGrp]]$qCells]))
		rownames(popLbls) <- colnames(subSCE)
		
		# Remove cell-type labels of Query cells:
		popLbls[intersect(colnames(sceQ), colnames(subSCE)), ] <- FALSE
		
		# Run co-clustering SOM
		subSCE <- CATALYST::cluster(subSCE, features = subProbes, seed = 1234,
				  xdim = pops2Clust3[[iGrp]]$somXY, 
				  ydim = pops2Clust3[[iGrp]]$somXY, maxK = 3)
		
		predRes <- getPredictions(subSCE, popLbls, purity, pops2Clust3, iGrp)
		
		# Map result indices from subSCE back to sce: 
		qCellsInRes <- intersect(colnames(sceQ), rownames(predRes))
		predTbl3[qCellsInRes, colnames(predRes)] <- predRes[qCellsInRes, ]
		
		print(paste0(iGrp, " Pred3.   ", Sys.time()))
	}

	return(predTbl3[colnames(sceQ), ])
}


#-------------------------------------------------------------------------------#
#								Post-processing									#
#-------------------------------------------------------------------------------#

####			 FlowJo T subset gates are rectangular						 ####
####			 Get rectangular SOM bounds for subTs.						 ####

getThreshT <- function(sceQ) {

	mat <- assay(sceQ, "exprs") # Query intensity matrix

	#### CD45RA boundary threshold

	# CD8 cells
	M1 <- max(mat[CCR7, lblTbl$pred1.CD8.TEMRA])
	M2 <- max(mat[CCR7, lblTbl$pred1.CD8.EM])
	M <- min(M1, M2)

	selCells <- mat[CCR7, lblTbl$pred.CD8] < M

	resMclust <- Mclust(mat[CD45RA, selCells], G = 2)

	centroid1 <- 
		mean(mat[CD45RA, selCells][resMclust$classification == 1])
	centroid2 <- 
		mean(mat[CD45RA, selCells][resMclust$classification == 2])
	N  <- ifelse(centroid1 < centroid2, 1, 2)

	M1 <- quantile(mat[CD45RA, selCells][resMclust$classification 
				   == N], probs= 0.99)
	M2 <- quantile(mat[CD45RA, selCells][resMclust$classification 
				   == setdiff(1:2, N)], probs=0.01)

	tCD45RA <- mean(M1, M2)


	# CD4 cells
	M1 <- max(mat[CCR7, lblTbl$pred1.CD4.TEMRA])
	M2 <- max(mat[CCR7, lblTbl$pred1.CD4.EM])
	M <- min(M1, M2)

	selCells <- mat[CCR7, lblTbl$pred.CD4] < M

	resMclust <- Mclust(mat[CD45RA, selCells], G = 2)

	centroid1 <- 
		mean(mat[CD45RA, selCells][resMclust$classification == 1])
	centroid2 <- 
		mean(mat[CD45RA, selCells][resMclust$classification == 2])
	N  <- ifelse(centroid1 < centroid2, 1, 2)

	M1 <- quantile(mat[CD45RA, selCells][resMclust$classification 
				   == N], probs= 0.99)
	M2 <- quantile(mat[CD45RA, selCells][resMclust$classification 
				   == setdiff(1:2, N)], probs=0.01)

	tCD45RA <- mean( c(tCD45RA, mean(M1, M2)) )


	#### CCR7 boundary threshold

	# CD8 cells
	M1 <- max(mat[CD45RA, lblTbl$pred1.CD8.CM])
	M2 <- max(mat[CD45RA, lblTbl$pred1.CD8.EM])
	M <- max(M1, M2)
	
	selCells <- mat[CD45RA, lblTbl$pred.CD8] < M

	resMclust <- Mclust(mat[CCR7, selCells], G = 2)

	centroid1 <- 
		mean(mat[CCR7, selCells][resMclust$classification == 1])
	centroid2 <- 
		mean(mat[CCR7, selCells][resMclust$classification == 2])
	N  <- ifelse(centroid1 < centroid2, 1, 2)

	M1 <- quantile(mat[CCR7, selCells][resMclust$classification 
				   == N], probs= 0.99)
	M2 <- quantile(mat[CCR7, selCells][resMclust$classification 
				   == setdiff(1:2, N)], probs=0.01)

	tCCR7 <- mean(M1, M2)

	M1 <- min(mat[CD45RA, lblTbl$pred1.CD8.naive])
	M2 <- min(mat[CD45RA, lblTbl$pred1.CD8.TEMRA])
	M <- min(M1, M2)

	selCells <- mat[CD45RA, lblTbl$pred.CD8] > M

	resMclust <- Mclust(mat[CCR7, selCells], G = 2)

	centroid1 <- 
		mean(mat[CCR7, selCells][resMclust$classification == 1])
	centroid2 <- 
		mean(mat[CCR7, selCells][resMclust$classification == 2])
	N  <- ifelse(centroid1 < centroid2, 1, 2)

	M1 <- quantile(mat[CCR7, selCells][resMclust$classification 
				   == N], probs= 0.99)
	M2 <- quantile(mat[CCR7, selCells][resMclust$classification 
				   == setdiff(1:2, N)], probs=0.01)

	tCCR7 <- mean( c(tCCR7, mean(M1, M2)) )

	# CD4 cells
	M1 <- max(mat[CD45RA, lblTbl$pred1.CD4.CM])
	M2 <- max(mat[CD45RA, lblTbl$pred1.CD4.EM])
	M <- max(M1, M2)

	selCells <- mat[CD45RA, lblTbl$pred.CD4] < M

	resMclust <- Mclust(mat[CCR7, selCells], G = 2)

	centroid1 <- 
		mean(mat[CCR7, selCells][resMclust$classification == 1])
	centroid2 <- 
		mean(mat[CCR7, selCells][resMclust$classification == 2])
	N  <- ifelse(centroid1 < centroid2, 1, 2)

	M1 <- quantile(mat[CCR7, selCells][resMclust$classification 
				   == N], probs= 0.99)
	M2 <- quantile(mat[CCR7, selCells][resMclust$classification 
				   == setdiff(1:2, N)], probs=0.01)

	tCCR7 <- c(tCCR7, mean(M1, M2))

	M1 <- min(mat[CD45RA, lblTbl$pred1.CD4.naive])
	M2 <- min(mat[CD45RA, lblTbl$pred1.CD4.TEMRA])
	M <- min(M1, M2)

	selCells <- mat[CD45RA, lblTbl$pred.CD4] > M

	resMclust <- Mclust(mat[CCR7, selCells], G = 2)

	centroid1 <- 
		mean(mat[CCR7, selCells][resMclust$classification == 1])
	centroid2 <- 
		mean(mat[CCR7, selCells][resMclust$classification == 2])
	N  <- ifelse(centroid1 < centroid2, 1, 2)

	M1 <- quantile(mat[CCR7, selCells][resMclust$classification 
				   == N], probs= 0.99)
	M2 <- quantile(mat[CCR7, selCells][resMclust$classification 
				   == setdiff(1:2, N)], probs=0.01)

	tCCR7 <- mean( c(tCCR7, mean(M1, M2)) )

	return( c(tCD45RA = tCD45RA, tCCR7 = tCCR7) )
}


####				Post-process CD4/CD8 subset predictions					 ####

refinePreds.subT <- function(sceQ) {

	# Stretch/trim sub-CD4/8 SOM clusters into rectangular regions.
	# Boundary detection uses a GMM with K=2 in SOM-defined areas.
	# Threshold previously calculated in file: "sqSubTs.R"

	Ts <- getThreshT(sceQ)
	thresholdCCR7 <- Ts["tCCR7"]
	thresholdCD45RA <- Ts["tCD45RA"]

	# Only post-process subsets that show benefit:
	lblTbl$lbl.CD4.naive <<- FALSE
	lblTbl$lbl.CD4.CM <<- lblTbl$pred1.CD4.CM
	lblTbl$lbl.CD4.EM <<- lblTbl$pred1.CD4.EM
	lblTbl$lbl.CD4.TEMRA <<- lblTbl$pred1.CD4.TEMRA
	lblTbl$lbl.CD8.naive <<- lblTbl$pred1.CD8.naive
	lblTbl$lbl.CD8.CM <<- FALSE
	lblTbl$lbl.CD8.EM <<- lblTbl$pred1.CD8.EM
	lblTbl$lbl.CD8.TEMRA <<- lblTbl$pred1.CD8.TEMRA

	iCD4 <- lblTbl[ , "pred.CD4"]
	iCD8 <- lblTbl[ , "pred.CD8"]

	hiCCR7 <- assay(sceQ, "exprs")[CCR7, ] > 
			(1.1 * thresholdCCR7)
	loCCR7 <- assay(sceQ, "exprs")[CCR7, ] < 
			(0.9 * thresholdCCR7)

	hi45RA <- assay(sceQ, "exprs")[CD45RA, ] > 
			(1.1 * thresholdCD45RA)
	lo45RA <- assay(sceQ, "exprs")[CD45RA, ] < 
			(0.9 * thresholdCD45RA)

	# CD4.naive	
	i <- iCD4 & hiCCR7 & hi45RA
	lblTbl$lbl.CD4.naive[i] <<- TRUE

	#### Below, unused predictions are commented out:
	
	# # CD4.CM
	# i <- iCD4 & hiCCR7 & lo45RA
	# lblTbl$lbl.CD4.CM[i] <<- TRUE

	# # CD4.EM	
	# i <- iCD4 & loCCR7 & lo45RA
	# lblTbl$lbl.CD4.EM[i] <<- TRUE

	# # CD4.EMRA
	# i <- iCD4 & loCCR7 & hi45RA
	# lblTbl$lbl.CD4.TEMRA[i] <<- TRUE


	# # CD8.naive	
	# i <- iCD8 & hiCCR7 & hi45RA
	# lblTbl$lbl.CD8.naive[i] <<- TRUE

	# CD8.CM
	i <- iCD8 & hiCCR7 & lo45RA
	lblTbl$lbl.CD8.CM[i] <<- TRUE

	# # CD8.EM	
	# i <- iCD8 & loCCR7 & lo45RA
	# lblTbl$lbl.CD8.EM[i] <<- TRUE

	# # CD8.EMRA
	# i <- iCD8 & loCCR7 & hi45RA
	# lblTbl$lbl.CD8.TEMRA[i] <<- TRUE
	
	return()
} # NB. no return value. Global lblTbl update.
 

####	Use SOM clusters to set rectangular B cell subset probe thresholds	 ####

getThreshB <- function(sceQ) {

	# Note: Using ROI SOMs (pred3)
	mat <- assay(sceQ, "exprs") # Query intensity matrix

	postSW <- quantile(
		mat[CD27, lblTbl$pred1.post.sw.memory.B.cells | 
				  lblTbl$pred3.post.sw.memory.B.cells ],
		probs = 0.1)
	preSW <- quantile(
		mat[CD27, lblTbl$pred1.pre.sw.memory.B.cells | 
				  lblTbl$pred3.pre.sw.memory.B.cells ],
		probs = 0.1)

	RHS <- min(preSW, postSW)
	LHS <- quantile(mat[CD27, lblTbl$pred1.naive.B.cells | 
							  lblTbl$pred3.naive.B.cells ],
		probs = 0.9)

	thresholdCD27 <- mean( c(LHS + RHS) )


	naiveB <- quantile(
		mat[IgD, lblTbl$pred1.naive.B.cells | 
				 lblTbl$pred3.naive.B.cells ],
		probs = 0.1)
	preSW <- quantile(
		mat[IgD, lblTbl$pred1.pre.sw.memory.B.cells | 
				 lblTbl$pred3.pre.sw.memory.B.cells],
		probs = 0.1)

	top <- min(naiveB, preSW)
	bottom <- quantile(
		mat[IgD, lblTbl$pred1.post.sw.memory.B.cells | 
				 lblTbl$pred3.post.sw.memory.B.cells ],
		probs = 0.9)

	thresholdIgD <- mean( c(top + bottom) )

	return( c(tCD27 = thresholdCD27, tIgD = thresholdIgD) )
}


####	 Convert B cell subset SOMs to rectangular gates					 ####

refinePreds.subB <- function(sceQ) {

	# Set rectangular boundaries tessalating around SOM clusters.

	Ts <- getThreshB(sceQ)
	thresholdCD27 <- Ts["tCD27"]
	thresholdIgD <- Ts["tIgD"]

	# pre & post sw counts did not improve significantly. So use SOM values:
	lblTbl$lbl.pre.sw.memory.B.cells <<- lblTbl$pred1.pre.sw.memory.B.cells
	lblTbl$lbl.post.sw.memory.B.cells <<- lblTbl$pred1.post.sw.memory.B.cells
	lblTbl$lbl.naive.B.cells <<- FALSE

	iB <- lblTbl[ , "pred.CD19p"]

	hiCD27 <- assay(sceQ, "exprs")[CD27, ] > 
			(1.1 * thresholdCD27)
	loCD27 <- assay(sceQ, "exprs")[CD27, ] < 
			(0.9 * thresholdCD27)

	hiIgD <- assay(sceQ, "exprs")[IgD, ] > 
			(1.1 * thresholdIgD)
	loIgD <- assay(sceQ, "exprs")[IgD, ] < 
			(0.9 * thresholdIgD)

	# # B.preSW
	# i <- iB & hiCD27 & hiIgD
	# lblTbl$lbl.pre.sw.memory.B.cells[i] <<- TRUE

	# # B.postSW	
	# i <- iB & hiCD27 & loIgD
	# lblTbl$lbl.post.sw.memory.B.cells[i] <<- TRUE

	# B.naive
	i <- iB & loCD27 & hiIgD
	lblTbl$lbl.naive.B.cells[i] <<- TRUE
	
	return()
	
} # Note lblTbl is revised globally with "<<-"


####				Post-process eosinophil/basophil predictions			 ####

# Get cells expressing high levels of 'probe'
getHi <- function(probe, thresh, fracTot = TRUE) {

	if (fracTot) {
		H = hist(assay(sceQ, "exprs")[probe, ], breaks = 100, plot = FALSE)
		cumH <- cumsum(H$count)
		
		knee <- min(which(cumH > thresh * max(cumH)))
		
		selCells <- names(which(assay(sceQ, "exprs")[probe, ] > 
					H$breaks[knee]))
	} else {
		probeMax <- max(assay(sceQ, "exprs")[probe, ])
		selCells <- names(which(assay(sceQ, "exprs")[probe, ] > 
								thresh * probeMax))
	}
	return(selCells)
}

refinePreds.EosBaso <- function(sceQ) {

	# Pre-load with pred0 SOM predictions
	lblTbl$lbl.eosinophil.CCR3.Siglec8 <<- lblTbl$pred.eosinophil.CCR3.Siglec8
	lblTbl$lbl.basophils <<- lblTbl$pred0.basophils

	CCR3p <- getHi(probe = CCR3, thresh = 0.75, fracTot = FALSE)
	CD123p <- getHi(CD123, thresh = 0.75, fracTot = FALSE)
	FceRIp <- getHi(FceRI, thresh = 0.75, fracTot = FALSE)
	Siglec8p <- getHi(Siglec8, thresh = 0.75, fracTot = FALSE)

	# Add cells with high values of selected probes.
	predBaso <- intersect(intersect(CCR3p, CD123p), FceRIp)
	predEos <- setdiff(intersect(CCR3p, Siglec8p), predBaso)

	# Add high-probe cells. 
	# Eosin counts did not improve. Left as is.
	# lblTbl[predEos, "lbl.eosinophil.CCR3.Siglec8"] <<- TRUE
	lblTbl[predBaso, "lbl.basophils"] <<- TRUE

	return()
}


####						Post-process pDC predictions					 ####

refinePreds.DCs <- function(sceQ) {

	DCs <- rownames(lblTbl)[lblTbl$lbl.DC]
	monos <- rownames(lblTbl)[lblTbl$lbl.monocytes]
	
	# tm.14 marks the low-CD14 boundary of monos
	tm.14 <- quantile(assay(sceQ, "exprs")[CD14, monos], probs=0.01)
	
	# Candidate DCs with CD14 > tm.14 are probably monos, not DCs:
	notDCs <- names(which(assay(sceQ, "exprs")[CD14, DCs] > tm.14))
	selDCs <- setdiff(DCs, notDCs)
	
	lblTbl[selDCs, "lbl.DC"] <<- TRUE
	
	
	pCells <- rownames(lblTbl)[lblTbl$lbl.pDC] 
	cCells <- rownames(lblTbl)[lblTbl$lbl.cDC] 
	
	tc.11 <- quantile(assay(sceQ, "exprs")[CD11c, cCells], probs=0.01)
	tp.11 <- quantile(assay(sceQ, "exprs")[CD11c, pCells], probs=0.99)
	tc.123 <- quantile(assay(sceQ, "exprs")[CD123, cCells], probs=0.99)
	tp.123 <- quantile(assay(sceQ, "exprs")[CD123, pCells], probs=0.01)
	
	cDCs <- names(assay(sceQ, "exprs")[CD11c, selDCs] > tc.11 &
				  assay(sceQ, "exprs")[CD123, selDCs] < tc.123)
	pDCs <- names(assay(sceQ, "exprs")[CD11c, selDCs] < tp.11 &
				  assay(sceQ, "exprs")[CD123, selDCs] > tp.123)
	
	lblTbl[cDCs, "lbl.cDC"] <<- TRUE
	lblTbl[pDCs, "lbl.pDC"] <<- TRUE
	
	return()
} # Note lblTbl is revised globally with "<<-"


####						Post-process plasmablast predictions			 ####

refinePreds.Plasmablasts <- function(sceQ) {

	CD19p <- rownames(lblTbl)[lblTbl$lbl.CD19p]
	
	plasmablasts <- rownames(lblTbl)[lblTbl$lbl.plasmablasts]
	
	if (length(plasmablasts) > 0) {
		tp.20 <- max(assay(sceQ, "exprs")[CD20, plasmablasts])
		tp.38 <- min(assay(sceQ, "exprs")[CD38, plasmablasts])
		
		plas <- names(assay(sceQ, "exprs")[CD20, plasmablasts] < tp.20 &
					  assay(sceQ, "exprs")[CD38, plasmablasts] > tp.38)
		plas <- intersect(plas, CD19p)
		
		lblTbl[plasmablasts, "lbl.plasmablasts"] <<- TRUE
	}
	return()
} # Note lblTbl is revised globally with "<<-"


#-------------------------------------------------------------------------------#
#				Select prediction method per cell type, merge.					#
#-------------------------------------------------------------------------------#

labelQueryCells <- function() {

	lblTbl$lbl.CD4 <<- lblTbl$pred.CD4
	lblTbl$lbl.CD8 <<- lblTbl$pred.CD8
	lblTbl$lbl.NK <<- lblTbl$pred.NK | lblTbl$pred0.NK 
	lblTbl$lbl.DC <<- lblTbl$pred.DC 
	lblTbl$lbl.CD11bp.CD16p.neutrophils <<- 
		lblTbl$pred.CD11bp.CD16p.neutrophils

	lblTbl$lbl.cDC <<- lblTbl$pred1.cDC 			   
	lblTbl$lbl.pDC <<- lblTbl$pred1.pDC 			   

	lblTbl$lbl.monocytes <<- lblTbl$pred0.monocytes
	lblTbl$lbl.CD19p <<- lblTbl$pred0.CD19p

	# lblTbl$lbl.CD8.TEMRA <- lblTbl$pred2.CD8.TEMRA

	lblTbl$lbl.CXCR5p <<- lblTbl$pred1.CXCR5p
	lblTbl$lbl.Treg <<- lblTbl$pred1.Treg

	lblTbl$lbl.plasmablasts <<- lblTbl$pred1.plasmablasts & 
								lblTbl$pred2.plasmablasts 

	return()
}


#-------------------------------------------------------------------------------#
#			 Make correlation plots for SELECTED predictions 					#
#-------------------------------------------------------------------------------#

pltPreds <- function(sceL, allLabels, pops) {

	plts = list()
	for (iPop in 1:length(pops)) {
		
		print(paste0(iPop, "   ", Sys.time()))
		
		flo = pred = NULL
		for (qSample in names(sceL)) {
			# Option: use CD45+ cells labeled by both FlowJo and auto-clean:
			shared <- intersect(colnames(sceL[[qSample]]), 
								rownames(allLabels[[qSample]]))
			
			floLabs <- colData(sceL[[qSample]])[ , pops[iPop]] # shared
			flo <- c( flo, sum(floLabs) ) # / length(shared) 
			
			predLabs <- allLabels[[qSample]][ , # shared, 
						paste0("lbl.", pops[iPop])]
			pred <- c( pred, sum(predLabs) ) # / length(shared)
		}
		
		corrRes <- rcorr(flo, pred)
		R <- sprintf("%0.2f", signif(corrRes$r[1, 2], 3))
		P <- sprintf("%0.1E", signif(corrRes$P[1, 2], 3))
		
		titleString <- gsub("CD11bp.CD16p.|.CCR3.Siglec8|s$", "", pops[iPop])
		titleString <- gsub("memory.B.cell", "Bmem", titleString)
		titleString <- gsub(".sw", "SW", titleString)
		titleString <- gsub("TEMRA", "EMRA", titleString)
		titleString <- paste0(titleString, 
							  ", r = ", R,
							  ", p = ", P)
		
		plotDF <- data.frame(FlowJo = flo, predicted = pred)
		
		# If defined, FlowJo outlier counts will be excluded from the fit
		P <- ggplot(plotDF, aes(x = FlowJo, y = predicted)) +
			geom_point(alpha = 0.8, size = 2, color = "coral") +
			geom_smooth(data = subset(plotDF, FlowJo < outliers[iPop]), 
				method = lm, fill = "lightskyblue") +
			ggtitle(titleString) +
			theme(plot.title = element_text(size = 2)) +
			theme_bw() +
			theme(plot.title = element_text(size = 10))
		plts[[iPop]] <- P
	}
	return(plts)
}


#-------------------------------------------------------------------------------#
#			 Make precision-recall plots for SELECTED predictions 				#
#-------------------------------------------------------------------------------#

pltStats <- function(sceL, allLabels, pops2Plt) {

	statsList = list()
	for (qSample in names(sceL)) {
		matchStats = NULL
		for (i in pops2Plt) {
			# Use only CD45+ cells labeled by both FlowJo and auto-clean:
			shared <- intersect(colnames(sceL[[qSample]]), 
								rownames(allLabels[[qSample]]))
			
			resTbl <- table(flo = colData(sceL[[qSample]])[shared, i], 
							pred = allLabels[[qSample]][shared, 
											 paste0("lbl.", i)])
			resTbl <- as.data.frame(resTbl)
			
			TP <- resTbl$Freq[resTbl$flo == TRUE &  resTbl$pred == TRUE]
			FP <- resTbl$Freq[resTbl$flo == FALSE &  resTbl$pred == TRUE]
			FN <- resTbl$Freq[resTbl$flo == TRUE &  resTbl$pred == FALSE]
			
			precision <- signif(TP / (TP + FP), 3) 
			recall <- signif(TP / (TP + FN), 3)
			
			N <- sum(resTbl$Freq[resTbl$flo == TRUE])
			
			matchStats <- rbind(matchStats, 
						  c(precision = precision, 
						  recall = recall, 
						  N = N))
		}
		rownames(matchStats) <- pops2Plt
		statsList[[qSample]] <- as.data.frame(matchStats)
		
		statsList[[qSample]]$Sample <- qSample
		statsList[[qSample]]$population <- pops2Plt
	}

	plotTbl <- do.call(rbind, statsList)
	rownames(plotTbl) <- NULL
	plotTbl$Sample <- factor(plotTbl$Sample, levels = unique(plotTbl$Sample))
	plotTbl$population <- gsub(".CCR3.Siglec8|CD11bp.CD16p.", 
							   "", plotTbl$population)

	saveRDS(plotTbl, file = "matchStats.RDS")
	
	plt = ggplot(plotTbl, aes(x = precision, y = recall, color = population)) +
		geom_point(size = 3) + 
		scale_color_manual(values=alpha(colors30, 0.75)) +
		scale_y_continuous(breaks = seq(0,10)/10, limits = c(0.5, 1)) +
		scale_x_continuous(breaks = seq(0,10)/10, limits = c(0.5, 1)) +
		xlab("precision ( == 1 - FDR )") +
		theme_bw()

	return(plt)
}


#-------------------------------------------------------------------------------#
#			 Make precision-recall plots for LOO predictions 					#
#-------------------------------------------------------------------------------#

pltStatsLOO <- function(allSCE, allLabels, pops2Plt) {

	statsList = list()
	for (qSample in IDs) {

		sceQ <- allSCE[ , 
				grep(paste0(qSample, "\\."), colnames(allSCE))]
		
		matchStats = NULL
		for (i in pops2Plt) {
			resTbl <- table(flo = colData(sceQ)[ , i], 
							pred = allLabels[[qSample]][ , 
											 paste0("lbl.", i)])
			resTbl <- as.data.frame(resTbl)
			
			TP <- resTbl$Freq[resTbl$flo == TRUE &  resTbl$pred == TRUE]
			FP <- resTbl$Freq[resTbl$flo == FALSE &  resTbl$pred == TRUE]
			FN <- resTbl$Freq[resTbl$flo == TRUE &  resTbl$pred == FALSE]
			
			precision <- signif(TP / (TP + FP), 3) 
			recall <- signif(TP / (TP + FN), 3)
			
			N <- sum(resTbl$Freq[resTbl$flo == TRUE])
			
			matchStats <- rbind(matchStats, 
						  c(precision = precision, 
						  recall = recall, 
						  N = N))
		}
		rownames(matchStats) <- pops2Plt
		statsList[[qSample]] <- as.data.frame(matchStats)
		
		statsList[[qSample]]$Sample <- qSample
		statsList[[qSample]]$population <- pops2Plt
	}

	plotTbl <- do.call(rbind, statsList)
	rownames(plotTbl) <- NULL
	plotTbl$Sample <- factor(plotTbl$Sample, levels = unique(plotTbl$Sample))
	plotTbl$population <- gsub(".CCR3.Siglec8|CD11bp.CD16p.", 
							   "", plotTbl$population)

	saveRDS(plotTbl, file = "matchStats.RDS")
	
	plt = ggplot(plotTbl, aes(x = precision, y = recall, color = population)) +
		geom_point(size = 3) + 
		scale_color_manual(values=alpha(colors30, 0.75)) +
		scale_y_continuous(breaks = seq(0,10)/10, limits = c(0.5, 1)) +
		scale_x_continuous(breaks = seq(0,10)/10, limits = c(0.5, 1)) +
		xlab("precision ( == 1 - FDR )") +
		theme_bw()

	return(plt)
}















