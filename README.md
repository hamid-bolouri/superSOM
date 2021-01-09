# superSOM
## CyTOF sample gating using supervised Self Organizing Maps

### Overview

The scripts in this repository use various R and Bioconductor packages (most notably `flowSOM`, `CATALYST`, and most of the Bioc flow libraries!) to enable automated gating of raw mass cytometery (CyTOF) FCS files. They were developed specifcally for whole-blood CyTOF immune profiling and have not been tested on other sample types.The overalll aproach is summarized in the figure below:

Briefly, a query (ungated) sample is SOM-clustered with N manually gated and labeled samples. In each SOM cluster with > P% cells of the same label, unlabeled cells are given the same label as the majority label. The number of training samples (N) and the cluster purity thrshold will vary depending on the data characteristics. In our data, we find 15 < N < 20 and 0.5 < P < 0.8 work well.

### Pre-requisites

superSOM requires the following packages. 

For `cleanUP` gating:

```
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
library(destiny)
library(scales)	
library(mclust)
library(MASS)
library(IDPmisc)
library(Hmisc)
```

For supervised SOME clustering:

```
library(SingleCellExperiment)
library(flowCore)
library(flowWorkspace)
library(ggcyto)
library(CytoML)
library(FlowSOM)
library(CATALYST)
library(mclust)
library(MASS)
library(batchelor)
library(scMerge)
library(scater)
library(cowplot)
library(readxl)
library(reshape2)
library(tidyverse)
library(BiocParallel) 
```

For downstream analysis:
```
library(Hmisc)
library(pheatmap)
```

The current implementation of superSOM has not been optimized for speed or memory use. We typically use `MulticoreParam(workers = 16)` on a Unix platform with > 64 GB of RAM.

### Automated 'clean up' gating

The superSOM pipeline has 3 parts. `CleanUP' gating takes in raw FCA files, performs Gaussian Mixture Modeling to extract CD45+ singlets, and generates a Single Cell Experiment list for auto-labeling. See example images below:

![cleanUPGating](https://user-images.githubusercontent.com/46689973/104108847-17e1dd80-527d-11eb-990b-7ba650d1bdf3.png)

Since the results of clean up will affect downstream prediction accuracy, it is worth optimizing this step fpr each use case. In genertal, we find clean up identifies CD45+ singlets with > 90% Precision  and Recall, as in the example figure below. 

![cleanUp_accuracy](https://user-images.githubusercontent.com/46689973/104108903-714a0c80-527d-11eb-8f30-b3f3ec241565.png)

### Supervised SOM clustering

To label the 9 major (mutually-exclusive) immune populations, we simply cluster all cells using all probes. WE find that this approach generally produces highly `pure' clusters of cells, as illustrated in the figures below.

![majPopClusters](https://user-images.githubusercontent.com/46689973/104108917-87f06380-527d-11eb-832f-edf84fdf9406.png)

![purity_profiles](https://user-images.githubusercontent.com/46689973/104108919-8a52bd80-527d-11eb-8605-b72548a5417f.png)

In addition to using all probes for clustering, superSOM also clusters the major population using the specific sets of features that define each population, as specified in the gating hierarchy. Example selected features are given below:

![strategy](https://user-images.githubusercontent.com/46689973/104108921-8c1c8100-527d-11eb-80e4-19665af258e3.png)

To gate/label sub-populations of cells, superSOM simply repeats the above procedure using the newly-labeled parent cell populations for the query sample(s) and feature-sets defined by the gating hierarchy. For cells defined by quadrant gates, we also use the collection of all gate boundaries in the training data to define `Regions of Interest` to define the parent cell populations more tightly, as illustrated in the example below:

![ROIs](https://user-images.githubusercontent.com/46689973/104108926-92126200-527d-11eb-92ab-5be7293bcf31.png)

After superSOM has perfomred all the different types of clustering/labeling described above, the user can select the particular combination of inputs and parameters that lead to the best (application-specific) precision and recall values for each cell population.

### Post-processing

After completion of SOM-clustering and label transfer, superSOM performs a post-processing step in which we adjust the boundaries of the clustering based gates to more closely resemble those defined by rectangular gates in FlowJo, as in the example below:

![postProcessing](https://user-images.githubusercontent.com/46689973/104108928-976fac80-527d-11eb-8fbc-3c0a081beb8a.png)

### How to use superSOM

superSOM describes an approach to automated gating developed for a particular use-case. It is not a general-purpose R/Bioc package. The code in this repository is intended as an example implementation of the underlying idea and can be used as a strating point for future applications.


