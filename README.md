# superSOM
## CyTOF sample gating using supervised Self Organizing Maps

The scripts in this repository use various R and Bioconductor packages (most notably `flowSOM' and `CATALYST`) to enable automated gating of raw mass cytometery (CyTOF) FCS files. They were developed specifcally for whole-blood CyTOF immune profiling and have not been tested on other sample types.The overalll aproach is summarized in the figure below:

Briefly, a query (ungated) sample is SOM-clustered with N manually gated and labeled samples. In each SOM cluster with > P% cells of the same label, unlabeled cells are given the same label as the majority label. The number of training samples (N) and the cluster purity thrshold will vary depending on the data characteristics. In our data, we find 15 < N < 20 and 0.5 < P < 0.8 work well.

### Automated 'clean up' gating

The superSOM pipeline has 3 parts. `CleanUP' gating takes in raw FCA files, performs Gaussian Mixture Modeling to extract CD45+ singlets, and generates a Single Cell Experiment list for auto-labeling. See example images below:

![cleanUPGating](https://user-images.githubusercontent.com/46689973/104108847-17e1dd80-527d-11eb-990b-7ba650d1bdf3.png)

Since the results of clean up will affect downstream prediction accuracy, it is worth optimizing this step fpr each use case. In genertal, we find clean up identifies CD45+ singlets with > 90% Precisio  and Recall, as in the example figure below. 

![cleanUp_accuracy](https://user-images.githubusercontent.com/46689973/104108903-714a0c80-527d-11eb-8f30-b3f3ec241565.png)

### Supervised SOM clustering

To label the 9 major (mutually-exclusive) immune populations, we simply cluster all cells using all probes. WE find that this approach generally produces highly `pure' clusters of cells, as illustrated in the figures below.

![majPopClusters](https://user-images.githubusercontent.com/46689973/104108917-87f06380-527d-11eb-832f-edf84fdf9406.png)

![purity_profiles](https://user-images.githubusercontent.com/46689973/104108919-8a52bd80-527d-11eb-8605-b72548a5417f.png)

![ROIs](https://user-images.githubusercontent.com/46689973/104108926-92126200-527d-11eb-92ab-5be7293bcf31.png)

![postProcessing](https://user-images.githubusercontent.com/46689973/104108928-976fac80-527d-11eb-8fbc-3c0a081beb8a.png)
