## Data repository of paper titled: Cardinality constrained weight based PCA

### Abstract of the paper

Principal component analysis (PCA) (Jolliffe, 1986) is a widely used analysis technique for dimension reduction. The components resulting from a PCA are linear combinations of all variables. This can make their interpretation difficult, especially when the number of variables in a data set is large. To that end, so-called regularized sparse PCA methods have been developed that aim at reducing the number of non-zero coefficients by relying on shrinkage penalties implying that only a subset of the variables makes up the components. The problem that shrinkage methods solve is not that of finding the best subset of variables optimizing the PCA criterion. In this paper we present cardinality constraint PCA (CCPCA) to solve the best subset PCA problem by fixing the number of non-zero coefficients through a cardinality constraint approach. For this purpose, we propose using the cardinality constraints regression algorithm from Adachi and Kiers (2017) and Bertsimas (2016). Consistent with results obtained for regression analysis, we found that CCPCA outperforms sparse PCA based on shrinkage penalties (Zou, 2006) when noise levels are low but not when noise levels are high.

Authors: Niek C. de Schipper, A Tonne, Katrijn Van Deun

#### Description of files:

- simstudy.R: Script to perform the simulation study
- analyze_results.R: Script to create the boxplots from the simulated data 
- analyze_results_bias.R: Script to create the tables in the paper 
- tucker.R: Script to create the tables in the paper 
- pca_ccreg.cpp: the CCPCA function as used in the paper, users should use the sparseWeightBasedPCA package 
- pca_ccreg.cpp: the CCPCA function as used in the paper, users should use the `sparseWeightBasedPCA` package 
- pcaccregPackage: contains the pca_ccreg.cpp function wrapped in a package, users should use the `sparseWeightBasedPCA` package 

#### Not included (send an email if you want them)
- sim_res: Raw results from the simulation study
- sim_res_rawdata: Raw data from the simulation study 

