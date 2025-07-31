metanr_packages <- function(){
  
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea","devtools","crmn","httr","qs")
  
  list_installed <- installed.packages()
  
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  
  if(length(new_pkgs)!=0){
    
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}

metanr_packages()

library(devtools)

# Step 2: Install MetaboAnalystR with documentation
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = TRUE, build_manual =T)

# Load MetaboAnalystR
library(MetaboAnalystR)

# Clean global environment
rm(list = ls())
mSet<-InitDataObjects("conc", "stat", FALSE);

mSet<-Read.TextData(mSet, "final_positive_negative_comb.csv", "rowu", "disc");
df = read.csv("final_positive_negative_comb.csv")

mSet<-SanityCheckData(mSet);

print(unique(mSet$dataSet$cls))

print(table(mSet$dataSet$cls))

mSet<-ReplaceMin(mSet);
mSet<-PreparePrenormData(mSet);
mSet<-Normalization(mSet, "NULL", "LogNorm", "MeanCenter", "S10T0", ratio=FALSE, ratioNum=20);

mSet<-PlotNormSummary(mSet, "norm_0_", format ="png", dpi=72, width=NA);
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "png", dpi=72, width=NA);


levels(mSet$dataSet$cls)
levels(mSet$dataSet$cls) <- c("Infected", "Control")
mSet <- FC.Anal.unpaired(mSet, 2.0)   # 2.0 is the fold-change threshold



# Plot fold-change analysis
mSet<-PlotFC(mSet, "fc_0_", "png", 72, width=NA)

# To view fold-change 
mSet$analSet$fc$fc.log

# Run unpaired parametric t-test with threshold 0.05 (p-value cutoff)
mSet <- Ttests.Anal(mSet, nonpar = FALSE, threshp = 0.05, paired = FALSE, equal.var = TRUE)

# Manually adjust p-values (FDR - Benjamini-Hochberg)
mSet$analSet$tt$adj.p.value <- p.adjust(mSet$analSet$tt$p.value, method = "fdr")

# Plot of the T-test results
mSet<-PlotTT(mSet, imgName = "tt_0_", format = "png", dpi = 72, width=NA)


# Perform the volcano analysis
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 0, F, 0.1, TRUE, "raw")

# Create the volcano plot
mSet<-PlotVolcano(mSet, "volcano_0_",1, format ="png", dpi=72, width=NA)

# Perform ANOVA
mSet <- ANOVA.Anal(mSet, F, 0.05, "fisher")

# Plot ANOVA
mSet <- PlotANOVA(mSet, "aov_0_", "png", 72, width=NA)

### OPTION 1 - Heatmap specifying pearson distance and an overview
mSet<-PlotCorrHeatMap(mSet, "corr_0_", "png", 72, width=NA, "col", "pearson", "bwm", "overview", F, F, 0.0)
### OPTION 2 - Heatmap specifying pearson correlation and a detailed view
mSet<-PlotCorrHeatMap(mSet, "corr_1_", format = "png", dpi=72, width=NA, "col", "spearman", "bwm", "detail", F, F, 100, 100)

# Perform PCA analysis
mSet<-PCA.Anal(mSet)

# Create PCA overview
mSet<-PlotPCAPairSummary(mSet, "pca_pair_0_", format = "png", dpi = 72, width=NA, 5)

# Create PCA scree plot
mSet<-PlotPCAScree(mSet, "pca_scree_0_", "png", dpi = 72, width=NA, 5)

# Create a 2D PCA score plot
mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_", format = "png", dpi=72, width=NA, 1, 2, 0.95, 1, 0)
mSet <- PlotPCA2DScore(
  mSet,
  imgName = "pca_score2d_0_",
  format = "png",
  dpi = 72,
  width = NA,
  style = 1,
  pcx = 1,
  pcy = 2,
  reg = 0.95,
  show = 1,
  grey.scale = 0
)

library(crosstalk)
# Create a 3D PCA score plot
mSet<-PlotPCA3DScoreImg(mSet, "pca_score3d_0_", "png", 72, width=NA, 1,2,3, 40)

# Create a PCA loadings Plots
mSet<-PlotPCALoading(mSet, "pca_loading_0_", "png", 72, width=NA, 1,2);

# Create a PCA Biplot
mSet<-PlotPCABiplot(mSet, "pca_biplot_0_", format = "png", dpi = 72, width=NA, 1, 2)

mSet$imgSet$pca.3d
# Partial Least Squares - Discriminant Analysis (PLS-DA)

mSet<-PLSR.Anal(mSet, reg=TRUE)

mSet<-PlotPLSPairSummary(mSet, "pls_pair_0_", "png", 72, width=NA, 5)

mSet<-PlotPLS2DScore(mSet, "pls_score2d_0_", "png", 72, width=NA, 1,2,0.95,1,0)

mSet<-PlotPLS3DScoreImg(mSet, "pls_score3d_0_", "png", 72, width=NA, 1,2,3, 40)

mSet<-PlotPLSLoading(mSet, "pls_loading_0_", "png", 72, width=NA, 1, 2);

mSet <- PLSDA.CV(mSet, methodName = "T", compNum = 5, choice = "Q2")

mSet<-PlotPLS.Classification(mSet, "pls_cv_0_", "png", 72, width=NA)

mSet<-PlotPLS.Imp(mSet, "pls_imp_0_", "png", 72, width=NA, "vip", "Comp. 1", 15, FALSE)

mSet<-PLSDA.Permut(mSet, 100, "accu")

mSet<-PlotPLS.Permutation(mSet, "pls_perm_1_", "png", 72, width=NA)

# Perform sPLS-DA analysis
mSet<-SPLSR.Anal(mSet, 5, 10, "same", "Mfold", 5, T)

#Sparse Partial Least Squares - Discriminant Analysis (sPLS-DA)

# Plot sPLS-DA overview
#mSet<-PlotSPLSPairSummary(mSet, "spls_pair_0_", format = "png", dpi=72, width=NA, 5)
mSet <- SPLSR.Anal(mSet, 5, 10, "same", "Mfold")

# Create 2D sPLS-DA Score Plot
mSet<-PlotSPLS2DScore(mSet, "spls_score2d_0_", format = "png", dpi=72, width=NA, 1, 2, 0.95, 1, 0)

# Create 3D sPLS-DA Score Plot
mSet<-PlotSPLS3DScoreImg(mSet, "spls_score3d_0_", format = "png", 72, width=NA, 1, 2, 3, 40)

# Create sPLS-DA loadings plot
mSet<-PlotSPLSLoading(mSet, "spls_loading_0_", format = "png", dpi=72, width=NA, 1,"overview")

# Perform cross-validation and plot sPLS-DA classification
mSet<-PlotSPLSDA.Classification(mSet, "spls_cv_0_", format = "png", dpi=72, width=NA)

#Orthogonal Partial Least Squares - Discriminant Analysis (orthoPLS-DA)

# Perform oPLS-DA analysis
mSet<-OPLSR.Anal(mSet, reg=TRUE)
mSet <- OPLSR.Anal(mSet)
str(mSet$analSet$opls.reg)

# Create a 2D oPLS-DA score plot
mSet<-PlotOPLS2DScore(mSet, "opls_score2d_0_", format = "png", dpi=72, width=NA, 1,2,0.95,1,0)

# Create a significant features plot
#mSet<-PlotOPLS.Splot(mSet, "opls_splot_0_", "all", "png", 72, width=NA);
mSet <- PlotOPLS.Splot(mSet, "opls_splot_0_", "all", "png", 72, width = NA)


# Create a plot of features ranked by importance
#mSet<-PlotOPLS.Imp(mSet, "opls_imp_0_", "png", 72, width=NA, "vip", "tscore", 15,FALSE)
#mSet <- PlotOPLS.Imp(mSet, "opls_imp_0_", "png", 72, width = NA, "vip", "tscore", 15, FALSE)

# Create a plot of the model overview
mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", format = "png", dpi=72, width=NA)

# Perform and plot oPLS-DA permutation 
mSet<-OPLSDA.Permut(mSet, num =100)

mSet<-PlotOPLS.Permutation(mSet, "opls_perm_2_", format = "png", dpi=72, width=NA)


#Significance Analysis of Microarrary (and Metabolites) (SAM)

# Perform SAM analysis
mSet<-SAM.Anal(mSet, "d.stat", FALSE, TRUE, 0.0, "sam_imp_0_")



# Create a SAM plot of FDR values
mSet<-PlotSAM.FDR(mSet, "sam_view_0_", format = "png", dpi=72, width=NA)

# Create a SAM plot of results
mSet<-PlotSAM.Cmpd(mSet, "sam_imp_0_", format = "png", dpi=72, width=NA)

# Empirical Bayesian Analysis of Microarray (and Metabolites) (EBAM)


#Empirical Bayesian Analysis of Microarray (and Metabolites) (EBAM)
# Perform EBAM analysis, plot EBAM analysis and create the EBAM matrix of significant features
mSet<-EBAM.Init(mSet, FALSE, TRUE, FALSE, -99.0, 0.9, "ebam_view_0_", "ebam_imp_0_")
# Create a EBAM plot of results
PlotEBAM.Cmpd(mSet, "ebam_imp_0_", "png", 72, width=NA)

#Hierarchical Clustering: Dendogram

# Perform hierarchical clustering and plot dendogram
#mSet<-PlotHCTree(mSet, "tree_0_", format = "png", dpi=72, width=NA, "euclidean", "ward.D")

# Perform hierarchical clustering and plot heat map
mSet<-PlotHeatMap(mSet, "heatmap_0_", "png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8, "overview", T, T, NULL, T, F, T, T, T)
# Create a vector with the indices of the 50 features you want to show (e.g., the first 50)
mSet <- PlotHeatMap(
  mSet,
  "heatmap_0_",
  "png",
  72,
  width = NA,
  "norm",
  "row",
  "euclidean",
  "ward.D",
  "bwm",
  "overview",
  TRUE,
  TRUE,
  NA,
  TRUE,
  FALSE
)

feature_indices <- 1:50

# Call PlotHeatMap with var.inx set to feature_indices
mSet <- PlotHeatMap(
  mSet,
  "heatmap_0_",
  "png",
  72,
  width = NA,
  "norm",
  "row",
  "euclidean",
  "ward.D",
  "bwm",
  "overview",
  TRUE,      # rowV
  TRUE,      # colV
  feature_indices,  # var.inx: display these features
  TRUE,      # border
  FALSE      # grp.ave
)


# Perform K-means analysis
mSet<-Kmeans.Anal(mSet, 3)

# Plot K-means analysis 
mSet<-PlotKmeans(mSet, "km_0_", format = "png", dpi=72, width=NA)

# Perform SOM analysis
mSet<-SOM.Anal(mSet, 1, 3,"linear","gaussian")

# Plot SOM analysis
mSet<-PlotSOM(mSet, "som_0_", format = "png", dpi=72, width=NA)

# Perform random forest analysis
mSet<-RF.Anal(mSet, 500, 7, 1)

# Plot random forest classification
mSet<-PlotRF.Classify(mSet, "rf_cls_0_", format = "png", dpi=72, width=NA)

# Plot random forest variables of importance
mSet<-PlotRF.VIP(mSet, "rf_imp_0_", format = "png", dpi=72, width=NA)

# Plot random forest outliers 
mSet<-PlotRF.Outlier(mSet, "rf_outlier_0_", format = "png", dpi=72, width=NA)


# Perform SVM 
mSet<-RSVM.Anal(mSet, 10)

mSet<-PlotRSVM.Classification(mSet, "svm_cls_0_", format = "png", dpi=72, width=NA)

mSet<-PlotRSVM.Cmpd(mSet, "svm_imp_0_", format = "png", dpi=72, width=NA)
