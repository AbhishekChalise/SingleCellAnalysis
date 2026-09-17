Single-Cell RNA-Seq Analysis Pipeline: COVID-19 Lung Atlas

Overview

This repository contains a complete single-cell RNA sequencing (scRNA-seq) analysis pipeline written in Python. The dataset used is a COVID-19 lung tissue sample (Melms et al., 2021). The pipeline progresses from raw RNA count matrices to clustered and annotated cell types using Deep Learning for doublet removal and graph-based clustering.

Tech Stack: Python, Scanpy, scvi-tools (PyTorch), Pandas, NumPy, Seaborn.

Pipeline Steps

Phase 1: Data Loading & Basic Cleaning

Goal: Remove empty droplets and unexpressed genes to reduce memory usage.

     Data Loading: Loaded raw integer mRNA counts from a CSV file and transposed it so cells are rows (n_obs) and genes are columns (n_vars).
     Basic Filtering: 
         sc.pp.filter_cells(adata, min_genes=200): Removed cells with fewer than 200 active genes (empty droplets).
         sc.pp.filter_genes(adata, min_cells=3): Removed genes appearing in fewer than 3 cells (technical noise).

Phase 2: Quality Control (QC)

Goal: Remove dead cells and extreme outliers.

     Mitochondrial QC: Calculated the percentage of Mitochondrial RNA per cell (pct_counts_mt). Cells with > 20% MT RNA are dead/dying and were removed.
     Ribosomal QC: Downloaded a curated list of 88 Ribosomal genes from the Broad Institute (MSigDB). Calculated pct_counts_ribo to check for technical noise.
     Outlier Removal: Used the 98th percentile of n_genes_by_counts to mathematically remove extreme high-count outliers (potential doublets) instead of guessing a hard cutoff number.

Phase 3: Doublet Removal (Deep Learning)

Goal: Identify and remove "fake" cells (two cells stuck in one droplet).

     Feature Extraction (scVI): Subsetted the data to 2,000 Highly Variable Genes using flavor='seurat_v3' (which requires raw counts). Trained a Variational Autoencoder (scvi.model.SCVI) to compress the 2,000 genes into a dense mathematical latent space.
     Classification (SOLO): Passed the trained VAE to scvi.external.SOLO. SOLO simulated fake doublets and trained a Neural Network classifier to distinguish between singlets and doublets.
     Manual Thresholding: Extracted raw logits (return_logits=True) from the predictions. Calculated the difference (dif = doublet - singlet). Plotted a histogram and kept only cells where dif > 1 to ensure high-confidence doublet removal.

Phase 4: Normalization & Scaling

Goal: Make the data mathematically fair for Machine Learning algorithms.

     Normalization: sc.pp.normalize_total(target_sum=1e4). Scaled every cell so they all have exactly 10,000 total mRNA molecules, removing the bias of cell size/capture efficiency.
     Log Transform: sc.pp.log1p(). Applied a natural logarithm to flatten extreme numbers, converting sparse count data into continuous decimals.
     Feature Selection: sc.pp.highly_variable_genes(n_top_genes=2000, subset=True, flavor='seurat'). Selected the top 2,000 genes with the highest variance. Discarded the remaining ~17,000 "housekeeping" genes.
     Scaling: sc.pp.scale(max_value=10). Converted all gene expression to Z-scores (mean=0, variance=1). Clipped extreme outliers to ±10 to prevent them from breaking the PCA math.

Phase 5: Dimensionality Reduction (PCA)

Goal: Compress the 2,000 genes into a smaller, denser mathematical space.

     PCA: sc.tl.pca(n_comps=50). Used linear algebra to compress the 2,000 dimensions into 50 Principal Components (PCs).
     Elbow Plot: sc.pl.pca_variance_ratio(n_pcs=50). Visually inspected the plot to find where the biological signal turns into noise. Selected the first 29 PCs for downstream clustering.

Phase 6: Graph Construction & Clustering

Goal: Group the cells by biological type.

     KNN Graph: sc.pp.neighbors(n_neighbors=20, n_pcs=29). Calculated the Euclidean distance between all cells in the 29-dimensional space. Drew 20 lines (edges) connecting each cell to its 20 closest neighbors, creating a giant interconnected web.
     Leiden Clustering: sc.tl.leiden(resolution=1.0). Ran a community detection algorithm on the KNN web. The algorithm found tightly connected "knots" and assigned a Cluster ID (Cluster 0, Cluster 1, etc.) to every cell.

Phase 7: Visualization

Goal: Allow humans to see the high-dimensional data.

     UMAP Calculation: sc.tl.umap(). Compressed the 29-dimensional space down to exactly 2 dimensions (X and Y coordinates) for every cell.
     UMAP Plotting: sc.pl.umap(color='leiden'). Plotted the 2D coordinates as a scatter plot. Colored the dots based on their Leiden Cluster ID. Dots close together have similar RNA.

Phase 8: Biological Discovery (Differential Expression)

Goal: Translate math into biology.

     Rank Genes: sc.tl.rank_genes_groups(groupby='leiden', method='wilcoxon'). Ran a Wilcoxon statistical test to find genes uniquely high in each cluster compared to the rest.
     Cell Type Annotation: Created a Pandas DataFrame of the top marker genes. By searching biological databases (e.g., SFTPB = Lung cells, CD247 = T-cells, MRC1 = Macrophages), manually mapped cluster IDs to real biological cell types.
     
