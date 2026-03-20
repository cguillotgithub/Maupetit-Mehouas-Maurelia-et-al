# Maupetit-Mehouas-Maurelia-et-al- 2026

# Figure 5 — No integration (Seurat + scDblFinder)

This repository contains the R Markdown workflow used to generate **Figure 5 (no integration)** from three scRNA-seq count matrices (CSV), using a standard Seurat pipeline (LogNormalize → HVGs → scaling → PCA → neighbors/clusters → UMAP) followed by **doublet detection with scDblFinder** and downstream plotting/exports.

## What this does

From three input CSV matrices (genes x cells):

1. Creates a Seurat object per sample (`SBC1`, `SBC3`, `SBC7`) and merges them.
2. Computes mitochondrial fraction (`percent.mt`) using features matching `^MT.`.
3. Filters cells with:
   - `percent.mt < 25`
   - `1000 < nFeature_RNA < 10000`
4. Joins layers (Seurat v5) after subsetting.
5. Removes features whose names contain `SBC` (barcode-like artifacts) from the RNA assay.
6. Runs the “no integration” Seurat workflow:
   - `NormalizeData` (LogNormalize; default)
   - `FindVariableFeatures` (vst; 1500 features)
   - `ScaleData` on **VariableFeatures** (memory-friendly)
   - `RunPCA` (50 PCs)
   - `FindNeighbors` (dims 1:50)
   - `FindClusters` (resolution = 0.3; `random.seed = 1`)
   - `RunUMAP` (dims 1:20; `n.neighbors = 50`; `seed.use = 1`)
7. Runs **scDblFinder** on the Seurat object converted to SingleCellExperiment, stores classification + score in metadata, and keeps **singlets only**.
8. Generates UMAPs (pre- and post-doublet filtering), a layered sample-colored UMAP, marker FeaturePlots, cluster composition test plot, and saves RDS + markers table.
