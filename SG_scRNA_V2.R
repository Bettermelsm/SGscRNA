#############################################
# Integrated scRNA-seq Pipeline (SG + FAN)
# Batch-Processing + Plot Auto-Export + SingleR
# EDIT: Sengene / 20251210
#############################################

############################################################
# E0. ENVIRONMENT
############################################################
options(stringsAsFactors = FALSE)
rm(list = ls()); gc()

############################################################
# E1. BASIC SETTINGS
############################################################
base_dir <- "/mnt/sda1/LSM_scRNAseq/20250527IM"
setwd(base_dir)
dir.create("Ana_res", showWarnings = FALSE)
outdir <- file.path(base_dir, "Ana_res")

# 自动保存绘图函数
save_plot <- function(fname, plt, w=2400, h=1600, res=300){
  png(filename = file.path(outdir, fname), width=w, height=h, res=res)
  print(plt)
  dev.off()
}

############################################################
# E2. LOAD PACKAGES
############################################################
library(Seurat)
library(hdf5r)
library(harmony)
library(Matrix)
library(irlba)
library(clustree)
library(clusterProfiler)
library(org.Hs.eg.db)
library(monocle3)
library(CellChat)
library(SingleR)
library(celldex)

library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)

############################################################
# PART1. DATA IMPORT (BATCH) & MERGING
############################################################
sample_folders <- c("sham1","sham2","sham3",
                    "model1","model2","model3",
                    "SLTNF1","SLTNF2","SLTNF3")

data_dir <- "/mnt/sda1/LSM_scRNAseq/202512MCAO/MCAOmtx/workflow_results/01_BasicAnalysis/Dnbc4tools"
seurat_list <- list()

for (smp in sample_folders){
  message("Loading sample: ", smp)
  
  spath <- file.path(data_dir, smp, "filter_matrix")
  mat <- Read10X(data.dir = spath, gene.column = 1)
  
  obj <- CreateSeuratObject(counts = mat,
                            project = smp,
                            min.cells = 3,
                            min.features = 200)
  obj$sample <- smp
  obj$group  <- case_when(
    grepl("^sham", smp)  ~ "sham",
    grepl("^model", smp) ~ "model",
    grepl("^SLTNF", smp) ~ "SLTNF",
    TRUE ~ "other"
  )
  
  seurat_list[[smp]] <- obj
}

# Normalize before integration (FAN style)
seurat_list <- lapply(seurat_list, function(x){
  NormalizeData(x) %>% FindVariableFeatures(selection.method="vst", nfeatures=2000)
})

anchors <- FindIntegrationAnchors(seurat_list, dims=1:30)
obj_int  <- IntegrateData(anchorset = anchors, dims=1:30)

DefaultAssay(obj_int) <- "integrated"

############################################################
# PART2. QC
############################################################
DefaultAssay(obj_int) <- "RNA"

# mt %
obj_int[["percent.mt"]] <- PercentageFeatureSet(obj_int, pattern="^mt-")

# RBC 基因
HB.genes.mouse <- c("Hba-a1","Hba-a2","Hbb-bs","Hbb-bt","Hbb-bh1","Hbb-bh2","Hba-x")
hb_present <- HB.genes.mouse[HB.genes.mouse %in% rownames(obj_int)]
obj_int[["percent.HB"]] <- PercentageFeatureSet(obj_int, features = hb_present)

# QC 图
save_plot("P2.1_QC_violin.png",
          VlnPlot(obj_int, features=c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"), ncol=4))

# 筛选（参数化，可供 Django 前端传入）
obj_qc <- subset(obj_int, subset = 
                   nFeature_RNA > 300 &
                   nFeature_RNA < 6000 &
                   percent.mt < 10)

############################################################
# PART3. NORMALIZATION (FAN STYLE)
############################################################
obj_norm <- NormalizeData(obj_qc)
obj_hvg  <- FindVariableFeatures(obj_norm, selection.method="vst", nfeatures=3000)
obj_scl  <- ScaleData(obj_hvg, features = rownames(obj_hvg))

############################################################
# PART4. PCA
############################################################
obj_pca <- RunPCA(obj_scl, features = VariableFeatures(obj_scl))

save_plot("P4.1_PCA_plot.png", DimPlot(obj_pca, reduction="pca"))
save_plot("P4.2_ElbowPlot.png", ElbowPlot(obj_pca))

############################################################
# PART5. CLUSTERING + UMAP
############################################################
obj_pca <- FindNeighbors(obj_pca, dims=1:30)
obj_pca <- FindClusters(obj_pca, resolution = 0.8)
obj_pca <- RunUMAP(obj_pca, dims=1:30)

save_plot("P5.UMAP_group.png", DimPlot(obj_pca, group.by="group"))
save_plot("P5.UMAP_cluster.png", DimPlot(obj_pca, group.by="seurat_clusters"))

############################################################
# PART6. Marker
############################################################
markers <- FindAllMarkers(obj_pca, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file.path(outdir, "Markers_all_clusters.csv"))

############################################################
# PART7. SingleR 注释（FAN 版保留）
############################################################
ref <- MouseRNAseqData()   # 小鼠参考（如果人类数据可改为 HPCA）
obj_sce <- as.SingleCellExperiment(obj_pca)

pred <- SingleR(test = obj_sce,
                ref  = ref,
                labels = ref$label.main)

obj_pca$SingleR <- pred$labels

save_plot("P7.SingleR_annotation.png",
          DimPlot(obj_pca, group.by="SingleR", label=TRUE))

############################################################
# PART8. SAVE OBJECT
############################################################
saveRDS(obj_pca, file.path(outdir, "Final_Integrated_Object.rds"))

message("==== Pipeline finished ====")
