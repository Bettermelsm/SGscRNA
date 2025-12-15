#############################################
# SG_scRNA_Seq_FLOW_V2
# Batch-Processing + Plot Auto-Export + SingleR
# EDIT_by_Sengene/Ribosome
# Update_DATE202512111752
#############################################

############################################################
# E0. ENVIRONMENT
############################################################
options(stringsAsFactors = FALSE)
rm(list = ls()); gc()

############################################################
# E1. BASIC SETTINGS
############################################################
base_dir <- "/ur/local/PATH"
setwd(base_dir)
dir.create("Ana_res2", showWarnings = FALSE)
outdir <- file.path(base_dir, "Ana_res2")

# saving_plots_automatically
save_plot <- function(fname, plt, w=2600, h=1600, res=300){
  png(filename = file.path(outdir, fname), width=w, height=h, res=res)
  print(plt)
  dev.off()
}

############################################################
# E2. LOAD PACKAGES
############################################################
#E1.1data_preprocessing_clustering
library(Seurat)
packageVersion("Seurat")##Seurat_5.3.0
#library(hdf5r)
library(harmony)
library(Matrix)

#E2.2PCA&Normalization
library(irlba)
library(clustree)

#E2.3Feature enrichment/annotation
library(clusterProfiler)
library(org.Hs.eg.db)

#E2.4Pseudotime&intercellular_communication
#library(monocle3)
library(ggalluvial)
#library(CellChat)

#E2.5Visualization
library(dplyr)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(stringr)
library(ggpubr)
library(data.table)

#E2.6PCA_ARI
library(mclust) 
library(purrr)

############################################################
# PART1. DATA IMPORT (BATCH) & MERGING
############################################################
sample_folders <- c("sample1","sham2","sham3",
                    "model1","model2","model3",
                    "Treat1","Treat2","Treat3")

data_dir <- base_dir
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

# Normalize before integration
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
obj_int[["percent.mt"]] <- PercentageFeatureSet(
  obj_int, 
  pattern = "^[mM][tT][ -]?"
)

# RBC_genes
get_HB_genes <- function(species = c("human", "rabbit", "mouse", "rat")) {  
  species <- match.arg(species)
  
  HB.list <- list(
    human = c(
      "HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ"
    ),
    rabbit = c(
      "HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ","LOC100344179"
    ),
    mouse = c(
      "Hba-a1","Hba-a2","Hbb-bs","Hbb-bt","Hbb-bh1","Hbb-bh2","Hba-x"
    ),
    rat = c(
      "Hba-a1","Hba-a2","Hbb-b1","Hbb-b2","Hbb-bh1","Hbb-bh2","Hba-x"
    )
  )
  
  return(HB.list[[species]])
}

HB.genes_total <- get_HB_genes("mouse")

hb_present <- HB.genes_total[HB.genes_total %in% rownames(obj_int)]
obj_int[["percent.HB"]] <- PercentageFeatureSet(obj_int, features = hb_present)

# QC_plots
save_plot("P2.1_QC_violin.png",
          VlnPlot(obj_int, features=c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"), ncol=4))



#Filtering (parameterized, available for Django frontend input)#
obj_qc <- subset(obj_int, subset = 
                   nFeature_RNA > 300 &
                   nFeature_RNA < 6000 &
                   percent.mt < 10)

############################################################
# PART3. NORMALIZATION & HVG
############################################################
obj_norm <- NormalizeData(obj_qc)
obj_hvg  <- FindVariableFeatures(obj_norm, selection.method="vst", nfeatures=3000)
obj_scl  <- ScaleData(obj_hvg, features = VariableFeatures(obj_hvg))

############################################################
# PART4. PCA+ElblowPlot_anno
############################################################
obj_pca <- RunPCA(obj_scl, 
                  features = VariableFeatures(obj_scl),
                  npcs=50,
                  verbose = FALSE)

pct <- obj_pca[["pca"]]@stdev / sum(obj_pca[["pca"]]@stdev) * 100 ; cumu <- cumsum(pct)
co1 <- which(cumu >= 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
pc.use <- min(co1, co2, na.rm = TRUE)

Elbowplot <- ElbowPlot(obj_pca)$data %>% ggplot() +
  geom_point(aes(x = dims,y = stdev)) +
  geom_vline(xintercept = pc.use, color = "darkred",linetype = "dashed", size = 1) +
  theme_bw() + labs(title = "Elbow plot: quantitative approach",
                    x = "Principal Components",
                    y = "Standard Deviation")

Elbowplot1 <- ElbowPlot(obj_pca,ndims = min(30, n_pcs))$data %>% ggplot() +
  geom_point(aes(x = dims,y = stdev)) +
  geom_vline(xintercept = pc.use, color = "darkred",linetype = "dashed", size = 1) +
  theme_bw() + labs(title = "Elbow plot: quantitative approach",
                    x = "Principal Components",
                    y = "Standard Deviation")
n_pcs <- length(obj_pca[["pca"]]@stdev)
print(paste("Available PCs:", n_pcs))
Elbowplot_annotated <- Elbowplot1+
  annotate("text", x = pc.use + 3, y = max(obj_pca[["pca"]]@stdev)*0.9,
           label = paste("Use", pc.use, "PCs"), color = "darkred", size = 5)

save_plot("P4.2_ElbowPlot_annotated.png",Elbowplot_annotated,w=4000)
save_plot("P4.1_PCA_plot.png", DimPlot(obj_pca, reduction="pca",dims = 1:2))
save_plot("P4.2_ElbowPlot.png", Elbowplot)

############################################################
# PART5. CLUSTERING + UMAP+ARI
##ARI (Adjusted Rand Index)
##ref_PMID:39294367       
############################################################

#Auxiliary_tool: Cluster_tree assists in selecting resolution values
obj_pca <- FindClusters(object = obj_pca,resolution = c(seq(.1,1.6,.2))) 
res_sel <-clustree(obj_pca@meta.data, prefix = "RNA_snn_res.")

save_plot("P5.1_Resolution_selecting2.png",res_sel,h=3200)

P5_clustering_umap_module <- function(
    seurat_obj,
    pc.use,
    resolution_prefix = "RNA_snn_res.",
    choose_res_fun = choose_resolution,
    save_prefix = "P5",
    check_identical = TRUE
) {
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
  })
  
  message(">>> [P5] Start clustering + UMAP module")
  
  ## -----------------------------
  ## S1. Resolution selection
  ## -----------------------------
  res_choice <- choose_res_fun(seurat_obj)
  main_res <- res_choice$main_resolution
  sub_res  <- res_choice$sub_resolution
  
  message(paste0(
    ">>> [P5] main_res = ", main_res,
    " | sub_res = ", sub_res
  ))
  
  ## -----------------------------
  ## S2. Neighbors + clustering
  ## -----------------------------
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:pc.use)
  
  main_col <- paste0(resolution_prefix, main_res)
  sub_col  <- paste0(resolution_prefix, sub_res, "_sub")
  
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution   = main_res,
    cluster.name = main_col
  )
  
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution   = sub_res,
    cluster.name = sub_col
  )
  
  ## -----------------------------
  ## S3. UMAP (only once)
  ## -----------------------------
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:pc.use)
  
  ## -----------------------------
  ## S4. Safety checks
  ## -----------------------------
  if (!all(c(main_col, sub_col) %in% colnames(seurat_obj@meta.data))) {
    stop(">>> [P5][ERROR] Clustering columns not found in meta.data")
  }
  
  if (check_identical) {
    if (identical(
      as.character(seurat_obj@meta.data[[main_col]]),
      as.character(seurat_obj@meta.data[[sub_col]])
    )) {
      stop(
        ">>> [P5][ERROR] Main and sub clustering are IDENTICAL. ",
        "Sub-resolution clustering is invalid."
      )
    }
  }
  
  ## -----------------------------
  ## S5. Plotting
  ## -----------------------------
  p_main <- DimPlot(
    seurat_obj,
    group.by = main_col,
    label = TRUE
  ) +
    ggtitle(
      paste0(
        "Main clustering | res = ", main_res,
        " | clusters = ",
        length(unique(seurat_obj@meta.data[[main_col]]))
      )
    )
  
  p_sub <- DimPlot(
    seurat_obj,
    group.by = sub_col,
    label = TRUE
  ) +
    ggtitle(
      paste0(
        "Sub clustering | res = ", sub_res,
        " | clusters = ",
        length(unique(seurat_obj@meta.data[[sub_col]]))
      )
    )
  
  save_plot(paste0(save_prefix, ".3_UMAP_cluster_main.png"), p_main)
  save_plot(paste0(save_prefix, ".4_UMAP_cluster_sub.png"),  p_sub)
  
  message(">>> [P5] Done")
  
  return(list(
    seurat_obj  = seurat_obj,
    main_res    = main_res,
    sub_res     = sub_res,
    main_column = main_col,
    sub_column  = sub_col,
    metrics     = res_choice$metrics_table
  ))
}

P5_res <- P5_clustering_umap_module(
  seurat_obj = obj_pca,
  pc.use     = pc.use,
  save_prefix = "P5"
)

obj_pca <- P5_res$seurat_obj

############################################################
# PART6. Marker
############################################################
markers <- FindAllMarkers(obj_pca, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file.path(outdir, "Markers_all_clusters.csv"))

############################################################
# PART7. SingleR_anno
############################################################
ref <- MouseRNAseqData()   #for mouse. should_be_change_to_HPCA_for_homo
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
