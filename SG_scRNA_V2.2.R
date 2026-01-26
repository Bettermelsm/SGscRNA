#############################################
# SG_scRNA_Seq_FLOW_V2.1
# Batch-Processing + Plot Auto-Export + SingleR
# EDIT_by_Sengene/Ribosome
# Update_DATE202601100055
# Proj_for_UCBF
#############################################

############################################################
# E0. ENVIRONMENT
############################################################
options(stringsAsFactors = FALSE)
rm(list = ls()); gc()

############################################################
# E1. BASIC SETTINGS
############################################################
base_dir <- "/mnt/sda1/LSM_scRNAseq/20241118BFlin_UC"
setwd(base_dir)
dir.create("Ana_res3", showWarnings = FALSE)
outdir <- file.path(base_dir, "Ana_res3")

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

#2.7DIY_choose_resolution
choose_resolution <- function(
    seurat_obj,
    res_prefix = "integrated_snn_res.",
    min_main_ARI = 0.9,
    min_sub_ARI  = 0.75,
    min_main_size = 50,
    min_sub_size  = 20,
    smooth_k = 3
) {
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(mclust)
    library(zoo)
  })
  
  #-----------------------------
  # S1. Extract resolution info
  #-----------------------------
  meta <- seurat_obj@meta.data
  res_cols <- grep(paste0("^", res_prefix), colnames(meta), value = TRUE)
  
  resolutions <- sort(as.numeric(sub(res_prefix, "", res_cols)))
  
  #-----------------------------
  # S2. ARI between adjacent resolutions
  #-----------------------------
  ari_df <- map_dfr(
    1:(length(resolutions) - 1),
    function(i) {
      r1 <- paste0(res_prefix, resolutions[i])
      r2 <- paste0(res_prefix, resolutions[i + 1])
      
      data.frame(
        res_low  = resolutions[i],
        res_high = resolutions[i + 1],
        ARI = adjustedRandIndex(
          as.vector(meta[[r1]]),
          as.vector(meta[[r2]])
        )
      )
    }
  )
  
  # smooth ARI to suppress noise
  ari_df$ARI_smooth <- rollmean(
    ari_df$ARI,
    k = smooth_k,
    fill = NA,
    align = "center"
  )
  
  #-----------------------------
  # S3. Cluster statistics
  #-----------------------------
  cluster_stats <- map_dfr(
    resolutions,
    function(r) {
      cl <- meta[[paste0(res_prefix, r)]]
      tab <- table(cl)
      
      data.frame(
        resolution = r,
        n_clusters = length(tab),
        min_size   = min(tab),
        median_size = median(tab)
      )
    }
  )
  
  df <- left_join(
    ari_df,
    cluster_stats,
    by = c("res_high" = "resolution")
  )
  
  #-----------------------------
  # S4. Main resolution (stable plateau)
  #-----------------------------
  main_candidates <- df %>%
    filter(
      ARI_smooth >= min_main_ARI,
      min_size >= min_main_size
    )
  
  if (nrow(main_candidates) == 0) {
    main_res <- NA
  } else {
    # choose median of stable region (NOT the first one)
    main_res <- median(main_candidates$res_high)
  }
  
  #-----------------------------
  # S5. Subcluster resolution (before ARI drop)
  #-----------------------------
  df$ARI_drop <- c(NA, diff(df$ARI_smooth))
  
  sub_candidates <- df %>%
    filter(
      ARI_smooth >= min_sub_ARI,
      min_size >= min_sub_size
    )
  
  if (nrow(sub_candidates) == 0) {
    sub_res <- NA
  } else {
    # last stable resolution before clear drop
    sub_res <- sub_candidates$res_high[
      which.max(sub_candidates$res_high)
    ]
  }
  
  #-----------------------------
  # Output
  #-----------------------------
  return(list(
    main_res = main_res,
    sub_res  = sub_res,
    metrics  = df
  ))
}

#E2.8 sample_clor
sample_color_base <- c("#1D74E7", "#4cb049", "#7DAEE0", "#80d7e1", "#a65527",
                       "#b781d2", "#bf5046", "#b395bd", "#d9e3f5", "#e4cbf2",
                       "#ece7a3", "#EA8379", "#f5cbe1", "#ffb7ba", "#003366",
                       "#336699", "#FF6699", "#CC3399", "#996633", "#FF9933",  
                       "#339966", "#33CCFF")

############################################################
# PART1. DATA IMPORT (BATCH) & MERGING
############################################################
sample_folders <- c("NC1","NC3","NC4",
                    "UC1","UC2","UC3","UC4")

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
    grepl("^NC", smp)  ~ "NC",
    grepl("^UC", smp) ~ "UC",
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

#When calculating QC indicators, must specify assay = "RNA" forcibly
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

HB.genes_total <- get_HB_genes("human")

hb_present <- HB.genes_total[HB.genes_total %in% rownames(obj_int)]
obj_int[["percent.HB"]] <- PercentageFeatureSet(obj_int, features = hb_present)

# QC_plots
save_plot("P2.1_QC_violin.png",
          VlnPlot(obj_int, features=c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"), ncol=4))

############################################################
# PART3.PRE-PROCESSING & SCALE
############################################################
DefaultAssay(obj_qc) <- "integrated"
obj_scl  <- ScaleData(obj_hvg, features = VariableFeatures(obj_hvg))

############################################################
# PART4. PCA+ElblowPlot_anno
############################################################
DefaultAssay(obj_scl) <- "integrated"
obj_pca <- RunPCA(obj_scl,npcs=50,verbose = FALSE)

#Inflection point identification
pct <- obj_pca[["pca"]]@stdev / sum(obj_pca[["pca"]]@stdev) * 100 
cumu <- cumsum(pct)
co1 <- which(cumu >= 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
pc.use <- min(co1, co2, na.rm = TRUE)

#if(pc.use < 10) { pc.use <- 20 }

n_pcs <- length(obj_pca[["pca"]]@stdev)
print(paste("Available PCs:", n_pcs))

Elbowplot <- ElbowPlot(obj_pca,ndims = min(30, n_pcs))$data %>% ggplot() +
  geom_point(aes(x = dims,y = stdev)) +
  geom_vline(xintercept = pc.use, color = "darkred",linetype = "dashed", linewdith = 1) +
  annotate("text", x = pc.use + 3, y = max(obj_pca[["pca"]]@stdev)*0.9,
           label = paste("Use", pc.use, "PCs"), color = "darkred", size = 5)+
  theme_bw() + labs(title = "Elbow plot: quantitative approach",
                    x = "Principal Components",
                    y = "Standard Deviation")


save_plot("P4.1_PCA_plot.png", DimPlot(obj_pca, reduction="pca",dims = 1:2))
save_plot("P4.2_ElbowPlot_annotated.png",Elbowplot,w=2000)

message(paste0(">>> [PART4] Selected PC count: ", pc.use))

############################################################
# PART5. CLUSTERING + UMAP + TSNE + ARI
##ARI (Adjusted Rand Index,ref_PMID:39294367)
############################################################
DefaultAssay(obj_pca) <- "integrated"

obj_pca <- ScaleData(obj_pca, verbose = FALSE)
obj_pca <- RunPCA(obj_pca, npcs = 30, verbose = FALSE)

# P5.1.Precomputed neighbors (for clustre selection resolution)
obj_pca <- FindNeighbors(obj_pca, reduction = "pca", dims = 1:pc.use)

# P5.2Auxiliary_tool: Cluster_tree assists in selecting resolution values
obj_pca <- FindClusters(object = obj_pca,resolution = c(seq(0.1,1.6,0.2))) 
res_sel <-clustree(obj_pca@meta.data, prefix = "integrated_snn_res.")
save_plot("P5.2_Resolution_selecting.png",res_sel,h=3200,w=2000)

# P5.3Define the main function of clustering + dimension reduction 
P5_clustering_module <- function(
    seurat_obj,
    pc.use,
    resolution_prefix = "integrated_snn_res.",
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
  
# P5.4Clustering + UMAP module 
  ## -----------------------------
  ## M1.Resolution selection
  ## -----------------------------
  res_choice <- choose_res_fun(seurat_obj, res_prefix = resolution_prefix)
  
  main_res <- res_choice$main_res
  sub_res  <- res_choice$sub_res
  
  if (is.null(main_res) || is.na(main_res) || length(main_res) == 0) {
    warning(">>> [P5] Auto-selection failed or returned NA. Defaulting main_res to 0.5")
    main_res <- 0.5
  }
  if (is.null(sub_res) || is.na(sub_res) || length(sub_res) == 0) {
    sub_res <- 0.8
  }
  
  message(paste0(">>> [P5] main_res = ", main_res, " | sub_res = ", sub_res))
  
  ## -----------------------------
  ## M2.Neighbors + clustering
  ## -----------------------------
 
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:pc.use)
  
  main_col <- paste0(resolution_prefix, main_res)
  sub_col  <- paste0(resolution_prefix, sub_res)
  
  #Run selected resolution
  seurat_obj <- FindClusters(seurat_obj, resolution = main_res, cluster.name = main_col, verbose = FALSE)
  
  if (sub_res != main_res) {
    seurat_obj <- FindClusters(seurat_obj, resolution = sub_res, cluster.name = sub_col, verbose = FALSE)
  }
  
  ## -----------------------------
  ## M3.UMAP + t-SNE
  ## -----------------------------
  message(">>> [P5] Running UMAP...")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:pc.use, reduction = "pca")
  
  message(">>> [P5] Running t-SNE...")
  perp_val <- min(30, floor((ncol(seurat_obj)-1)/3))
  seurat_obj <- RunTSNE(seurat_obj, dims = 1:pc.use, reduction = "pca", 
                        perplexity = perp_val, check_duplicates = FALSE)
  
  ## -----------------------------
  ## M4.Safety checks
  ## -----------------------------
  if (!main_col %in% colnames(seurat_obj@meta.data)) {
    #Try to fall back to the default column name lookup (to prevent Seurat version differences)
    default_col <- paste0("integrated_snn_res.", main_res)
    if(default_col %in% colnames(seurat_obj@meta.data)) {
      main_col <- default_col
    } else {
      stop(">>> [P5][ERROR] Clustering columns not found in meta.data")
    }
  }
  
  ## -----------------------------
  ## M5.Plotting
  ## -----------------------------
  # UMAP Plot
  p_main_umap <- DimPlot(seurat_obj, group.by = main_col, reduction = "umap", label = TRUE) +
    ggtitle(paste0("UMAP | res = ", main_res))
  
  # t-SNE Plot
  p_main_tsne <- DimPlot(seurat_obj, group.by = main_col, reduction = "tsne", label = TRUE) +
    ggtitle(paste0("t-SNE | res = ", main_res))
  
  save_plot(paste0(save_prefix, ".4_UMAP_cluster.png"), p_main_umap)
  save_plot(paste0(save_prefix, ".5_tSNE_cluster.png"), p_main_tsne)
  
  message(">>> [P5] Done")
  
  return(list(
    seurat_obj  = seurat_obj,
    main_res    = main_res,
    sub_res     = sub_res,
    main_column = main_col,
    metrics     = res_choice$metrics
  ))
}

P5_res <- P5_clustering_module(
  seurat_obj = obj_pca,
  pc.use     = pc.use,
  save_prefix = "P5",
  resolution_prefix = "integrated_snn_res."
)

obj_pca <- P5_res$seurat_obj

#####P5.6Preliminary visualization view results#####
DimPlot(obj_pca, label = T, reduction = "tsne")
DimPlot(obj_pca_clean, label = T, reduction = "umap")
DimPlot(obj_pca_clean, label = F, reduction = "tsne", group.by = "group")
DimPlot(obj_pca_clean, label = F, reduction = "tsne", group.by = "orig.ident")

# P5.7get the original object and remove UC4
# ---------------------------------------------------------

obj_pca_clean <- subset(obj_pca, subset = orig.ident != "UC4")
DefaultAssay(obj_pca_clean) <- "RNA"

# P5.8re-normalize, normalize and run PCA (required steps)
# ---------------------------------------------------------
obj_pca_clean <- NormalizeData(obj_pca_clean)
obj_pca_clean <- FindVariableFeatures(obj_pca_clean)
obj_pca_clean <- ScaleData(obj_pca_clean)
# 50 pcs are calculated to provide enough data for screening algorithm
obj_pca_clean <- RunPCA(obj_pca_clean, npcs = 50) 

pct <- obj_pca_clean[["pca"]]@stdev / sum(obj_pca_clean[["pca"]]@stdev) * 100 
cumu <- cumsum(pct)
co1 <- which(cumu >= 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
pc.use_auto <- min(co1, co2, na.rm = TRUE)

if(is.na(pc.use_auto) || pc.use_auto < 10) { 
  message("Calculated pc.use is too small or NA, defaulting to 15.")
  pc.use_auto <- 15 
}

message(paste(">>> Auto-selected PCs:", pc.use_auto))

# ---------------------------------------------------------
obj_pca_clean <- RunUMAP(obj_pca_clean, dims = 1:pc.use_auto)
obj_pca_clean <- RunTSNE(obj_pca_clean, dims = 1:pc.use_auto)

# ---------------------------------------------------------
Ptem6 <- DimPlot(obj_pca_clean, label = T, reduction = "tsne") + ggtitle(paste0("t-SNE (No UC4, PCs=", pc.use_auto, ")"))
Ptem7 <- DimPlot(obj_pca_clean, label = T, reduction = "umap") + ggtitle(paste0("UMAP (No UC4, PCs=", pc.use_auto, ")"))
Ptem8 <- DimPlot(obj_pca_clean, label = F, reduction = "tsne", group.by = "group")
Ptem9 <- DimPlot(obj_pca_clean, label = F, reduction = "tsne", group.by = "orig.ident")

save_plot("P5.8.1t-SNE_NoUC4.png", Ptem6)
save_plot("P5.8.2umap_NoUC4.png", Ptem7)
save_plot("P5.8.3Group_tsne_NoUC4.png",Ptem8)
save_plot("P5.8.4Ident_NoUC4.png", Ptem9)

#P5.9consolidate the data layer to facilitate subsequent mapping
obj_pca_clean <- JoinLayers(obj_pca_clean)

#umap_before_anno
my_palette_func <- colorRampPalette(sample_color_base)
n_colors_needed <- length(unique(Idents(obj_pca)))
final_colors <- my_palette_func(n_colors_needed)

P5.9 <- DimPlot(obj_pca, reduction = "umap", label = TRUE, label.size = 4.5) +
  scale_color_manual(values = final_colors) +
  theme(text = element_text(family = "Times New Roman")) 
save_plot("P5.9umap_before_anno.png",P5.9)

P5.10 <- DimPlot(obj_pca_clean, reduction = "umap", label = TRUE, label.size = 4.5) +
  scale_color_manual(values = final_colors) +
  theme(text = element_text(family = "Times New Roman")) 
save_plot("P5.10umap_removeUC4.png",P5.10)

##
obj_pca <- obj_pca_clean

############################################################
# PART6. Marker
############################################################
markers <- FindAllMarkers(obj_pca, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file.path(outdir, "Markers_all_clusters.csv"))

############################################################
# PART7. SINGLE-R ANNOTATION (SPECIES-FLEXIBLE)
############################################################

# Define supported references (extendable)
get_SingleR_reference <- function(
    species = c("human", "mouse", "rat", "rabbit"),
    ref_type = c("HPCA", "Blueprint", "MouseRNAseq", "DICE", "Monaco")
) {
  species <- match.arg(species)
  ref_type <- match.arg(ref_type)
  
  # Check required packages
  if (!requireNamespace("celldex", quietly = TRUE)) {
    stop("Package 'celldex' required for reference data. Install via:\n",
         "BiocManager::install('celldex')")
  }
  if (!requireNamespace("SingleR", quietly = TRUE)) {
    stop("Package 'SingleR' required. Install via:\n",
         "BiocManager::install('SingleR')")
  }
  
  message(">>> [PART7] Loading reference for species = '", species, 
          "', ref_type = '", ref_type, "'")
  
  ref <- switch(
    paste0(species, "_", ref_type),
    
    # ─── Human ───────────────────────────────────────────────
    "human_HPCA"        = celldex::HumanPrimaryCellAtlasData(),
    "human_Blueprint"   = celldex::BlueprintEncodeData(),
    "human_DICE"        = celldex::DICEData(),
    "human_Monaco"      = celldex::MonacoImmuneData(),
    
    # ─── Mouse ───────────────────────────────────────────────
    "mouse_MouseRNAseq" = celldex::MouseRNAseqData(),  # default mouse
    "mouse_Blueprint"   = celldex::BlueprintEncodeData(species = "mouse"),
    
    # ─── Rat / Rabbit (fallback to closest or custom) ────────
    "rat_MouseRNAseq"   = {
      warning("No native rat reference in celldex; using mouse as proxy (caution!).")
      celldex::MouseRNAseqData()
    },
    "rabbit_MouseRNAseq" = {
      warning("No native rabbit reference; using mouse as proxy (caution!).")
      celldex::MouseRNAseqData()
    },
    
    # ─── Defaults ────────────────────────────────────────────
    "human"  = celldex::HumanPrimaryCellAtlasData(),   # HPCA = default human
    "mouse"  = celldex::MouseRNAseqData(),             # default mouse
    "rat"    = celldex::MouseRNAseqData(),             # proxy
    "rabbit" = celldex::MouseRNAseqData(),             # proxy
    
    # ─── Error ───────────────────────────────────────────────
    stop("Unsupported species/ref combination: ", species, "/", ref_type)
  )
  
  # Extract label column (common names across refs)
  label_col <- switch(
    class(ref)[1],
    "SummarizedExperiment" = {
      if ("label.main" %in% colnames(colData(ref))) "label.main"
      else if ("label.fine" %in% colnames(colData(ref))) "label.fine"
      else if ("label" %in% colnames(colData(ref))) "label"
      else stop("No known label column in reference.")
    },
    "data.frame" = {
      if ("label.main" %in% names(ref)) "label.main" else "label"
    },
    stop("Unknown reference object class.")
  )
  
  return(list(ref = ref, label_col = label_col))
}

# ───────────────────────────────────────────────────────────────
# USER CONFIGURATION ZONE (EDIT HERE)
# ───────────────────────────────────────────────────────────────
#************************************************************************
#*#**********************************************************************
SPECIES     <- "human" # ← type_here：human / mouse / rat / rabbit
REF_TYPE    <- "HPCA"  # ← human: HPCA/Blueprint/DICE/Monaco; mouse: MouseRNAseq/Blueprint
ASSAY_USED  <- "RNA"   # do not change.["RNA"（log-normalized data）]
# ───────────────────────────────────────────────────────────────
#************************************************************************
#*#**********************************************************************

# Load reference
ref_info <- get_SingleR_reference(species = SPECIES, ref_type = REF_TYPE)
ref      <- ref_info$ref
label_col <- ref_info$label_col

# ───────────────────────────────────────────────────────
# SEURAT v5.3+ → SCE CONVERSION (Robust Method)
# ───────────────────────────────────────────────────────
message(">>> [PART7] Converting to SingleCellExperiment (Seurat v5 fix)")

# merged Layers to one whole data
if (inherits(obj_pca[["RNA"]], "Assay5")) {
  message("Detected Assay5. Joining layers to create unified 'data' matrix...")
  obj_pca <- JoinLayers(obj_pca, assay = "RNA")
}

# get the Log-normalized matrix (Using LayerData in Seurat V5+)
# make sure that the name of layer is "data" after JoinLayers
norm_mat <- LayerData(obj_pca, assay = "RNA", layer = "data")

message("Extracted matrix dimensions: ", nrow(norm_mat), " genes × ", ncol(norm_mat), " cells")

#Check and convert sparse matrices (single may process dense matrices faster, but pay attention to memory)
#If the memory is insufficient (error: cannot allocate vector of size...), # if block and directly use the sparse matrix
if (inherits(norm_mat, "dgCMatrix")) {
  mem_est <- object.size(norm_mat) * 3 / 1024^2 # Estimate MB after conversion to dense
  message("Converting sparse matrix to dense... (Est. RAM usage: ", round(mem_est, 1), " MB)")
  
  # If the matrix is very large (>10gb), skip as.matrix. The current version of single supports sparse matrix
  if(mem_est < 10000) { 
    norm_mat <- as.matrix(norm_mat)
  } else {
    message("Matrix is too large, keeping as sparse matrix to save RAM.")
  }
}

# establish SCE object
library(SingleCellExperiment)
obj_sce <- SingleCellExperiment(
  assays = list(logcounts = norm_mat),
  colData = obj_pca@meta.data
)

rownames(obj_sce) <- rownames(norm_mat)

message(">>> SCE conversion successful.")

# Proceed to SingleR

pred <- SingleR(
        test = obj_sce,
        ref  = ref,
        labels = ref[[label_col]],
        assay.type.test = "logcounts"
)

# Store results
anno_col_name <- paste0("SingleR_", SPECIES, "_", REF_TYPE)
obj_pca[[anno_col_name]] <- pred$labels

# Log summary
cat(">>> [PART7] Annotation completed:\n")
cat("    - Species:", SPECIES, "\n")
cat("    - Reference:", REF_TYPE, "\n")
cat("    - Column name:", anno_col_name, "\n")
cat("    - Unique cell types:", length(unique(pred$labels)), "\n")

# Plot
p_singleR <- DimPlot(
  obj_pca,
  group.by = anno_col_name,
  label = TRUE,
  label.size = 3,
  repel = TRUE
) +
  ggtitle(paste("SingleR Annotation:", SPECIES, "-", REF_TYPE)) +
  theme(plot.title = element_text(hjust = 0.5))

save_plot("P7.SingleR_annotation.png", p_singleR, w = 3000, h = 2000)

# Optional: Save prediction scores for QC
pred_df <- data.frame(
  cell = rownames(pred), 
  label = pred$labels,
  delta = pred$delta.next
)
write.csv(pred_df, file.path(outdir, "SingleR_prediction_scores.csv"), row.names = FALSE)

############################################################
# PART8. SAVE OBJECT
############################################################
saveRDS(obj_pca, file.path(outdir, "Final_Integrated_Object.rds"))

message("==== Pipeline finished ====")


############################################################
# PART 9. CELL PROPORTION & STATISTICAL COMPARISON (FIXED)
############################################################
message(">>> [PART 9] Running Fixed Proportion Analysis...")

# 9.1 样本级比例计算 (显式定义 Cluster 列名)
prop_sample <- obj_pca@meta.data %>%
  group_by(orig.ident, Cluster = Idents(obj_pca)) %>% # 直接在这里命名 Cluster
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup()

# 合并组别信息 (使用 dplyr::select 避免冲突)
sample_info <- obj_pca@meta.data %>% 
  dplyr::select(orig.ident, group) %>% 
  distinct()

prop_sample <- left_join(prop_sample, sample_info, by = "orig.ident")

# 9.2 全局比例计算 (用于堆叠图)
prop_global <- prop_sample %>%
  group_by(group, Cluster) %>%
  summarise(mean_perc = mean(percentage), .groups = "drop")

# 9.3 绘图 P9.1: 堆叠柱状图
p9.1 <- ggplot(prop_global, aes(x = group, y = mean_perc, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.6, color = "white") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  scale_fill_manual(values = final_colors) +
  labs(title = "Global Composition", x = "", y = "Relative Abundance") +
  theme_classic() + theme(text = element_text(family = "Times New Roman"))

# 9.4 绘图 P9.2: 样本级箱线图 (修正版)
p9.2 <- ggplot(prop_sample, aes(x = group, y = percentage, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 1.2) +
  facet_wrap(~Cluster, scales = "free_y", ncol = 4) + # 此时 Cluster 变量必然存在
  scale_fill_manual(values = c("NC" = "#4cb049", "UC" = "#bf5046")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
  theme_bw() + theme(text = element_text(family = "Times New Roman"), legend.position = "none")

# 9.5 最终大图整合 (Figure 1 & 2)
library(patchwork)
fig_final_1 <- ((VlnPlot(obj_int, features=c("nFeature_RNA","percent.mt"), ncol=2) + NoLegend()) /
                  (Elbowplot | res_sel) / (P5.9 | Ptem8)) + 
  plot_annotation(tag_levels = 'A') & theme(text = element_text(family = "Times New Roman"))

fig_final_2 <- ((p_singleR) / (p9.1 | p9.2) + plot_layout(heights = c(1.2, 1))) + 
  plot_annotation(tag_levels = 'A') & theme(text = element_text(family = "Times New Roman"))

ggsave(file.path(outdir, "SCI_Figure1_Clustering.png"), fig_final_1, width = 12, height = 16)
ggsave(file.path(outdir, "SCI_Figure2_Quantitation.png"), fig_final_2, width = 14, height = 14)

message("==== PART 9 Analysis Finished: All plots generated successfully ====")


############################################################
# PART 10. SCI FIGURE LEGEND GENERATION (BILINGUAL)
############################################################

# 提取当前环境中的关键参数
pc_val   <- pc.use_auto
res_val  <- P5_res$main_res
n_cells  <- ncol(obj_pca)
n_cl     <- length(unique(Idents(obj_pca)))

# --- Figure 1: 流程与质控 ---
f1_legend <- paste0(
  "Figure 1. Single-cell landscape and clustering quality control of the UCBF dataset.\n",
  "（图 1. UCBF 数据集的单细胞图谱与质量控制。）\n\n",
  "(A) Quality Control: Initial quality control (QC) was performed to filter out low-quality cells. Violin plots show the distribution of gene counts (nFeature), UMI counts (nCount), and mitochondrial gene percentages (percent.mt) per cell. ",
  "（A. 初步质控过滤了低质量细胞。小提琴图展示了每个细胞的基因数、UMI 数和线粒体基因比例。）\n",
  "(B) PC Selection: An Elbow plot showing the quantitative approach for PC selection; the red dashed line indicates the selected ", pc_val, " PCs. ",
  "（B. 碎石图展示了确定主成分数的定量方法；红虚线表示选定的 ", pc_val, " 个主成分。）\n",
  "(C) Resolution Optimization: Clustering tree (Clustree) visualization illustrating the relationship between clusters at different resolutions, supporting the selection of resolution ", res_val, ". ",
  "（C. 聚类树状图展示了不同分辨率下集群间的演变关系，支持了分辨率参数 ", res_val, " 的选择。）\n",
  "(D-E) Global Clustering: UMAP and t-SNE visualizations identifying diverse cell populations, with NC and UC groups showing consistent sample integration. ",
  "（D-E. UMAP 和 t-SNE 可视化识别出多样的细胞群，NC 和 UC 组显示出良好的样本整合效果。）"
)

# --- Figure 2: 鉴定与比例定量 ---
f2_legend <- paste0(
  "Figure 2. Cell type identity and compositional alterations in UC.\n",
  "（图 2. 溃疡性结肠炎中的细胞身份鉴定与比例改变。）\n\n",
  "(A) Cell Annotation: Automated cell type identification performed by SingleR using the ", REF_TYPE, " reference dataset, identifying ", n_cl, " major cell lineages. ",
  "（A. 细胞注释：利用 SingleR 结合 ", REF_TYPE, " 参考数据集进行的自动化细胞类型鉴定，识别出 ", n_cl, " 类主要的细胞系谱。）\n",
  "(B) Compositional Overview: A stacked bar chart showing the relative abundance of each cell lineage across the NC and UC groups. ",
  "（B. 比例总览：堆叠柱状图展示了 NC 组和 UC 组中各细胞系谱的相对丰度。）\n",
  "(C) Statistical Comparison: Boxplots illustrating the percentage of cell clusters in individual samples. P-values were calculated using the Wilcoxon rank-sum test. ",
  "（C. 统计比较：箱线图展示了单个样本中各细胞群的百分比。P 值通过 Wilcoxon 秩和检验计算。）"
)

# --- 打印输出到控制台 ---
cat("\n\n==================== SCI FIGURE LEGENDS ====================\n\n")
cat(f1_legend, "\n\n")
cat("------------------------------------------------------------\n\n")
cat(f2_legend, "\n\n")
cat("============================================================\n")

# --- 同时导出为文本文件 ---
writeLines(c(f1_legend, "\n\n", f2_legend), file.path(outdir, "Figure_Legends_Bilingual.txt"))
