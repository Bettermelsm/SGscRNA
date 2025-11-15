#####SG_scRNA_Seq_FLOW
####EDIT_by_Ribosome_Sengene
##Update_DATE202510151524

#######1.Creat_Seurta_project
library(dplyr)
library(Seurat)
library(patchwork)
library(png)
library(ggplot2) 
library(SingleR)
library(celldex)

#1.Data_get

setwd("/mnt/sda1/LSM_scRNAseq/20241118BFlin_UC")
#dir.create("Ana_res")
#pbmc.data <- Read10X(data.dir =DIRNC1,gene.column = 1)
#pbmc <- CreateSeuratObject(counts = pbmc.data,project = "pbmc3k",min.cells = 3,min.features = 200)
#counts_matrix <- GetAssayData(pbmc, assay = "RNA", slot = "counts") %>% data.frame()

##1.1set_dir_path
sample_folders <- c("NC1", "NC3", "NC4", "UC1", "UC2", "UC3","UC4")
data_dir <- "/mnt/sda1/LSM_scRNAseq/20241118BFlin_UC/All_filter_data_raw" 

seurat_list <- list()# a_new_list_to_store_Seurat_objects

## 1.2Read each sample in a loop
for (folder in sample_folders) {
  sample_path <- file.path(data_dir, folder)
  data_matrix <- Read10X(data.dir = sample_path, gene.column = 1)
  seurat_obj <- CreateSeuratObject(counts = data_matrix, 
                                   project = folder, 
                                   min.cells = 3, 
                                   min.features = 200)  
  seurat_obj$sample <- folder
  # 添加组信息：NC or UC
  seurat_obj$group <- ifelse(grepl("^NC", folder), "NC", "UC")
  seurat_list[[folder]] <- seurat_obj
}

# Normalize_each_sample
seurat_list <- lapply(seurat_list, NormalizeData, normalization.method = "LogNormalize", scale.factor = 10000)

# find_anchors_and_interateData-Sample_as_batch
anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:20)
pbmc_integrated <- IntegrateData(anchorset = anchors, dims = 1:20)

# swfit_to_integrated assay
DefaultAssay(pbmc_integrated) <- "integrated"

pbmc_flt <- pbmc_integrated

#2.QC
##2.1Check_the_ratio_of_mitochondrial_cells_to_red_blood_cells
###Remove_cells_with_a_high_proportion_of_mitochondrial_gene_expression_levels
###Cells_undergoing_autophagic and apoptosis typically exhibit high expression of mitochondrial genes.Such cells can be removed by using this indicator. 

# 关键修改：切换到RNA assay来计算QC指标
DefaultAssay(pbmc_flt) <- "RNA"

pbmc_flt [["percent.mt"]]<-PercentageFeatureSet(pbmc_flt, pattern = "^mt-") 

HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") # Common RBC genes in human blood.
#HB.genes_total <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", 
#                    "HBG2", "HBM", "HBQ1", "HBZ", "LOC100344179") # Common RBC genes in rabbit blood.
HB_m <- match(HB.genes_total,rownames(pbmc_flt@assays$RNA))

HB.genes <- rownames(pbmc_flt@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
pbmc_flt[["percent.HB"]]<-PercentageFeatureSet(pbmc_flt,features=HB.genes)

head(pbmc_flt@meta.data)[,c(2,3,4,5)]

##2.2Feature、count、线粒体基因、红细胞基因占比可视化
##nCount_RNA就是UMI的读数，nFeature_RNA就是Gene

DIR00 <- c("/mnt/sda1/LSM_scRNAseq/20241118BFlin_UC")
DIR01 <- paste0(DIR00,"/Ana_res")
setwd(DIR01)
png("P2.2QC1.png",width = 2700,height = 1600,res=200)
VlnPlot(pbmc_flt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 4)
dev.off()

##2.3几个指标之间的相关性
##FeatureScatter通常用于可视化特征-特征关系
plot2.3.1 <- FeatureScatter(pbmc_flt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2.3.2 <- FeatureScatter(pbmc_flt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2.3.3 <- FeatureScatter(pbmc_flt, feature1 = "nCount_RNA", feature2 = "percent.HB")

# 使用patchwork组合图形（水平排列）
plot2.3 <- plot2.3.1 + plot2.3.2 + plot2.3.3 + plot_layout(nrow = 1)  # 确保三个图排成一行

# 保存组合图形
png('P2.3FeatureScatter_plot.png', 
    width = 2700, height = 1600, res = 300)
print(plot2.3)  # 显式打印组合图
dev.off()


##2.4数据均一化与标准化
pbmc_flt<-subset(pbmc_flt, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) 

###2.4.1数据归一化,去除样本/细胞效应
pbmc_flt1<-NormalizeData(pbmc_flt,normalization.method = "LogNormalize", scale.factor = 10000)


#Warning messages:
#  1: In scale_x_log10() :
#  log-10 transformation introduced infinite values.
#2: Removed 1555 rows containing missing values or values outside the scale range (`geom_point()`). 
#3: In scale_x_log10() :
# log-10 transformation introduced infinite values.
#4: Removed 1555 rows containing missing values or values outside the scale range (`geom_point()`). 
#在画P2.4时有如上提示，可能是有未过滤完的低表达量基因，故再去除一次


# 方法1：使用 LayerData() 获取数据并添加伪计数
#tryCatch({
# 获取 counts 层数据
#  counts_layer <- LayerData(pbmc_flt1, layer = "counts")

# 添加伪计数解决0值问题
#  counts_layer[counts_layer == 0] <- 1

# 更新对象
#  pbmc_flt1 <- SetAssayData(pbmc_flt1, layer = "counts", new.data = counts_layer)
#}, error = function(e) {
#  message("使用 LayerData 方法失败，尝试备选方案: ", e$message)
#})

# 方法2：如果上述方法失败，直接过滤零表达基因
#if (!exists("counts_layer")) {
# 计算每个基因的表达细胞数
#  n_cells_expressed <- Matrix::rowSums(GetAssayData(pbmc_flt1, assay = "RNA") > 0)

# 保留在至少1个细胞中表达的基因
#  genes_to_keep <- names(n_cells_expressed[n_cells_expressed > 0])
#  pbmc_flt1 <- subset(pbmc_flt1, features = genes_to_keep)
#}


###特征选择：高变基因，FindVariableFeatures，默认方法”vst”
pbmc_flt2<- FindVariableFeatures(pbmc_flt1,selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc_flt2),10)# 提取差异最大的 top10 基因

plot2.4.1 <- VariableFeaturePlot(pbmc_flt2,log = NULL)
plot2.4.2 <- LabelPoints(plot = plot2.4.1, points = top10, repel = TRUE)
#plot2.4 <- CombinePlots(plots = list(plot2.4.1, plot2.4.2),legend="bottom")
plot2.4 <- plot2.4.1+ plot2.4.2+plot_layout(nrow = 1)& theme(legend.position = "bottom") 
png('P2.4Top10_variable_genes.png',width = 3500,height = 1900,res=300)
plot2.4
dev.off()

###2.4.2数据标准化
pbmc_flt3<-ScaleData(pbmc_flt2)

#all.genes <- rownames(pbmc)
#pbmc <- ScaleData(pbmc, features = all.genes)

#3.降维分析
##3.1特征提取：PCA降维

pbmc_flt4<-RunPCA(pbmc_flt3,features = VariableFeatures(object=pbmc_flt3))
print(pbmc_flt4[["pca"]], dims = 1:5, nfeatures = 5) 

##3.2结果查看：每个细胞在PC轴上的坐标
head(pbmc_flt4@reductions$pca@cell.embeddings)

##3.3结果查看：每个细胞在PC轴上的坐标每个基因对每个PC轴的贡献度（loading值）
head(pbmc_flt4@reductions$pca@feature.loadings)

##3.4结果查看：loading值分析
t(Loadings(object =pbmc_flt4[["pca"]])[1:5,1:5])
t(Loadings(object =pbmc_flt4, reduction = "pca")[1:5,1:5])

new.loadings <- Loadings(object = pbmc_flt4[["pca"]])
new.loadings <- new.loadings + 0.01
Loadings(object = pbmc_flt4[["pca"]]) <- new.loadings
plot3.4 <- VizDimLoadings(pbmc_flt4)
png("P3.4VizDimloadings_plot.png",width = 2000,height = 4500,res=300)
plot3.4
dev.off()

#查看PCA结果的另一种写法
print(pbmc_flt4[["pca"]], dims = 1:5, nfeatures = 5)

##3.5PCA结果可视化
plot3.5 <- DimPlot(pbmc_flt4, reduction = "pca")
png("P3.5PCA_1and2_plot.png",width = 1800,height = 1000,res=300)
plot3.5
dev.off()

##3.6PCA热图分析

png("P3.6PCA_1pheatmap15.png",width = 4000,height = 3600,res=300)
DimHeatmap(pbmc_flt4, dims = 1:20, cells = 500, balanced = TRUE)
dev.off()


##3.7PCA维度选择
pbmc_flt4_pca <- JackStraw(pbmc_flt4, num.replicate = 100)
pbmc_flt4_pca <- ScoreJackStraw(pbmc_flt4_pca, dims = 1:20)
plot3.7.1<-JackStrawPlot(pbmc_flt4_pca, dims = 1:20)
plot3.7.2<-ElbowPlot(pbmc_flt4_pca)
#plot3.7 <- plot3.7.1+ plot3.7.2+plot_layout(nrow = 1)& theme(legend.position = "bottom") 
plot3.7 <- plot3.7.1 + plot3.7.2 + plot_layout(nrow = 1) & theme(plot.title = element_text(hjust = 0.5))
png("P3.7PCA_dimetion_selection.png",width = 4000,height = 3600,res=300)
plot3.7
dev.off()

############################################[optional_in_this_project]
#4.聚类
# Cluster the cells 
#Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm. First calculate k-nearest neighbors and construct the SNN graph. Then optimize the modularity function to determine clusters. For a full description of the algorithms, see Waltman and van Eck (2013) The European Physical Journal B. Thanks to Nigel Delaney (evolvedmicrobe@github) for the rewrite of the Java modularity optimizer code in Rcpp!
pbmc_flt4 <- FindNeighbors(pbmc_flt4, dims = 1:20)
pbmc_flt4 <- FindClusters(pbmc_flt4, resolution = 0.8)#resolution 越大，群体越细。通常 0.6~1.0 是一个比较好的范围。

#table(pbmc@active.ident) # 查看每一类有多少个细胞

# 提取某一类细胞。
#head(subset(as.data.frame(pbmc@active.ident),pbmc@active.ident=="2"))

##4.1提取某一cluster细胞。
#subpbmc<-subset(x = pbmc,idents="2")

#head(WhichCells(pbmc,idents="2"))
#head(Idents(pbmc), 5)

#提取部分细胞

#head(colnames(pbmc@assays$RNA@counts)[1:30])
#subbset<-subset(x=pbmc,cells=colnames(pbmc@assays$RNA@counts)[1:30])

##4.2构建系统发育分析（Phylogenetic Analysis of Identity Classes）
#pbmc_flt4_tree<-BuildClusterTree(pbmc_flt4)
#Tool(object = pbmc_flt4, slot = 'BuildClusterTree')
#PlotClusterTree(pbmc_flt4)
#pbmc<-CalculateBarcodeInflections(pbmc_flt4)
#SubsetByBarcodeInflections(pbmc_flt4)

############################################



##4.3可视化降维
###4.3.1UMAP
#非线性降维——这个目的是为了可视化，而不是特征提取（PCA），虽然它也可以用来做特征提取。

#pbmc_flt5<- RunUMAP(pbmc_flt4,dims=1:10)
#pbmc_flt5<- RunUMAP(pbmc_flt4,dims=1:20)

#因为分群不理想，所以作出新调整
pbmc_flt5 <- RunUMAP(pbmc_flt4, dims = 1:20, 
                     n.neighbors = 20,        # 调整邻居数
                     min.dist = 0.15)          # 调整最小距离，原为0.3
head(pbmc_flt5@reductions$umap@cell.embeddings)

###4.3.2tSNE
#pbmc_flt6 <- RunTSNE(pbmc_flt4, dims = 1:10)
#pbmc_flt6 <- RunTSNE(pbmc_flt4, dims = 1:20)
pbmc_flt6 <- RunTSNE(pbmc_flt4, dims = 1:20, perplexity = 40)  # 调整perplexity,原为30
head(pbmc_flt6@reductions$tsne@cell.embeddings)

#比较一下两个可视化的结果
plot4.3.1<-DimPlot(pbmc_flt5, reduction = "umap",label = TRUE)
plot4.3.2<-DimPlot(pbmc_flt6, reduction = "tsne",label = TRUE)
plot4.3 <- plot4.3.1 + plot4.3.2 + plot_layout(nrow = 1) & theme(plot.title = element_text(hjust = 0.5))
png("P4.3UMAP_tSNE3.png",width = 5700,height = 2600,res=300)
plot4.3
dev.off()

pbmc_flt8_data<-as.data.frame(Idents(pbmc_flt6))
pbmc_flt8_data<-cbind(Cell=rownames(pbmc_flt8_data),pbmc_flt8_data)

#5.Annotation_with_singleR
##5.1 确保对象处理完整

pbmc_anno <- pbmc_flt6
pbmc_anno[["umap"]] <- pbmc_flt5[["umap"]]
pbmc_flt9 <- pbmc_anno
pbmc_flt9 <- NormalizeData(pbmc_flt9)
pbmc_flt9 <- JoinLayers(pbmc_flt9)

##5.2 提取表达矩阵
pbmc_flt9_matrix <- LayerData(pbmc_flt9, layer = "data")

##5.3 获取参考数据集（修改版本）
get_reference <- function() {
  # 尝试不同的加载方式
  tryCatch({
    # 方式1：直接加载，不设置localHub
    ref <- HumanPrimaryCellAtlasData()
    message("成功加载 HumanPrimaryCellAtlasData")
    return(ref)
  }, error = function(e) {
    message("方式1失败: ", e$message)
    
    tryCatch({
      # 方式2：尝试BlueprintEncodeData
      ref <- BlueprintEncodeData()
      message("成功加载 BlueprintEncodeData")
      return(ref)
    }, error = function(e) {
      message("方式2失败: ", e$message)
      
      tryCatch({
        # 方式3：尝试DatabaseImmuneCellExpressionData
        ref <- DatabaseImmuneCellExpressionData()
        message("成功加载 DatabaseImmuneCellExpressionData")
        return(ref)
      }, error = function(e) {
        message("方式3失败: ", e$message)
        
        # 方式4：尝试MonacoImmuneData
        tryCatch({
          ref <- MonacoImmuneData()
          message("成功加载 MonacoImmuneData")
          return(ref)
        }, error = function(e) {
          message("方式4失败: ", e$message)
          stop("无法加载任何参考数据集，请检查网络连接和包安装")
        })
      })
    })
  })
}

# 尝试加载参考数据集
hpca.se <- get_reference()

##5.4 运行SingleR（添加并行处理）
library(BiocParallel)
pred <- SingleR(
  test = pbmc_flt9_matrix,
  ref = hpca.se,
  labels = hpca.se$label.main,
  BPPARAM = MulticoreParam(workers = 4)  # 使用4核心加速
)

##5.5 添加注释结果
pbmc_flt9$celltype <- pred$labels

##5.6 可视化
#library(patchwork)
#p1 <- DimPlot(pbmc_flt9, reduction = "tsne", group.by = "celltype", label = TRUE) + 
#  NoLegend() + ggtitle("t-SNE")
#p2 <- DimPlot(pbmc_flt9, reduction = "umap", group.by = "celltype", label = TRUE) + 
#  NoLegend() + ggtitle("UMAP")

# 设置绘图参数增加可读性
p1 <- DimPlot(pbmc_flt9, reduction = "tsne", group.by = "celltype", label = TRUE,
              pt.size = 0.8,        # 增大点大小
              raster = TRUE,        # 启用栅格化
              raster.dpi = c(1024, 1024),  # 提高栅格化分辨率
              shuffle = TRUE,       # 随机打乱点顺序
              label.size = 3        # 调整标签大小
) + 
  NoLegend() + 
  ggtitle("t-SNE") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(pbmc_flt9, reduction = "umap", group.by = "celltype", label = TRUE,
              pt.size = 0.8,        # 增大点大小
              raster = TRUE,        # 启用栅格化
              raster.dpi = c(1024, 1024),  # 提高栅格化分辨率
              shuffle = TRUE,       # 随机打乱点顺序
              label.size = 2        # 调整标签大小
) + 
  NoLegend() + 
  ggtitle("UMAP") +
  theme(plot.title = element_text(hjust = 0.5))


png("P5.SingleR2.png",width = 5700,height = 2600,res=300)
p2+ p1
dev.off()

write.csv(pred,"/home/LSM/LSM_scRNAseq/20241118BFlin_UC/Ana_res/singleR_anno.csv")



#6.Target_cells_analysis
##6.1.提取成纤维细胞子集
fibroblast_subset <- subset(pbmc_flt9, subset = celltype == "Fibroblasts")

##6.2.重新计算高变基因（仅针对成纤维细胞）
fibroblast_subset <- FindVariableFeatures(
  fibroblast_subset,
  selection.method = "vst",   # 使用方差稳定变换方法
  nfeatures = 2000            # 选择前2000个高变基因
)

## 6.3.提取并查看高变基因
var_genes <- VariableFeatures(fibroblast_subset)
head(var_genes, 20)  # 查看前20个高变基因
write.csv(var_genes,"/home/LSM/LSM_scRNAseq/20241118BFlin_UC/Ana_res/var_genes_fib.csv")

## 6.4. 可视化高变基因
###6.4.1 高变基因识别图
png("/home/LSM/LSM_scRNAseq/20241118BFlin_UC/Ana_res/P6.4.1VariableFeaturePlot_fib.png",width = 5700,height = 2600,res=300)
var_plot <- VariableFeaturePlot(fibroblast_subset)
LabelPoints(
  plot = var_plot,
  points = head(var_genes, 10),  # 标记前10个基因
  repel = TRUE
)
dev.off()

###6.4.2 高变基因热图
#### 先进行缩放处理
fibroblast_subset <- ScaleData(fibroblast_subset, features = var_genes)
#### 绘制热图
png("/home/LSM/LSM_scRNAseq/20241118BFlin_UC/Ana_res/P6.4.2Var_gene_Plot_fib_heatmap.png",width = 5700,height = 2600,res=300)
DoHeatmap(
  fibroblast_subset,
  features = head(var_genes, 30),  # 显示前30个高变基因
  group.by = "orig.ident",         # 按样本分组
  size = 3                         # 基因名字体大小
) + 
  theme(axis.text.y = element_text(size = 8))  # 调整Y轴字体大小
dev.off()


#7. 功能富集分析（了解高变基因的生物学意义）
library(clusterProfiler)
library(org.Hs.eg.db)

##7.1 将基因名转换为ENTREZ ID
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = var_genes,
  keytype = "SYMBOL",
  column = "ENTREZID"
)

##7.2去除NA值
entrez_ids <- na.omit(entrez_ids)

##7.3GO富集分析
go_enrich <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # 生物过程
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

##7.4可视化GO结果
png("/home/LSM/LSM_scRNAseq/20241118BFlin_UC/Ana_res/P7.4var_gene_fib_GO.png",width = 5700,height = 2600,res=300)
dotplot(go_enrich, showCategory = 15) + 
  ggtitle("GO Enrichment of Fibroblast Variable Genes")
dev.off()

##7.5. 在UMAP上可视化特定高变基因的表达
# 选择前5个高变基因
top_genes <- head(var_genes, 5)

#7.6 绘制UMAP图
png("/home/LSM/LSM_scRNAseq/20241118BFlin_UC/Ana_res/P7.6Var_gene_fib_UMAP.png",width = 5700,height = 2600,res=300)
FeaturePlot(
  fibroblast_subset,
  features = top_genes,
  reduction = "umap",
  ncol = 3,
  order = TRUE,  # 高表达细胞显示在上层
  pt.size = 0.5
) & 
  scale_color_gradientn(colors = c("blue", "white", "red"))  # 自定义颜色
dev.off()

# 7. 保存结果
write.csv(
  data.frame(Gene = var_genes),
  file = "/home/LSM/LSM_scRNAseq/20241118BFlin_UC/Ana_res/fibroblast_variable_genes.csv",
  row.names = FALSE
)



#6.Target_cells_analysis
##6.1-1.提取巨噬细胞子集
mcp_subset <- subset(pbmc_flt9, subset = celltype == "Macrophage")

##6.2-1.重新计算高变基因（仅针对巨噬细胞）
mcp_subset <- FindVariableFeatures(
  mcp_subset,
  selection.method = "vst",   # 使用方差稳定变换方法
  nfeatures = 2000            # 选择前2000个高变基因
)

## 6.3.提取并查看高变基因
var_genes_mcp <- VariableFeatures(mcp_subset)
head(var_genes_mcp, 20)  # 查看前20个高变基因
write.csv(var_genes_mcp,"/home/LSM/LSM_scRNAseq/20241118BFlin_UC/Ana_res/var_genes_mcp.csv")

## 6.4. 可视化高变基因
###6.4.1 高变基因识别图
png("/home/LSM/LSM_scRNAseq/20241118BFlin_UC/Ana_res/P6.4.1-1VariableFeaturePlot_mcp.png",width = 3700,height = 1600,res=300)
var_plot <- VariableFeaturePlot(mcp_subset)
LabelPoints(
  plot = var_plot,
  points = head(var_genes, 10),  # 标记前10个基因
  repel = TRUE
)
dev.off()

###6.4.2 高变基因热图
#### 先进行缩放处理
mcp_subset <- ScaleData(mcp_subset, features = var_genes)
#### 绘制热图
png("/home/LSM/LSM_scRNAseq/20241118BFlin_UC/Ana_res/P6.4.2-1Var_gene_Plot_mcp_heatmap.png",width = 5700,height = 2600,res=300)
DoHeatmap(
  mcp_subset,
  features = head(var_genes_mcp, 30),  # 显示前30个高变基因
  group.by = "orig.ident",         # 按样本分组
  size = 3                         # 基因名字体大小
) + 
  theme(axis.text.y = element_text(size = 8))  # 调整Y轴字体大小
dev.off()


#7. 功能富集分析（了解高变基因的生物学意义）
library(clusterProfiler)
library(org.Hs.eg.db)

##7.1 将基因名转换为ENTREZ ID
entrez_ids_mcp <- mapIds(
  org.Hs.eg.db,
  keys = var_genes_mcp,
  keytype = "SYMBOL",
  column = "ENTREZID"
)

##7.2去除NA值
entrez_ids_mcp <- na.omit(entrez_ids)

##7.3GO富集分析
go_enrich_mcp <- enrichGO(
  gene = entrez_ids_mcp,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # 生物过程
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

##7.4可视化GO结果
png("/home/LSM/LSM_scRNAseq/20241118BFlin_UC/Ana_res/P7.4-1var_gene_mcp_GO.png",width = 5700,height = 2600,res=300)
dotplot(go_enrich_mcp, showCategory = 15) + 
  ggtitle("GO Enrichment of Macrophage Variable Genes")
dev.off()

##7.5. 在UMAP上可视化特定高变基因的表达
# 选择前5个高变基因
top_genes_mcp <- head(var_genes_mcp, 5)

#7.6 绘制UMAP图
png("/home/LSM/LSM_scRNAseq/20241118BFlin_UC/Ana_res/P7.6-1Var_gene_mcp_UMAP.png",width = 5700,height = 2600,res=300)
FeaturePlot(
  mcp_subset,
  features = top_genes_mcp,
  reduction = "umap",
  ncol = 3,
  order = TRUE,  # 高表达细胞显示在上层
  pt.size = 0.5
) & 
  scale_color_gradientn(colors = c("blue", "white", "red"))  # 自定义颜色
dev.off()

# 7. 保存结果
write.csv(
  data.frame(Gene = var_genes_mcp),
  file = "/home/LSM/LSM_scRNAseq/20241118BFlin_UC/Ana_res/mcp_variable_genes.csv",
  row.names = FALSE
)

save.image("pbmc_flt9.RData")



############################
##auto_SG_scRNAseq_report
setwd(DIR00)
projNum <- c("SG20251115EXP01")
rmd_content <- paste0("
---
title: Guangzhou Sengene Technology Co,ltd;.
output:
  html_document:
    toc: true
    toc_float: true
    toc_depth: 5 
---

![](https://upload-images.jianshu.io/upload_images/15945431-d8603d7f4cf68dda.jpg?imageMogr2/auto-orient/strip%7CimageView2/2/w/620)

<center><font size='10'>生物信息学分析结题报告</center></font size>

<div style='text-align: right;'><font size='4'>项目编号:",projNum,"</font></div>


## 1.scRNA分析简介</br>
### 1.1单细胞测序基本原理

### 1.2WGCNA分析操作步骤</br>
具体地说，参照上面的原理，本次分析依次进行了以下操作：</br>
* 数据输入、清洗和预处理
  + 载入数据
  + 检查缺失值和识别离群值（异常值）
  + 载入表型数据
</br>
* 构建表达网络
  + 参数设置
  + 一步法构建网络和模块检测
</br>
* 模块与表型数据关联并识别重要基因
  + 模块-表型数据关联
  + 基因与表型数据的关系(GS)和重要模块(MM):基因显著性和模块成员
  + 模块内分析：鉴定具有高GS和高MM的基因
  + 输出网络分析结果
</br>
* 网络交互分析（GO注释等）</br>

## 2.数据集信息</br>
本次分析使用数据集为GSE180857，分组信息及样本信息见GEO数据库。</br></br>

## 3.材料与方法</br>

WGCNA(weighted gene co-expression network analysis，权重基因共表达网络分析)
使用R4.3.1进行WGCNA分析。WGCNA进行共表达网络分析。构建'unsigned'网络，'pearsoon'法计算排序后的特征基因模块（MEs_col）与性状数据(traitData)间的相关性。在兴趣模块基因中筛选用于后续分析的相关基因。</br>

## 4.文章引用与致谢</br>
如果您的研究课题使用了三君科技的测序和分析服务，您在论文发表时，在Method部分或Acknowledgements部分引用或提及三君科技，是对我们工作的认可，我们表示衷心感谢。以下语句可供参考：

>方法部分：</br>
The cDNA libraries were sequenced on the Illumina sequencing platform by Guangzhou Sengene Technology Co., Ltd (Guangzhou, China).
</br>

>Acknowledgements部分：</br>
We are grateful to thank Guangzhou Sengene Technology Co., Ltd for assisting in sequencing and/or bioinformatics analysis.

</br>

## 5.分析结果</br>
使用测序后Count值数据作为输入文件，原始数据中共有",
                      dim(rawdat)[1],"个基因。原始进行数据清洗（注释，去重，过滤低表达基因（保留至少80%样本中Count>10））后，剩余基因数为",
                      nrow(count_matrix),"，随后进行log2_VST转置，作为WGCNA输入矩阵。WGCNA分析具体参数如下：flodChange:",
                      foldChange,",padj:",padj,",筛选中位绝对偏差前65%的基因,基因数量为:",dim(dataExprVar1)[1],
                      "构建无标度网络,软阈值为",power,",共分为",ncol(net$MEs),"个模块,根据模块与表型相关分析，选取表型为",pheno,"的模块ME",module,
                      "作为后续研究内容。

###5.1样本聚类分析</br>
```{r fig1, echo=FALSE}
knitr::include_graphics(file.path(DIR05,'S1sampletree.png'))
```

>Fig.1 Sample clustering and outliersc detection .</br>
图一.样本聚类图。通过对",dim(exprDat)[2],"个样本进行聚类分析，未见异常离群样品，继续进行后续分析。</br>


###4.2软阈值筛选</br>
```{r fig2, echo=FALSE}
knitr::include_graphics(file.path(DIR05,'S2Soft threshold.png'))
```

>Fig.2 Soft Threshold screening.</br>
图二.软阈值筛选。本次分析通过绘制样品聚类查看分组信息和有无异常样品后，获得经验性软阈值为",power,"，符合无向网络取值，可继续后续分析。</br>

###4.3基因筛选及模块构建</br>
```{r fig3, echo=FALSE}
knitr::include_graphics(file.path(DIR05,'S3Module colors.png'))
```

>Fig.3 Cluster Dendrogram and module construction.</br>
图三.基因筛选及模块构建。将表达谱中的所有基因进行筛选，选取中位绝对偏差前65%的基因，得到",dim(dataExprVar1)[1],"个差异表达基因，共构建52个模块。其中grey灰色为未分入任何模块的基因。


###4.4层级聚类树及特征基因邻接分析</br>
```{r fig4, echo=FALSE}
knitr::include_graphics(file.path(DIR05,'S4Eigengene adjacency heatmap.png'))
```

>Fig.4 Eigengene adjacency heatmap</br>
图4.层级聚类树及特征基因邻接热图。根据基因间表达量进行聚类所得到的各模块间的相关性图，构建共表达网络，划分模块，合并相似模块。

###4.5模块与分组表型数据关联分析</br>
```{r fig5, echo=FALSE}
knitr::include_graphics(file.path(DIR05,'/S5Module-trait relationships.png'))
```

>Fig.5 Modules and phenotype association</br>
图五.模块与分组表型数据关联分析。根据数据集分类，将带“NC”标签的样品归为对照组（NC），将带“polysaccharide1”标签的样品归为干预组（polysaccharide1）将带“polysaccharide2”标签的样品归为干预组（polysaccharide2），进行表型与模块关联分析。将相关系数（negative correlation coefficient）为0.88, p = 0.002的ME",module,"模块，选为兴趣模块，重点关注",pheno,"组。


###4.6兴趣模块基因分布及基因显著性</br>
```{r fig6, echo=FALSE}
knitr::include_graphics(file.path(DIR05,'S7GSMM.png'))
```

>Fig.6 Module membership and gene significance.</br>
图6.兴趣模块基因分布及基因显著性。提取MEturquoise模块中polysaccharide2组基因，设制条件为GS > 0.75 & MM > 0.75筛选出与性状高度相关的基因，也是与性状相关的模型的关键基因进行分布展示，共有222个基因与“ polysaccharide2”表型显著相关，即与NC相比，222个基因变化与MCAO呈正相关。

#```{r fig7, echo=FALSE}
#knitr::include_graphics(file.path(DIR06,'/ALL/ALLstring_hires_image.png'))
#knitr::include_graphics(file.path(DIR06,'/MCL1/MCL1.png'))
#knitr::include_graphics(file.path(DIR06,'/MCL2/MCL2.png'))
#knitr::include_graphics(file.path(DIR06,'/MCL3/MCL3.png'))
#```

>Fig.7 STRING通路富集分析。将222个基因构建PPI网络，使用MCL算法聚类，最终聚为42个cluster，排名前三的Cluster如上所示。

#```{r fig8, echo=FALSE}
#knitr::include_graphics(file.path(DIR06,'/MCL2/enrichment_Process_sim0.8_graph.png')
#```

#>Fig.8 STRING通路富集分析。Cluster2富集于GO:0035694	Mitochondrial protein catabolic process，GO:0010506	Regulation of autophagy，
#GO:0097345	Mitochondrial outer membrane permeabilization，
#GO:0010917	Negative regulation of mitochondrial membrane potential，
#GO:0071456	Cellular response to hypoxia，其中Dram1,Bnip3l,Bmf,Bnip3参与GO:0010506	Regulation of autophagy，可作为后续研究通路。


[^WGCNAP1]:Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics. 2008;9:559. Published 2008 Dec 29. doi:10.1186/1471-2105-9-559IF:
        
        
        
         2.9 Q1

")

setwd(DIR00)
file_path <- paste0(projNum,"_auto_WGCNA_report.Rmd")
writeLines(rmd_content, con = file_path)
library(rmarkdown)
render(paste0(projNum,"_auto_WGCNA_report.Rmd"), output_format = "html_document")
