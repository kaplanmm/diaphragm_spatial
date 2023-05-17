library(Seurat)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(ggcorrplot)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(writexl)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(ggraph)
library(clustree)

set.seed(12345)

## Load 10X data
object <- Seurat::Load10X_Spatial(
  data.dir = "path to outputs/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", 
  slice = "batch1", 
  filter.matrix = TRUE, 
  to.upper = FALSE
)  ## do it for all samples

object$batch <- "batch1" ## add batch
object$age <- "E14 or E18" ## add age
object$condition <- "control, dysgenic or bCAT KO" ## add condition

experiment.list <- c(object1, object2, ....) ## make list of objects


# run normalization for each sample
for (i in 1:length(experiment.list)) {
  experiment.list[[i]] <- SCTransform(experiment.list[[i]],
                                      assay = "Spatial",
                                      verbose = FALSE)
}

# select the feeatures for downstream integration 
experiment.features <- SelectIntegrationFeatures(object.list = experiment.list, nfeatures = 3000)
experiment <- PrepSCTIntegration(object.list = experiment.list, anchor.features = experiment.features, 
                                 verbose = TRUE)

experiment.anchors <- FindIntegrationAnchors(object.list = experiment, normalization.method = "SCT", 
                                             anchor.features = experiment.features, verbose = TRUE, dims = 1:30)
experiment.integrated <- IntegrateData(anchorset = experiment.anchors, normalization.method = "SCT", 
                                       verbose = TRUE, dims = 1:30)

# set the assay to an integrated analysis 
DefaultAssay(experiment.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
experiment <- RunPCA(experiment.integrated, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
experiment <- RunUMAP(experiment, reduction = "pca", dims = 1:25)
experiment <- RunTSNE(experiment, reduction = "pca", dims = 1:25)
experiment <- FindNeighbors(experiment, reduction = "pca", dims = 1:25)

## run Clustree to select resolution
for(i in 1:10){
  experiment <- FindClusters(experiment, resolution = i*0.1, verbose = FALSE)
}

clustree(experiment)

# if 0.4 resolution selected
experiment <- FindClusters(experiment, resolution = 0.4) 

#######################################################################
DefaultAssay(experiment.integrated) <- "integrated"
DefaultAssay(experiment) <- "integrated"

DimPlot(experiment, reduction = "umap", group.by = "seurat_clusters", pt.size = 1, label = T, label.size = 10) & NoLegend()

SpatialDimPlot(experiment, label = F, ncol = 3, label.size = 3, pt.size.factor = 2, crop = T)+ NoLegend()


SpatialDimPlot(experiment, 
               cells.highlight = CellsByIdentities(object = experiment,
                                                   idents = "clusters to highlight"), 
               facet.highlight = T, ncol = 4)


DefaultAssay(experiment) <- "Spatial"
# normalize RNA data and scale them for heatmap and DGE
experiment <- NormalizeData(experiment, verbose = TRUE)
experiment <- ScaleData(experiment, verbose = TRUE)

#### Find Markers
Idents(experiment) <- experiment@meta.data$seurat_clusters
all_markers <- FindAllMarkers(experiment)

Idents(experiment) <- experiment@meta.data$seurat_clusters ## or "age" or "condition" instead of seurat_clusters
DoHeatmap(experiment, features = "genes of interests", raster=F) +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) +  guides(color=FALSE)


## DEG between clusters
DEGs <- FindMarkers(experiment, ident.1 = "cluster, condition or age of interest_1", 
                    ident.2 = "cluster, condition or age of interest_2")


## FeaturePlots
FeaturePlot(experiment,features = c("genes of interest"),
                 min.cutoff = "q10", max.cutoff = "q90")  & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

### FeatureScatter
FeatureScatter(experiment, "feature_1", "feature_2",
                    pt.size = 1.5, group.by = "orig.ident") & NoLegend()
## VlnPlots
VlnPlot(object = experiment,
        features = c("genes of interest"), ncol = 5,
        pt.size = 0, stack = T, flip = T)  & NoLegend()


## SpatialFeaturePlots
DefaultAssay(experiment) <- "SCT"
SpatialFeaturePlot(experiment, 
                        features = c("gene of interest"),
                        alpha = c(0.05, 1),
                        ncol = 3, pt.size.factor = 2) & theme(legend.position = "bottom",
                                                              legend.key.size = unit(0.5, "cm"))



## GENE ONTOLOGY
# run GO 
GO <- enrichGO("list of genes",OrgDb = "org.Mm.eg.db", 
                          ont = "BP", 
                          keyType = "SYMBOL")

# cnetplot
GO_graph <- cnetplot(GO, node_label="all", 
                                color_category='firebrick', 
                                color_gene='steelblue')

# treeplot
GO2 <- pairwise_termsim(GO)
treeplot(GO2)

#finding correlation with Chrna1
matrix<-experiment[["Spatial"]]@data
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Chrna1",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
correlations <- as.data.frame(correlations) ## -->> then write in excel









