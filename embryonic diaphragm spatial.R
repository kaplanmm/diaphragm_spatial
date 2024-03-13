## load packages
library(Seurat)
library(ggplot2)
library(tibble)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(ggcorrplot)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(writexl)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(ggraph)
library(clustree)
library(readxl)

set.seed(12345)


#### E14 Controls
control_e14_1 <- Seurat::Load10X_Spatial(
  # The directory contains the read count matrix H5 file and the image data in a subdirectory called `spatial`. 
  data.dir = "/Users/mehmetmahsumkaplan/Visium_data/1st Round/Visium_Mehmet_2/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  slice = "batch1", # specify name of the stored image
  filter.matrix = TRUE, 
  to.upper = FALSE
)
control_e14_2 <- Seurat::Load10X_Spatial(
  # The directory contains the read count matrix H5 file and the image data in a subdirectory called `spatial`. 
  data.dir = "/Users/mehmetmahsumkaplan/Visium_data/1st Round/Visium_Mehmet_4/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  slice = "batch2", # specify name of the stored image
  filter.matrix = TRUE, 
  to.upper = FALSE
)
control_e14_3 <- Seurat::Load10X_Spatial(
  # The directory contains the read count matrix H5 file and the image data in a subdirectory called `spatial`. 
  data.dir = "/Users/mehmetmahsumkaplan/Visium_data/2nd and 3rd Round/MK8_D1_result/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  slice = "batch3", # specify name of the stored image
  filter.matrix = TRUE, 
  to.upper = FALSE
)





## E18 Controls
control_e18_1 <- Seurat::Load10X_Spatial(
  # The directory contains the read count matrix H5 file and the image data in a subdirectory called `spatial`. 
  data.dir = "/Users/mehmetmahsumkaplan/Visium_data/2nd and 3rd Round/MK11/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  slice = "batch4", # specify name of the stored image
  filter.matrix = TRUE, 
  to.upper = FALSE
)
control_e18_2 <- Seurat::Load10X_Spatial(
  # The directory contains the read count matrix H5 file and the image data in a subdirectory called `spatial`. 
  data.dir = "/Users/mehmetmahsumkaplan/Visium_data/2nd and 3rd Round/MK9/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  slice = "batch5", # specify name of the stored image
  filter.matrix = TRUE, 
  to.upper = FALSE
)



## E14 Dysgenic
dysgenic_e14_1 <- Seurat::Load10X_Spatial(
  # The directory contains the read count matrix H5 file and the image data in a subdirectory called `spatial`. 
  data.dir = "/Users/mehmetmahsumkaplan/Visium_data/1st Round/Visium_Mehmet_1/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  slice = "batch6", # specify name of the stored image
  filter.matrix = TRUE, 
  to.upper = FALSE
)

## E18 Dysgenic
dysgenic_e18_1 <- Seurat::Load10X_Spatial(
  # The directory contains the read count matrix H5 file and the image data in a subdirectory called `spatial`. 
  data.dir = "/Users/mehmetmahsumkaplan/Visium_data/2nd and 3rd Round/MK12/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  slice = "batch7", # specify name of the stored image
  filter.matrix = TRUE, 
  to.upper = FALSE
)




## E18 bCat KO
bcat_e18_1 <- Seurat::Load10X_Spatial(
  # The directory contains the read count matrix H5 file and the image data in a subdirectory called `spatial`. 
  data.dir = "/Users/mehmetmahsumkaplan/Visium_data/2nd and 3rd Round/MK10/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  slice = "batch10", # specify name of the stored image
  filter.matrix = TRUE, 
  to.upper = FALSE
)

# add batch
control_e14_1$batch <- "batch1"
control_e14_2$batch <- "batch2"
control_e14_3$batch <- "batch3"
control_e18_1$batch <- "batch4"
control_e18_2$batch <- "batch5"


dysgenic_e14_1$batch <- "batch6"
dysgenic_e18_1$batch <- "batch7"

bcat_e18_1$batch <- "batch10"

#add age
control_e14_1$age <- "E14"
control_e14_2$age <- "E14"
control_e14_3$age <- "E14"
control_e18_1$age <- "E18"
control_e18_2$age <- "E18"
dysgenic_e14_1$age <- "E14"
dysgenic_e18_1$age <- "E18"
bcat_e18_1$age <- "E18"

#add condition
control_e14_1$condition <- "control"
control_e14_2$condition <- "control"
control_e14_3$condition <- "control"
control_e18_1$condition <- "control"
control_e18_2$condition <- "control"
dysgenic_e14_1$condition <- "dysgenic"
dysgenic_e18_1$condition <- "dysgenic"
bcat_e18_1$condition <- "bCat_KO"

# make a list of objects to be integrated,  either of combination was used
experiment.list <- c(control_e14_1, 
                     control_e14_2, 
                     control_e14_3) # for E14.5 samples integrated

experiment.list <- c(control_e18_1,  # for E18.5 samples integrated
                     control_e18_2)

experiment.list <- c(control_e14_1, 
                     control_e14_2, 
                     control_e14_3,
                     control_e18_1, 
                     control_e18_2) # for E14.5 and E18.5 control integrated

experiment.list <- c(control_e14_1, 
                     control_e14_2, 
                     control_e14_3,
                     dysgenic_e14_1) # for E14.5 control vs dysgenic integrated

experiment.list <- c(control_e18_1, 
                     control_e18_2,
                     dysgenic_e18_1) # for E18.5 control vs dysgenic integrated


experiment.list <- c(control_e18_1, 
                     control_e18_2,
                     bcat_e18_1)  # for E18.5 control vs bCat cKO integrated



# run SCT normalization for each sample
for (i in 1:length(experiment.list)) {
  experiment.list[[i]] <- SCTransform(experiment.list[[i]],
                                      assay = "Spatial",
                                      verbose = FALSE)
}



# select the features for downstream integration 
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
experiment <- FindClusters(experiment, 
                           resolution = 0.4) 


###########################################################################
###########################################################################
###########################################################################
DefaultAssay(experiment.integrated) <- "integrated"
DefaultAssay(experiment) <- "integrated"

# SPATIALDIMPLOTS
SpatialDimPlot(experiment, label = F, 
               ncol = 3,
               label.size = 3, 
               pt.size.factor = 2, 
               crop = T, 
               cols = my_cols)+ NoLegend() # cols argument can be customized

### HUGHLIGHT CLUSTERS
plot(SpatialDimPlot(experiment, 
                    cells.highlight = CellsByIdentities(object = experiment,
                                                        idents = "names of the clusters to highlight"), 
                    facet.highlight = T, ncol = 4))

# DIMPLOTS
DimPlot(experiment, 
        reduction = "umap", 
        group.by = "CellType", 
        pt.size = 1, 
        cols = my_cols) # cols argument can be customized for each cluster

# BATCH EFFECTS
FeatureScatter(experiment, 
               "nCount_Spatial", "nFeature_Spatial",
                    pt.size = 1.5, 
               group.by = "batch or age or condition")

DimPlot(experiment, 
        reduction = "umap", 
        group.by = "batch or age or condition")



###########################################################################
###########################################################################
###########################################################################
## SPATIALFEATUREPLOTS
DefaultAssay(experiment) <- "SCT" 
plot(SpatialFeaturePlot(experiment, 
                        features = c( "feature of interest"),
                        alpha = c(0.05, 1),
                        ncol = 3, pt.size.factor = 2) & theme(legend.position = "bottom",
                                                              legend.key.size = unit(0.5, "cm")))



#######################################################################
###########################################################################
###########################################################################
DefaultAssay(experiment) <- "Spatial"
# normalize RNA data and scale them for heatmap and DGE
experiment <- NormalizeData(experiment, verbose = TRUE)
experiment <- ScaleData(experiment, verbose = TRUE)


## FIND ALL MARKERS
Idents(experiment) <- experiment@meta.data$seurat_clusters
all_markers <- FindAllMarkers(experiment)


all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

# HEATMAPS
DoHeatmap(experiment, features = top5$gene, raster=F) +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) +  guides(color=FALSE)


## DEG between clusters/ages/conditions
Idents(experiment) <- experiment@meta.data$seurat_clusters ## or "age" or "condition" instead of seurat_clusters
DEGs <- FindMarkers(experiment, ident.1 = "cluster, condition or age of interest_1", 
                    ident.2 = "cluster, condition or age of interest_2")


## FEATUREPLOTS
FeaturePlot(experiment,features = c("feature of interest"),
            min.cutoff = "q10", max.cutoff = "q90", 
            ncol = 4, split.by = "condition / age") & theme(legend.position = "right") 


### FEATURESCATTERS
FeatureScatter(experiment, "feature_1", "feature_2",
               pt.size = 1.5, group.by = "orig.ident") & NoLegend()


### SUBSETTING only muscle clusters
muscle <- subset(experiment, idents = c("list of clusters identified as muscle"),
                     invert = F)

## VIOLINPLOTS
Idents(experiment) <- experiment@meta.data$condition #or age or seurat cluster
VlnPlot(object = experiment, 
        features = c("feature of interest"), pt.size = 0, cols = my_cols, ncol = 1) & NoLegend() # cols argument can be costumize

## DOTPLOTS in Fig2D
c7_0_9_1 <- subset(experiment, idents = c(0, 7, 9, 1), #subset for clusters to be analyzed (Cluster 7 is NMJ, 9 is MTJ, 0 is default muscle, 1 is central tendon)
                 invert = F)
my_levels <- c(7,0,9,1)
c7_0_9_1@active.ident <- factor(x = c7_0_9_1@active.ident, levels = my_levels)

DefaultAssay(c7_0_9_1) <- "Spatial"
c7_0_9_1 <- NormalizeData(c7_0_9_1, verbose = TRUE)
c7_0_9_1 <- ScaleData(c7_0_9_1, verbose = TRUE)

DotPlot(object = c7_0_9_1, features =  "list of genes") + coord_flip()




## clusterProfiler GENE ONTOLOGY
# run GO 
GO <- enrichGO("list of genes to run GO for",OrgDb = "org.Mm.eg.db",  #list of genes are DEGs between cluster, ages or conditions
               ont = "BP", 
               keyType = "SYMBOL")

# cnetplot
GO_graph <- cnetplot(GO, node_label="all", 
                     color_category='firebrick', 
                     color_gene='steelblue', showCategory = 10)

# treeplot
GO2 <- pairwise_termsim(GO)
treeplot(GO2)



#finding correlation with Chrna1
matrix<-experiment[["Spatial"]]@data
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Chrna1",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
correlations <- as.data.frame(correlations) ## -->> then write in excel















