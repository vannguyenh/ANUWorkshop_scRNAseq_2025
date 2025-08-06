####################################
# 0. Load the libraries
####################################

library("Seurat")
library("ggplot2")
library("cowplot")
library("dplyr")
library("Matrix")
library("viridis")
library("gprofiler2")

####################################
# 1. Load one sample
####################################
# read the matrix file Rep1?ICBdT_data
## the matrix of the raw read counts for all cels in a sample
Rep1_ICBdT_data <- Read10X_h5("data/Rep1_ICBdT-sample_filtered_feature_bc_matrix.h5")
Rep1_ICB_data <- Read10X_h5("data/Rep1_ICB-sample_filtered_feature_bc_matrix.h5")

####################################
# 1. Create the Seurat Object
####################################

# counts: unnormalised data such as raw -- raw data with raw counts
# project is usually used as the sample name
# min.cells: features detected at least in this "number" of cells
# min.features: include cells where at least "this many" features are detected
Rep1_ICBdT_data_seurat_obj <- CreateSeuratObject(counts = Rep1_ICBdT_data,
                                                 project = "Rep1_ICBdT",
                                                 min.cells = 10,
                                                 min.features = 100)

Rep1_ICB_data_seurat_obj <- CreateSeuratObject(counts = Rep1_ICB_data,
                                               project = "Rep1_ICB",
                                               min.cells = 10,
                                               min.features = 100)

####################################
# 1.1 Understand the Seurat Object
####################################
# metadata table with 
# each row representing a cell
# each column representing an attribute of the cells
## org.ident: sample this cell come from
## nCount_RNA: # of reads in that cell
## nFeature_RNA: # of genes in that cell
head(Rep1_ICBdT_data_seurat_obj@meta.data)
# or we can access the medata like this
head(Rep1_ICBdT_data_seurat_obj[[]])

# view a column of the metadata
head(Rep1_ICBdT_data_seurat_obj@meta.data$orig.ident)
# access each column in the data
head(Rep1_ICBdT_data_seurat_obj[["orig.ident"]])

# add a column to metadata table
## name the new column as "percent.mt" - percentage of mitchondrial genes
## "percent.mt" is useful for removing low-quality cells later
# PercentageFeatureSet() -> calculate the percentage of all counts belonging to a given set of features
## pattern="^mt-": match features starting with 'mt-'
## assay="RNA" : calculate this on the assay RNA
Rep1_ICBdT_data_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(Rep1_ICBdT_data_seurat_obj,
                                                                   pattern = "^mt-",
                                                                   assay = "RNA")
Rep1_ICB_data_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(Rep1_ICB_data_seurat_obj,
                                                                 pattern = "^mt-",
                                                                 assay = "RNA")
####################################
# 1.2 Quick look at counts, genes (nFeature_RNA) and mt_genes (percent.mt)
####################################

# violin plot for the # of reads, # of genes and % of mit. genes
## 2 peaks: one is close to 0, another is around 6000
## meaning of pt.size: size of points (values)
## layer: layer ot pull expression data from e.g "counts" or "data"
p1 <- VlnPlot(Rep1_ICB_data_seurat_obj,
              layer = "counts",
              features = "nCount_RNA",
              pt.size = 0) +
  scale_y_continuous(breaks=c(0,10000, 20000, 30000)) # add y-axis ticks
p1

p2 <- VlnPlot(Rep1_ICB_data_seurat_obj,
              layer = "counts",
              features = c("nFeature_RNA"),
              pt.size = 0) +
  scale_y_continuous(breaks=c(0,1000,2000,3000,4000,5000,6000))
p2

p3 <- VlnPlot(Rep1_ICB_data_seurat_obj,
              layer = "counts",
              features = "percent.mt",
              pt.size = 0) +
  scale_y_continuous(breaks=c(0,5,10,15,20,25))

p3

# combine 3 plots
p <- plot_grid(p1, p2, p3, ncol=3)
p

####################################
# 2. all 6 Samples
####################################

####################################
# 2.1 Create Seurat Objects for all 6 Samples
####################################

sample_names <- c("Rep1_ICBdT", "Rep3_ICBdT", "Rep5_ICBdT",
                  "Rep1_ICB", "Rep3_ICB", "Rep5_ICB")
sample.data <- list()

for (sample in sample_names) {
  cat("\nWorking on sample: ", sample, "\n")
  path <- paste("data/", sample, "-sample_filtered_feature_bc_matrix.h5", sep="")
  cat("Input file path: ", path, "\n")
  data <- Read10X_h5(path)
  seurat_obj <- CreateSeuratObject(counts = data,
                                   project = sample,
                                   min.cells = 10,
                                   min.features = 100)
  sample.data[[sample]] <- seurat_obj
}
# access one time from the list
sample.data$Rep1_ICB
sample.data[["Rep1_ICBdT"]]

####################################
# 2.2 Calculate mt RNA percentage for all 6 Samples
####################################

for (sample in sample_names) {
  sample.data[[sample]][["percent.mt"]] <- PercentageFeatureSet(sample.data[[sample]],
                                                                pattern = "^mt-",
                                                                assay = "RNA")
}

####################################
# 2.3 Quality Assessment for all 6 Samples
####################################

# Violin plots
for (sample in sample_names) {
  jpeg(sprintf("%s/%s_unfilteredQC.jpg", "figures", sample),
       width = 16, height = 5, units = 'in', res = 150)
  p1 <- VlnPlot(sample.data[[sample]], features = "nCount_RNA", pt.size = 0)
  p2 <- VlnPlot(sample.data[[sample]], features = "nFeature_RNA", pt.size = 0) +
    scale_y_continuous(breaks=c(0,300,500,1000,2000,4000))
  p3 <- VlnPlot(sample.data[[sample]], features = "percent.mt", pt.size = 0) +
    scale_y_continuous(breaks=c(0, 12.5, 25, 50))
  p <- plot_grid(p1, p2, p3, ncol = 3)
  print(p)
  dev.off()
}

####################################
# 2.4 Filtering based on mtRNA percentage
####################################
#Rep1_ICBdT_data_seurat_obj[["keep_cell_percent.mt"]] <- ifelse(Rep1_ICBdT_data_seurat_obj[["percent.mt"]] <=12, TRUE, FALSE)
#Rep1_ICB_data_seurat_obj[["keep_nFeature_RNA"]] <- ifelse(Rep1_ICB_data_seurat_obj[["nFeature_RNA"]] >= 1000, TRUE, FALSE)
#Rep1_ICB_data_seurat_obj[["keep_nCount_RNA"]] <- ifelse(Rep1_ICB_data_seurat_obj[["nCount_RNA"]], TRUE, FALSE)

# apply this for all 6 samples
for (sample in sample_names) {
  cat("\nWorking on sample: ", sample, "\n")
  # call out the seurat object: sample.data[[sample]]

  sample.data[[sample]][["keep_cell_percent.mt"]] <- ifelse(sample.data[[sample]][["percent.mt"]] <=12, TRUE, FALSE)
  sample.data[[sample]][["keep_nFeature"]] <- ifelse(sample.data[[sample]][["nFeature_RNA"]] >= 1000, TRUE, FALSE)
  sample.data[[sample]][["keep_nCount_RNA"]] <- ifelse(sample.data[[sample]][["percent.mt"]] >= 2000, TRUE, FALSE)
}


}
# number of cells after filtering
for (sample in sample_names) {
  print(sample)
  print("Before filtering: ")
  print(sum(sample.data[[sample]][["nFeature_RNA"]] & 
              sample.data[[sample]][["percent.mt"]]))
  print("After filtering: ")
  print(sum(sample.data[[sample]][["keep_nFeature"]] & 
              sample.data[[sample]][["keep_cell_percent.mt"]]))
}
####################################
# 2.5 Merging samples & Filtering merged data
####################################

for (sample in sample_names) {
  cat("\nInformation of sample: ", sample, "\n")
  head(sample.data[[sample]][[]])
}
# merge samples
unfiltered_merged <- merge(x=sample.data$Rep1_ICBdT,
                           y=c(sample.data$Rep3_ICBdT, sample.data$Rep5_ICBdT, 
                               sample.data[['Rep1_ICB']], sample.data[['Rep3_ICB']], sample.data[['Rep5_ICB']]),
                           add.cell.ids=sample_names) #add sample names to the end of cell barcode
# filter merged data
# JoinLayers() -> merge them into one layer which is a single matrix for further analysis
merged <- subset(unfiltered_merged, nFeature_RNA > 1000 & percent.mt < 12) |> JoinLayers()


####################################
# 3. Normalise and Scale the Data
####################################

# NormalizeData -> gene expression levels comparable between cells
# --> remove technical differences between cells, primarily sequencing depth
# -- perform a "LogNormalize" transformation: calculate relative counts, scale the data, log-transform data
merged <- NormalizeData(merged, assay = "RNA", normalization.method = "LogNormalize",
                        scale.factor = 10000)

# ScaleData -> gene expression levels comparable between genes

