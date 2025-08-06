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
# 1. Create the Seurat Object
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

####################################
# 2. Calculate mt RNA percentage for all 6 Samples
####################################

for (sample in sample_names) {
  sample.data[[sample]][["percent.mt"]] <- PercentageFeatureSet(sample.data[[sample]],
                                                                pattern = "^mt-",
                                                                assay = "RNA")
}

####################################
# 3. Quality Assessment for all 6 Samples
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
# 4. Filtering based on mtRNA percentage
####################################

# apply this for all 6 samples
for (sample in sample_names) {
  cat("\nWorking on sample: ", sample, "\n")
  # call out the seurat object: sample.data[[sample]]
  
  sample.data[[sample]][["keep_cell_percent.mt"]] <- ifelse(sample.data[[sample]][["percent.mt"]] <=12, TRUE, FALSE)
  sample.data[[sample]][["keep_nFeature"]] <- ifelse(sample.data[[sample]][["nFeature_RNA"]] >= 1000, TRUE, FALSE)
  sample.data[[sample]][["keep_nCount_RNA"]] <- ifelse(sample.data[[sample]][["percent.mt"]] >= 2000, TRUE, FALSE)
}

# number of cells before and after filtering
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
# 5. Merging samples & Filtering merged data
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

