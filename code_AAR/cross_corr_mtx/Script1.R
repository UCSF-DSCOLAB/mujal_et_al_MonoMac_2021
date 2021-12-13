library(Seurat)
get_genes <- function(genes, sobj, drop=FALSE) {
  ###
  # DESCRIPTION
  # Get a dataframe of gene values per cell from the input seurat object
  #
  # INPUTS
  # genes: A character vector containing gene names
  # sobj: The Seurat object
  # OUTPUT
  # A dataframe of input genes as rows and all cells as columns
  ###
  gene_idx <- sapply(genes, function(x) { match(x, rownames(sobj)) })
  if (sum(is.na(gene_idx)) > 0) {
    fail_genes <- names(gene_idx[is.na(gene_idx)])
    print("The following genes are not in the gene list.")
    print(fail_genes)
    if (drop){
      print("Dropping genes and continuing.")
      genes <- names(gene_idx[!is.na(gene_idx)])
    } else {
      return(1)  
    }
  }
  genes <- as.data.frame(as.matrix(sobj@assays$RNA@data[genes,]))
}

saturate <- function(vec, sat=0, binary=FALSE){
  ###
  # DESCRIPTION
  # A Function to convert a vector of scores into a saturated vectore of scores. A saturated vector is one where all values below the 
  # provided "saturation" (percentile of data) are set to 0. If the binary flag is specified, all values greater than or equal to the
  # saturation will be set to 1.
  #
  # INPUTS
  # vec: A numeric vector of scores
  # sat: A value to saturate the vector at (float (0.0-1.0) or percent (1.0-100.0))
  # binary: A flag to indicate if we should make the output vector a binary one.
  #
  # OUTPUT
  # A vector of saturated scores
  ###
  sat = if (sat > 1.0) sat/100 else sat
  z <- quantile(vec, sat)
  for (i in 1:length(vec)){
    if (vec[i] < z) {
      vec[i] = 0
    } else if(binary) {
      vec[i] = 1
    }
  }
  vec
}


setwd('~/projects/adriana_paper/cross_corr_mtx/')

load('P1.Robj')
m1m2 <- read.table('M1M2_signatures.tsv', header = T, colClasses = c('character', 'factor'))
P1 <- UpdateSeuratObject(P1)
Idents(P1) <- P1@meta.data[["res.0.6"]]

data <- list(
  "all" = P1,
  "no_dcs" = subset(P1, idents=c("0", "1", "2", "3", "5"))
)

for (frac in names(data)){
  # Process M1
  m1_genes <- m1m2[m1m2$cluster=='M1',]$gene
  gene_df <- get_genes(m1_genes, data[[frac]], drop=TRUE)
  score <- colMeans(gene_df)
  data[[frac]]@meta.data$m1_score_raw <- score
  data[[frac]]@meta.data$m1_score_50 <- saturate(vec=score, sat=0.5, binary=FALSE)    # Saturate at 50%
  data[[frac]]@meta.data$m1_score_50b <- saturate(vec=score, sat=0.5, binary=TRUE)    # Saturate at 50% binary
  
  # Process M2
  m2_genes <- m1m2[m1m2$cluster=='M2',]$gene
  gene_df <- get_genes(m2_genes, data[[frac]], drop=TRUE)
  score <- colMeans(gene_df)
  data[[frac]]@meta.data$m2_score_raw <- score
  data[[frac]]@meta.data$m2_score_50 <- saturate(vec=score, sat=0.5, binary=FALSE)    # Saturate at 50%
  data[[frac]]@meta.data$m2_score_50b <- saturate(vec=score, sat=0.5, binary=TRUE)    # Saturate at 50% binary
  
  # Process Both
  m1m2_genes <- m1m2$gene
  gene_df <- get_genes(m1m2_genes, data[[frac]], drop=TRUE)
  score <- colMeans(gene_df)
  data[[frac]]@meta.data$m1m2_score_raw <- score
  data[[frac]]@meta.data$m1m2_score_50 <- saturate(vec=score, sat=0.5, binary=FALSE)    # Saturate at 50%
  data[[frac]]@meta.data$m1m2_score_50b <- saturate(vec=score, sat=0.5, binary=TRUE)    # Saturate at 50% binary
  
  # Make a featureplot 
  
  pdf(paste0("clusters_", frac, ".pdf"), paper='a4r')
  p <- TSNEPlot(object = data[[frac]], pt.size=5)
  print(p)
  dev.off()
  pdf(paste0("featureplot_", frac, ".pdf"), paper='a4r')
  p <- FeaturePlot(object = data[[frac]],
              features = c("m1_score_raw", "m1_score_50", "m1_score_50b",
                                "m2_score_raw", "m2_score_50", "m2_score_50b",
                                "m1m2_score_raw", "m1m2_score_50", "m1m2_score_50b"),
              cols=c("light grey", "red"), pt.size = 1.5)
  print(p)
  dev.off()
}

pdf(paste0("m1_score_50_nodcs.pdf"), paper='a4r')
p <- FeaturePlot(object = data[["no_dcs"]],
                 features = c("m1_score_50"),
                 cols=c("light grey", "red"), pt.size = 1.5)
print(p)
dev.off()

pdf(paste0("m2_score_50_nodcs.pdf"), paper='a4r')
p <- FeaturePlot(object = data[["no_dcs"]],
                 features = c("m2_score_50"),
                 cols=c("light grey", "red"), pt.size = 1.5)
print(p)
dev.off()
