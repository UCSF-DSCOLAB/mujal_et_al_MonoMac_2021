#P10
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scMerge)
  library(scater)
  library(DropletUtils)
  library(Seurat)
  library(reshape)
  library(dplyr)
  library(scran)
  library(monocle)
})


#setwd('~/projects/adriana_paper/monocle/mouse/P10')
setwd('/data/shared/krummellab/arrao/projects/adriana_paper/monocle/mouse/P10')

########################################
#  From Adriana
# P10:
#   -Cluster 0 - Ly6c+ H2-Ab+ Mono
#   -Cluster 1 - Ly6c+ Hp+ Mono
#   -Cluster 2 - Ms4a7+ TAM
#   -Cluster 4 - Arg+
#   -Cluster 7 - Ifn-responsive
########################################
# Define all variables

Robj <- 'P10.Robj'
Rvar <- 'P10'
resToUse <- 'res.0.6'

cluster_mapping <- read.table('P10_clusters.tsv', sep="\t", row.names=1, header=TRUE)
sample_name <- 'P10'
cluster_cols <- as.character(cluster_mapping$colors)
names(cluster_cols) <- as.character(cluster_mapping[[sample_name]])
clustersOfInterest <- rownames(cluster_mapping)
clustersOfInterest <- clustersOfInterest[clustersOfInterest != "999"]

mostDistant <- c('1', '2')
expectedRoot <- 'Ly6c+ Hp+ Mono'

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  num_top_to_use <- -1  
} else if (length(args) == 1){
  if (args[1] == 'all'){
    num_top_to_use <- -1  
    } else {
      num_top_to_use <- as.numeric(args[1])
      if (is.na(num_top_to_use) || num_top_to_use < 1){
        stop("num_genes_to_use can only be a number > 0 or 'all'")
      }
    }
} else{
  stop("only accepts one argument at max (num_genes_to_use) that can only be a number > 0 or 'all'")
}
########################################

if (!file.exists('data.Robj')){
  load(Robj)
  # Convert the seurat 2.x object into a Seurat 3.x object
  data = list(seurat_data=UpdateSeuratObject(get(Rvar)))
  rm(list=c(Rvar))
  
  Idents(data$seurat_data) <- data$seurat_data@meta.data[[resToUse]]
  data$seurat_data <- subset(data$seurat_data, idents=clustersOfInterest)
  
  data$seurat_data@meta.data$cluster <- Idents(data$seurat_data)
  data$seurat_data@meta.data$cell_type <- as.character(sapply(data$seurat_data@meta.data$cluster, function(x) {cluster_mapping[as.character(x), sample_name, drop=TRUE]}))
  if (!identical(sort(as.vector(table(data$seurat_data@meta.data$cell_type))), sort(as.vector(table(data$seurat_data@meta.data$cluster))))){
    stop("Something went very wrong")
  }
  Idents(data$seurat_data) <- data$seurat_data@meta.data$cell_type

  data$markers <- FindAllMarkers(data$seurat_data, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.35, test.use='MAST')
  data$top_marker_per_cell <- data$markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)

  data$most_distant <- subset(data$seurat_data, subset= cluster %in% mostDistant)
  data$most_distant_markers <- FindAllMarkers(data$most_distant, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.35, test.use='MAST')
  data$top_most_distant_marker <- data$most_distant_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)

  DimPlot(data$seurat_data, pt.size=1)
  ggsave(file=paste(sample_name, "seurat_tsne.png", sep="_"))

  FeaturePlot(data$seurat_data, features = data$top_marker_per_cell$gene, pt.size=0.5)
  ggsave(file=paste(sample_name, "seurat_top_genes_per_cluster.png", sep="_"))

  DimPlot(data$most_distant, pt.size=1)
  ggsave(file=paste(sample_name, "md_seurat_tsne.png", sep="_"))

  FeaturePlot(data$most_distant, features = data$top_most_distant_marker$gene, pt.size=1)
  ggsave(file=paste(sample_name, "md_seurat_top_genes_per_cluster.png", sep="_"))

  DoHeatmap(data$seurat_data, features=(data$markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))$gene) + scale_color_manual(values = cluster_cols)
  ggsave(file=paste(sample_name, "seurat_top5_heatmap.png", sep="_"))  

  data$raw <- as.matrix(data$seurat_data@assays$RNA@counts)
  # data$raw <- data$raw[,colnames(data$raw)%in%colnames(data$seurat_data@assays$RNA@scale.data)]
  #data$raw <- data$raw[order(rownames(data$raw)),]

#  clusters <- levels(data$top_marker_per_cell$cluster)

#  cluster_ids <- as.character(data$seurat_data@meta.data[colnames(data$raw), "cluster"])
#  cluster_ids <- data.frame(
#    row.names = colnames(data$raw),
#    cluster = cluster_ids,
#    CellType = factor(sapply(cluster_ids, function(x) {cluster_mapping[x, sample_name]}), levels=clusters)
#  )

  cluster_ids <- data.frame(
    row.names = colnames(data$raw),
    cluster =  as.character(data$seurat_data@meta.data$cluster),
    CellType = as.factor(data$seurat_data@meta.data$cell_type)
  )

  phenoData <- new("AnnotatedDataFrame", data = cluster_ids)
  featureData <- new("AnnotatedDataFrame", data = data.frame(row.names=rownames(data$raw), gene_short_name=rownames(data$raw)))

  data$HSMM <- newCellDataSet(as(data$raw, "sparseMatrix"),
                              phenoData = phenoData,
                              featureData = featureData,
                              expressionFamily = negbinomial.size() )

  data$HSMM <- estimateSizeFactors(data$HSMM)
  data$HSMM <- estimateDispersions(data$HSMM)
  data$HSMM <- detectGenes(data$HSMM, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(data$HSMM), num_cells_expressed >= 10))
  data$HSMM <- data$HSMM[expressed_genes,]
  valid_cells <- row.names(subset(pData(data$HSMM), num_genes_expressed >= 500))
  data$HSMM <- data$HSMM[,valid_cells]
  save(list=c('data'), file='data.Robj')
} else {
  load('data.Robj')
}

#data$HSMM.bk <- data$HSMM

if (TRUE) {
  if (num_top_to_use == -1) {
      top_genes <- as.data.frame(data$markers %>% group_by(cluster))
      num_top_to_use <- 'all'
  } else {
    top_genes <- as.data.frame(data$markers %>% group_by(cluster) %>% top_n(n = num_top_to_use, wt = avg_logFC))
  }
  
  suffix <- "genes_trajectories.pdf"
  clusters <- as.character(unique(data$top_marker_per_cell$cluster))
} else {
  if (num_top_to_use == -1) {
      top_genes <- as.data.frame(data$most_distant_markers %>% group_by(cluster))
      num_top_to_use <- 'all'
  } else {
    top_genes <- as.data.frame(data$most_distant_markers %>% group_by(cluster) %>% top_n(n = num_top_to_use, wt = avg_logFC))
  }
  suffix <- "md_genes_trajectories.pdf"
  clusters <- as.character(unique(data$most_distant_markers$cluster))
}
gene_sigs <- list()

for (cluster in clusters){
  gene_sigs[[cluster]] <- top_genes[top_genes$cluster == cluster, 'gene']
}

ordering_genes <- unique(as.character(unlist(gene_sigs, recursive=FALSE)))
data$HSMM <- setOrderingFilter(data$HSMM, ordering_genes)
data$HSMM <- reduceDimension(data$HSMM, max_components = 2,  reduction_method = 'DDRTree', verbose = F)
data$HSMM <- orderCells(data$HSMM)

pdf(paste(sample_name, num_top_to_use, suffix, sep="_"), paper = 'a4r', onefile = TRUE)
p1 <- plot_cell_trajectory(data$HSMM, color_by = "CellType", show_branch_points = F) + scale_color_manual(values = cluster_cols, name = "CellTypes")
p2 <- plot_cell_trajectory(data$HSMM, color_by = "CellType", show_branch_points = F) + scale_color_manual(values = cluster_cols, name = "CellTypes") + facet_wrap(~CellType) 
p3 <- ggplot(pData(data$HSMM) %>% count(State, CellType), aes(State, n, fill=CellType)) + geom_bar(stat="identity")
print(list(p1, p2, p3))
dev.off()


# Reorder based on the state with the max concentration of the root state
temp = pData(data$HSMM) %>% group_by(CellType) %>% count(State) %>% top_n(n=1, wt=n)
state = as.integer(temp[temp$CellType == expectedRoot, 'State'])
rm(temp)

pData(data$HSMM)$sample <- sample_name
data$HSMM <- orderCells(data$HSMM, root_state = state)
pData(data$HSMM)$RelPseudotime <- pData(data$HSMM)$Pseudotime/max(pData(data$HSMM)$Pseudotime)
pdf(paste(sample_name, num_top_to_use, suffix, sep="_"), paper = 'a4r', onefile = TRUE)
plots <- list(
  plot_cell_trajectory(data$HSMM, color_by = "CellType", show_branch_points = F) + scale_color_manual(values = cluster_cols, name = "CellTypes"),
  plot_cell_trajectory(data$HSMM, color_by = "State", show_branch_points = T),
  plot_cell_trajectory(data$HSMM, color_by = "CellType", show_branch_points = F) + scale_color_manual(values = cluster_cols, name = "CellTypes") + facet_wrap(~CellType),
  ggplot(pData(data$HSMM) %>% count(State, CellType), aes(State, n, fill=CellType)) + geom_bar(stat="identity") + scale_fill_manual(values = cluster_cols, name = "CellTypes"),
  plot_cell_trajectory(data$HSMM, color_by = "RelPseudotime", show_branch_points = F),
  ggplot(pData(data$HSMM), aes(y=RelPseudotime)) + geom_violin(aes(x=sample), color = "red", fill="red", alpha=0.3, trim=FALSE) + coord_flip()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()),
  ggplot(pData(data$HSMM), aes(x=CellType, y=RelPseudotime, fill=CellType)) + geom_violin(trim = FALSE) + scale_fill_manual(values = cluster_cols, name = "CellTypes") + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
)
print(plots)
dev.off()

HSMM <- data$HSMM
save(list=c('cluster_cols', 'num_top_to_use', 'suffix', 'HSMM'), file=paste0('HSMM_', num_top_to_use, '.Robj'))

orderin_exprs <- exprs(data$HSMM)[ordering_genes,]
orderin_taus <- apply(orderin_exprs, MARGIN=1, FUN=function(zz) {
  cor.test(zz, pData(data$HSMM)$RelPseudotime, method = 'kendall')$estimate
})

orderin_rs <- apply(orderin_exprs, MARGIN=1, FUN=function(zz) {
  cor.test(zz, pData(data$HSMM)$RelPseudotime, method = 'spearman')$estimate
})

ordering_cors <- data.frame(orderin_taus, orderin_rs)
ordering_cors$orderin_rs2 <- ordering_cors$orderin_rs**2
colnames(ordering_cors) <- c('tau', 'r', 'r^2')
write.table(file=paste0('pseudotime_vs_gene_exp_', num_top_to_use, '.tsv'), ordering_cors, sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
