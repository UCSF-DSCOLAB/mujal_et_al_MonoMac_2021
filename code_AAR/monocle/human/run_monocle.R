#HUMAN
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


setwd('~/projects/adriana_paper/monocle/human/')
#setwd('/data/shared/krummellab/arrao/projects/adriana_paper/monocle/human/')
cluster_mapping <- data.frame(
  rbind(
    c(0,   "LYPD3+ Monocytes", "deeppink2"),
    c(1,   "S100A12+ Monocytes", "lightpink2"),
    c(2,   "IFN-responsive Monocytes", "Dark blue"),
    c(3,   "APOE+ Macrophages", "deeppink4"),
    c(4,   "SEPP1+ Monos/Macs", "orange1"),
    c(999, "Rest", "black")  # remaining
  ),
  row.names = "X1"
)
colnames(cluster_mapping) <- c('MonoMac', 'colors')
cluster_cols <- as.character(cluster_mapping$colors)
names(cluster_cols) <- as.character(cluster_mapping$MonoMac)

########################################

if (!file.exists('data.Robj')){
  load('Mono_Mac.Robj')
  
  # Convert the seurat 2.x object into a Seurat 3.x object
  data = list(seurat_data= UpdateSeuratObject(H_KID_TUM_MY_0.4_Small_MONO_MACS))
  rm(H_KID_TUM_MY_0.4_Small_MONO_MACS)
  
  data$seurat_data@meta.data$cluster <- Idents(data$seurat_data)
  data$seurat_data@meta.data$cell_type <- as.character(sapply(data$seurat_data@meta.data$cluster, function(x) {cluster_mapping[x, 'MonoMac', drop=TRUE]}))
  Idents(data$seurat_data) <- as.factor(as.character(sapply(data$seurat_data@meta.data$cluster, function(x) {cluster_mapping[x, 'MonoMac', drop=TRUE]})))
  data$raw <- as.matrix(data$seurat_data@assays$RNA@counts)
  data$raw <- data$raw[,colnames(data$raw)%in%colnames(data$seurat_data@assays$RNA@scale.data)]
  data$markers <- FindAllMarkers(data$seurat_data, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.35, test.use='MAST')
  data$top_marker_per_cell <- data$markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
  
  data$most_distant <- SubsetData(data$seurat_data, subset.name="cluster", accept.value=c("1", "3"))
  data$most_distant_markers <- FindAllMarkers(data$most_distant, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.35, test.use='MAST')
  data$top_most_distant_marker <- data$most_distant_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
  
  DimPlot(data$seurat_data, pt.size=3)
  ggsave(file=paste('MonoMac', "seurat_tsne.png", sep="_"))
  
  FeaturePlot(data$seurat_data, features = data$top_marker_per_cell$gene, pt.size=1.5)
  ggsave(file=paste('MonoMac', "seurat_top_genes_per_cluster.png", sep="_"))
  
  DimPlot(data$most_distant, pt.size=3)
  ggsave(file=paste('MonoMac', "md_seurat_tsne.png", sep="_"))
  
  FeaturePlot(data$most_distant, features = data$top_most_distant_marker$gene, pt.size=1.5)
  ggsave(file=paste('MonoMac', "md_seurat_top_genes_per_cluster.png", sep="_"))
  
  data$raw <- data$raw[order(rownames(data$raw)),]
  
  clusters <- as.character(unique(data$top_marker_per_cell$cluster))
  
  cluster_ids <- as.character(data$seurat_data@meta.data[colnames(data$raw), "cluster"])
  cluster_ids <- data.frame(
    row.names = colnames(data$raw),
    cluster = cluster_ids,
    CellType = factor(sapply(cluster_ids, function(x) {cluster_mapping[x, 'MonoMac']}), levels=clusters)
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
  
  data$HSMM.bk <- data$HSMM
  
  data$HSMM <- data$HSMM.bk
  
  save(list=c('data'), file='data.Robj')
  quit(save = 'no')
} else {
  load('data.Robj')
}

if (TRUE){
  num_top_to_use <- 500
  top_genes <- as.data.frame(data$markers %>% group_by(cluster) %>% top_n(n = num_top_to_use, wt = avg_logFC))
  num_top_to_use <- 'all'
  suffix <- "genes_trajectories.pdf"
  clusters <- as.character(unique(data$top_marker_per_cell$cluster))
} else {
  num_top_to_use <- 500
  top_genes <- as.data.frame(data$most_distant_markers %>% group_by(cluster) %>% top_n(n = num_top_to_use, wt = avg_logFC))
  num_top_to_use <- 'all'
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

pdf(paste('MonoMac', num_top_to_use, suffix, sep="_"), paper = 'a4r', onefile = TRUE)
p1 <- plot_cell_trajectory(data$HSMM, color_by = "CellType", show_branch_points = F) + scale_color_manual(values = cluster_cols, name = "CellTypes")
p2 <- plot_cell_trajectory(data$HSMM, color_by = "CellType", show_branch_points = F) + scale_color_manual(values = cluster_cols, name = "CellTypes") + facet_wrap(~CellType) 
p3 <- ggplot(pData(data$HSMM) %>% count(State, CellType), aes(State, n, fill=CellType)) + geom_bar(stat="identity")
print(list(p1, p2, p3))
dev.off()


# Reorder based on the state with the max concentration of S100A12+ Monocytes
temp = pData(data$HSMM) %>% group_by(CellType) %>% count(State) %>% top_n(n=1, wt=n)
state = as.integer(temp[temp$CellType == 'S100A12+ Monocytes', 'State'])
rm(temp)

pData(data$HSMM)$sample <- 'MonoMac'
data$HSMM <- orderCells(data$HSMM, root_state = state)
pData(data$HSMM)$RelPseudotime <- pData(data$HSMM)$Pseudotime/max(pData(data$HSMM)$Pseudotime)
pdf(paste('MonoMac', num_top_to_use, suffix, sep="_"), paper = 'a4r', onefile = TRUE)
plots <- list(
  plot_cell_trajectory(data$HSMM, color_by = "CellType", show_branch_points = F) + scale_color_manual(values = cluster_cols, name = "CellTypes"),
  plot_cell_trajectory(data$HSMM, color_by = "State", show_branch_points = T),
  plot_cell_trajectory(data$HSMM, color_by = "CellType", show_branch_points = F) + scale_color_manual(values = cluster_cols, name = "CellTypes") + facet_wrap(~CellType),
  ggplot(pData(data$HSMM) %>% count(State, CellType), aes(State, n, fill=CellType)) + scale_fill_manual(values = cluster_cols, name = "CellTypes") + geom_bar(stat="identity"),
  plot_cell_trajectory(data$HSMM, color_by = "RelPseudotime", show_branch_points = F),
  ggplot(pData(data$HSMM), aes(y=RelPseudotime)) + geom_violin(aes(x=sample), color = "red", fill="red", alpha=0.3, trim=FALSE) + coord_flip()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()),
  ggplot(pData(data$HSMM), aes(y=RelPseudotime)) + geom_violin(aes(x=CellType), trim = FALSE) + coord_flip() + facet_wrap(~sample, ncol=1, drop=TRUE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
)
print(plots)
dev.off()

HSMM <- data$HSMM
save(list=c('HSMM'), file='HSMM.Robj')

ordering_exprs <- exprs(data$HSMM)[ordering_genes,]
ordering_taus <- apply(ordering_exprs, MARGIN=1, FUN=function(zz) {
  cor.test(zz, pData(data$HSMM)$RelPseudotime, method = 'kendall')$estimate
})

ordering_rs <- apply(ordering_exprs, MARGIN=1, FUN=function(zz) {
  cor.test(zz, pData(data$HSMM)$RelPseudotime, method = 'spearman')$estimate
})

ordering_cors <- data.frame(ordering_taus, ordering_rs)
ordering_cors$ordering_rs2 <- ordering_cors$ordering_rs**2
colnames(ordering_cors) <- c('tau', 'r', 'r^2')
write.table(file="pseudotime_vs_gene_exp.tsv", ordering_cors, sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


plot_genes_in_pseudotime(HSMM[c('S100A8', 'S100A12')], color_by = "CellType")
