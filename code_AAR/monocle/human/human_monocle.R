library(dplyr)
library(ggplot2)
library(Seurat)
library(monocle3)


# Mad props to https://github.com/zihengxuwu
# https://github.com/satijalab/seurat/issues/1658#issuecomment-513920073

setwd('/krummellab//data1/arrao/projects/adriana_paper_tipcc/human_monocle/with_SEPP')
load('../Final_object_mono_mac_merged_V2.Robj')

DimPlot(merged_data_subset_mono_Mac_4)
View(merged_data_subset_mono_Mac_4@meta.data)


colors <- setNames(
      c("lightpink2",
        "deeppink2",
        "Dark blue",
        "orange1",
        "turquoise4",
        "deeppink4"),
      c("CD14+ monocytes",
        "CD14+ mono-int",
        "IFNresp-mono",
        "Stress Response SPP1",
        "SEPP1 TAM",
        "C1Q+ TAM")
      )


cds_list <- list()

gene_annotation <- as.data.frame(rownames(merged_data_subset_mono_Mac_4@reductions[["pca"]]@feature.loadings), 
                                 row.names = rownames(merged_data_subset_mono_Mac_4@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

# part two, cell information
stopifnot(all(colnames(merged_data_subset_mono_Mac_4@assays$RNA) == rownames(merged_data_subset_mono_Mac_4@meta.data)))

cell_metadata <- merged_data_subset_mono_Mac_4@meta.data


# part three, counts sparse matrix
expression_matrix <- merged_data_subset_mono_Mac_4@assays[["RNA"]]@counts
expression_matrix <- expression_matrix[rownames(merged_data_subset_mono_Mac_4@reductions[["pca"]]@feature.loadings), ]

### Construct the basic cds object

cds_obj <- new_cell_data_set(expression_matrix,
                             cell_metadata = cell_metadata,
                             gene_metadata = gene_annotation)
rm(expression_matrix, cell_metadata, gene_annotation)

### Construct and assign the made up partition

recreate.partition <- c(rep(1, length(cds_obj@colData@rownames)))
names(recreate.partition) <- cds_obj@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_obj@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
rm(recreate.partition)

### Assign the cluster info

list_cluster <- as.vector(merged_data_subset_mono_Mac_4$Final_Annotations)
names(list_cluster) <- merged_data_subset_mono_Mac_4@assays[["RNA"]]@data@Dimnames[[2]]

cds_obj@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
rm(list_cluster)


### Could be a space-holder, but essentially fills out louvain parameters

cds_obj@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


### Assign UMAP coordinate

cds_obj@int_colData@listData[["reducedDims"]][["UMAP"]] <- merged_data_subset_mono_Mac_4@reductions[["umap"]]@cell.embeddings



### Assign feature loading for downstream module analysis

cds_obj@preprocess_aux$gene_loadings <- merged_data_subset_mono_Mac_4@reductions[["pca"]]@feature.loadings

### Learn graph, this step usually takes a significant period of time for larger samples

print("Learning graph, which can take a while depends on the sample")

cds_obj <- learn_graph(cds_obj, use_partition = T)

pdf("clusters.with.trajectory.pdf", width = 10, height = 10, useDingbats=F)
clus <- plot_cells(cds_obj, 
                   color_cells_by = 'cluster',
                   label_groups_by_cluster=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE)
print(clus)
dev.off()

pdf("orig.ident.with.trajectory.pdf", width = 10, height = 10, useDingbats=F)
clus <- plot_cells(cds_obj, 
                   color_cells_by = 'orig.ident',
                   label_groups_by_cluster=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE)
print(clus)
dev.off()

write.table(merged_data_subset_mono_Mac_4@reductions[["umap"]]@cell.embeddings,
            file='umap_coords.tsv',
            sep='\t',
            col.names=T, 
            row.names=T)

root_cells <- colnames(merged_data_subset_mono_Mac_4)[merged_data_subset_mono_Mac_4$Final_Annotations == "CD14+ monocytes"]
root_cells <- colnames(merged_data_subset_mono_Mac_4)[merged_data_subset_mono_Mac_4@reductions$umap@cell.embeddings[,1] < -5]

cds_obj <- order_cells(cds_obj, root_cells = root_cells)

### Now we plot pseudotime

print("Plotting pseudotime")

pdf("pseudotime.with.trajectory.pdf", width = 10, height = 10, useDingbats=F)
ptime <- plot_cells(cds_obj, 
                    color_cells_by = 'pseudotime',
                    label_groups_by_cluster=FALSE,
                    label_leaves=FALSE,
                    label_branch_points=FALSE)
print(ptime)
dev.off()

df <- data.frame(cds_obj@principal_graph_aux$UMAP$pseudotime,
                 factor(as.vector(cds_obj@clusters@listData[["UMAP"]][["clusters"]]), levels=rev(names(colors))),
                 pData(cds_obj)[["orig.ident"]])

colnames(df) <- c('pseudotime', 'clusters', 'orig.ident')

df$rel_pseudotime <- (df$pseudotime-min(df$pseudotime))/(max(df$pseudotime)-min(df$pseudotime))
df$indication <- gsub("^IPI", "", gsub("[0-9]{3}[.].*", "", df$orig.ident))

pdf("diff_pseudotime_cluster_violins.pdf", width = 10, height = 10, useDingbats=F)
dptime <- ggplot(df, aes(x=rel_pseudotime, y=clusters, fill=clusters)) + 
              geom_violin(scale="width") +
              scale_fill_manual(values=colors) +
              theme_bw()
print(dptime)
dev.off()

pdf("diff_pseudotime_cluster_violins_by_indication.pdf", width = 15, height = 10, useDingbats=F)
dptime <- ggplot(df, aes(x=rel_pseudotime, y=clusters, fill=clusters)) + 
              geom_violin(scale="width") +
              scale_fill_manual(values=colors) +
              theme_bw() + 
              facet_wrap(~indication)
print(dptime)
dev.off()

 
rm(merged_data_subset_mono_Mac_4, dptime, ptime, clus)
save(list=ls(), file='workspace.RData')