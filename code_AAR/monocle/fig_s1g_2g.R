#P1
suppressPackageStartupMessages({
  library(cowplot)
  library(gridExtra)
  library(Seurat)
  library(reshape)
  library(dplyr)
  library(monocle)
})


#setwd('/krummellab/data1/arrao/projects/adriana_paper/monocle/mouse/P1')
setwd('/Users/arjunarkalrao/projects/adriana_paper/monocle/mouse/P1')
name_translation <- read.table('new_names.tsv', header=T, sep='\t', row.names=1, colClasses = c('character'))
load('HSMM_all.Robj')
p1 <- plot_cell_trajectory(HSMM, markers = c('H2-Ab1'), use_color_gradient = TRUE, show_branch_points=F, cell_size=0.5) + scale_colour_gradient(low='blue', high='red')
p2 <- plot_cell_trajectory(HSMM, markers = c('Fcgr1'), use_color_gradient = TRUE, show_branch_points=F, cell_size=0.5) + scale_colour_gradient(low='blue', high='red')
p3 <- plot_cell_trajectory(HSMM, markers = c('Ly6c2'), use_color_gradient = TRUE, show_branch_points=F, cell_size=0.5) + scale_colour_gradient(low='blue', high='red')

pdf('fig_s1g.pdf', width = 10, height = 15)
print(CombinePlots(plots = list(p1, p2, p3), ncol=1))
dev.off()


load('data.Robj')
top_markers <- data$markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
x <- data$seurat_data@assays$RNA@scale.data[unique(top_markers$gene), colnames(HSMM)]
x <- rbind(x, pData(HSMM)$RelPseudotime)
rownames(x)[26] = 'RelPseudotime'

x <- as.data.frame(t(x))

Molten <- melt(x, id.vars = "RelPseudotime")
cluster_mapping <- read.table('P1_clusters.tsv', sep="\t", row.names=1, header=TRUE)

y <- as.data.frame((data$markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))[, c('gene', 'cluster')])
rownames(y)<- y$gene
y$lty=rep(c("solid", "dashed", "dotted", "dotdash", "twodash"), 5)
y$color=sapply(y$cluster, function(z){cluster_mapping[cluster_mapping$P1==as.character(z), 'colors', drop=T]})
#### FIX factors in y
y$color <- as.character(y$color)
y$cluster <- as.character(y$cluster)
y$cluster2 <- sapply(y$cluster, function(x){as.character(name_translation[x, 'new'])})


Molten$cluster <- sapply(Molten$variable, function(z){y[as.character(z), 'cluster2', drop=T]})
Molten$color <- sapply(Molten$variable, function(z){y[as.character(z), 'color', drop=T]})
Molten$lty <- sapply(Molten$variable, function(z){y[as.character(z), 'lty', drop=T]})
#p1 <-  ggplot(Molten, aes(x = RelPseudotime, y=value, colour = color, fill= color, linetype = lty)) + 
#ggplot(Molten, aes(x = RelPseudotime, y=value, lty=lty, color=color, fill=color)) + 
#    geom_smooth(method='loess') + 
#    scale_color_manual(values=levels(as.factor(Molten$color))) +
#    scale_fill_manual(values=levels(as.factor(Molten$color))) +
#    guides(color = guide_legend(order = 1),  # Hides the "color" legend because it has 4 values and we don't want that.
#           fill = guide_legend(order = 1),  # Forces the "Values" legend group to be on top
#           lty = guide_legend(order = 1)) + 
#    facet_wrap(~cluster)

dat <- split(Molten, f = Molten$cluster)
plots <- list()
for (ptype in c("Ly6c2+ Hp+ Mono", "IFN-responsive", "Stress-responsive", "Mono-Int", "C1qa+ TAM")) {
  plots[[ptype]] <- ggplot(dat[[ptype]]) +
      geom_smooth(method='loess', aes(x = RelPseudotime, y=value, color=color, fill=color, lty=lty)) +
      scale_color_manual("", values=y[y$cluster2==ptype, 'color'], guide = FALSE) +
      scale_fill_manual("", values=y[y$cluster2==ptype, 'color'], guide = FALSE) +
      scale_linetype_manual("", label = y[y$cluster2==ptype, 'gene'], values = y[y$cluster2==ptype, 'lty'], guide = FALSE) +
      guides(linetype = guide_legend(title=paste0(ptype, ' genes'), override.aes = list(color = y[y$cluster2==ptype, 'color'], fill = y[y$cluster2==ptype, 'color'], linetype = y[y$cluster2==ptype, 'lty'])))
}

both <- plot_grid(plotlist = plots, ncol=1, align = "v")
save_plot("fig_1i.pdf", both, base_height = 20, base_width = 10)



p1 <- plot_cell_trajectory(HSMM, markers = c('Arg1'), use_color_gradient = TRUE, show_branch_points=F, cell_size=0.5) + scale_colour_gradient(low='blue', high='red')
p2 <- plot_cell_trajectory(HSMM, markers = c('C1qa'), use_color_gradient = TRUE, show_branch_points=F, cell_size=0.5) + scale_colour_gradient(low='blue', high='red')
p3 <- plot_cell_trajectory(HSMM, markers = c('Mgl2'), use_color_gradient = TRUE, show_branch_points=F, cell_size=0.5) + scale_colour_gradient(low='blue', high='red')
p4 <- plot_cell_trajectory(HSMM, markers = c('Il7r'), use_color_gradient = TRUE, show_branch_points=F, cell_size=0.5) + scale_colour_gradient(low='blue', high='red')
p5 <- plot_cell_trajectory(HSMM, markers = c('Vcam1'), use_color_gradient = TRUE, show_branch_points=F, cell_size=0.5) + scale_colour_gradient(low='blue', high='red')
p6 <- plot_cell_trajectory(HSMM, markers = c('Mki67'), use_color_gradient = TRUE, show_branch_points=F, cell_size=0.5) + scale_colour_gradient(low='blue', high='red')
p7 <- plot_cell_trajectory(HSMM, markers = c('Birc5'), use_color_gradient = TRUE, show_branch_points=F, cell_size=0.5) + scale_colour_gradient(low='blue', high='red')
p8 <- plot_cell_trajectory(HSMM, markers = c('Cxcl10'), use_color_gradient = TRUE, show_branch_points=F, cell_size=0.5) + scale_colour_gradient(low='blue', high='red')
p9 <- plot_cell_trajectory(HSMM, markers = c('Adgre1'), use_color_gradient = TRUE, show_branch_points=F, cell_size=0.5) + scale_colour_gradient(low='blue', high='red')

pdf('fig_2g.pdf', width = 10, height = 45)
print(CombinePlots(plots = list(p1, p2, p3, p4, p5, p6 ,p7, p8, p9), ncol=1))
dev.off()

pdf('fig_2g_separate.pdf', width = 20, height = 10, onefile = T)
plot_cell_trajectory(HSMM, markers = c('Arg1'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red')
plot_cell_trajectory(HSMM, markers = c('C1qa'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red')
plot_cell_trajectory(HSMM, markers = c('Mgl2'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red')
plot_cell_trajectory(HSMM, markers = c('Il7r'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red')
plot_cell_trajectory(HSMM, markers = c('Vcam1'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red')
plot_cell_trajectory(HSMM, markers = c('Mki67'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red')
plot_cell_trajectory(HSMM, markers = c('Birc5'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red')
plot_cell_trajectory(HSMM, markers = c('Cxcl10'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red')
plot_cell_trajectory(HSMM, markers = c('Adgre1'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red')
dev.off()



HSMM$CellType2  = sapply(as.character(HSMM$CellType), function(x) {name_translation[x, 'new']})
HSMM$CellType2  = factor(HSMM$CellType2, levels = c("C1qa+ TAM", "Mono-Int", "Stress-responsive", "IFN-responsive", "Ly6c2+ Hp+ Mono" ))

cluster_cols2 = cluster_cols[names(cluster_cols)!='Rest']
names(cluster_cols2) <- sapply(as.character(names(cluster_cols2)), function(x) {name_translation[x, 'new']})

ggplot(pData(HSMM), aes(x=CellType2, y=RelPseudotime, fill=CellType2)) + 
    geom_violin(trim = FALSE) + 
    scale_fill_manual(values = cluster_cols2, name = "CellTypes") + 
    coord_flip() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

# unused attempt at contour marking HSMM plots
# plot_cell_trajectory(HSMM, color_by = "CellType", show_branch_points = F) + stat_density_2d(aes(fill = stat(level)), geom = "polygon") + scale_color_manual(values = cluster_cols, name = "CellTypes") + facet_wrap(~CellType)
