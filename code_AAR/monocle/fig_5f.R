#P10 P11
suppressPackageStartupMessages({
  library(cowplot)
  library(gridExtra)
  library(Seurat)
  #library(reshape)
  library(dplyr)
  library(monocle)
})


#setwd('/Users/arjunarkalrao/projects/adriana_paper/monocle/mouse/')
setwd('/data')

cluster_mapping <- read.table('mouse/P1/P1_clusters.tsv', sep="\t", row.names=1, header=TRUE)
cluster_cols <- as.character(cluster_mapping$colors)
names(cluster_cols) <-  as.character(cluster_mapping[["P1"]])
cluster_cols['UK'] <- 'black'

new_celltype_names = c(
  "Ly6c+ H2-Ab+ Mono"= "Ly6c2+ H2-Ab1+ Mono-Int",
  "Ifn-responsive" ="Ifn-responsive",
  "Arg1+" = "Stress-responsive",
  "Ms4a7+ TAM" = "C1qa+ TAM",
  "Ly6c+ Hp+ Mono" = "Ly6c2+ Hp+ Mono",
  "UK" = "Mgl2+ TAM",
  "Rest" = "Rest"
)
new_cluster_cols = cluster_cols
names(new_cluster_cols) <- new_celltype_names[names(new_cluster_cols)]

new_cluster_cols['Mgl2+ TAM'] = 'turquoise3'

load('mouse/P10_SCT_HSMM.Robj')
stopifnot(!'Rest' %in% levels(HSMM$CellType))
HSMM$CellType2 <- factor(new_celltype_names[as.vector(HSMM$CellType)],
                         levels=c("Ly6c2+ Hp+ Mono",
                                  "Ifn-responsive",
                                  "Ly6c2+ H2-Ab1+ Mono-Int",
                                  "Stress-responsive",
                                  "Mgl2+ TAM",
                                  "C1qa+ TAM"))

pdf('fig_5f_P10_separate.pdf', width = 20, height = 10, onefile = T)
plot_cell_trajectory(HSMM, color_by = 'CellType2', show_branch_points=F, cell_size=2.0) + scale_color_manual(values=new_cluster_cols) + scale_x_continuous(trans = "reverse") + theme_void() + theme(legend.position = "none")
plot_cell_trajectory(HSMM, markers = c('C1qc'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('C1qa'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Arg1'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Mgl2'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Apoe'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('H2-Ab1'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('H2-Aa'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Hp'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Ly6c2'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Chil3'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Ms4a7'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
dev.off()

pdf('fig_5f_P10_separate_grey_red.pdf', width = 20, height = 10, onefile = T)
plot_cell_trajectory(HSMM, color_by = 'CellType2', show_branch_points=F, cell_size=2.0) + scale_color_manual(values=new_cluster_cols) + scale_x_continuous(trans = "reverse") + theme_void() + theme(legend.position = "none")
plot_cell_trajectory(HSMM, markers = c('C1qc'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('C1qa'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Arg1'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Mgl2'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Apoe'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('H2-Ab1'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('H2-Aa'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Hp'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Ly6c2'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Chil3'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Ms4a7'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
dev.off()

pdf('fig_5f_P10_new_dpt.pdf', width = 10, height = 20, onefile = T)
ggplot(pData(HSMM), aes(y=RelPseudotime, x=CellType2, fill=CellType2)) + 
  geom_violin(trim = F) + 
  scale_x_discrete(limits = rev(levels(pData(HSMM)$CellType2))) + 
  scale_fill_manual(values=new_cluster_cols) +
  coord_flip() + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

load('mouse/P11_SCT_HSMM.Robj')
stopifnot(!'Rest' %in% levels(HSMM$CellType))
HSMM$CellType2 <- factor(new_celltype_names[as.vector(HSMM$CellType)],
                         levels=c("Ly6c2+ Hp+ Mono",
                                  "Ifn-responsive",
                                  "Ly6c2+ H2-Ab1+ Mono-Int",
                                  "Stress-responsive",
                                  "Mgl2+ TAM",
                                  "C1qa+ TAM"))

pdf('fig_5f_P11_separate.pdf', width = 20, height = 10, onefile = T)
plot_cell_trajectory(HSMM, color_by = 'CellType2', show_branch_points=F, cell_size=2.0) + scale_color_manual(values=new_cluster_cols) + scale_x_continuous(trans = "reverse") + theme_void() + theme(legend.position = "none")
plot_cell_trajectory(HSMM, markers = c('C1qc'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('C1qa'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Arg1'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Mgl2'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Apoe'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('H2-Ab1'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('H2-Aa'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Hp'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Ly6c2'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Chil3'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Ms4a7'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='blue', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
dev.off()

pdf('fig_5f_P11_separate_grey_red.pdf', width = 20, height = 10, onefile = T)
plot_cell_trajectory(HSMM, color_by = 'CellType2', show_branch_points=F, cell_size=2.0) + scale_color_manual(values=new_cluster_cols) + scale_x_continuous(trans = "reverse") + theme_void() + theme(legend.position = "none")
plot_cell_trajectory(HSMM, markers = c('C1qc'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('C1qa'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Arg1'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Mgl2'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Apoe'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('H2-Ab1'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('H2-Aa'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Hp'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Ly6c2'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Chil3'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
plot_cell_trajectory(HSMM, markers = c('Ms4a7'), use_color_gradient = TRUE, show_branch_points=F, cell_size=2.0) + scale_colour_gradient(low='lightgrey', high='red') + scale_x_continuous(trans = "reverse") + theme_void() 
dev.off()

pdf('fig_5f_P11_new_dpt.pdf', width = 10, height = 20, onefile = T)
ggplot(pData(HSMM), aes(y=RelPseudotime, x=CellType2, fill=CellType2)) + 
  geom_violin(trim = F) + 
  scale_x_discrete(limits = rev(levels(pData(HSMM)$CellType2))) + 
  scale_fill_manual(values=new_cluster_cols) +
  coord_flip() + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()