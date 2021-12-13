#P1_P10_P11
suppressPackageStartupMessages({
  library(cowplot)
  library(ggplot2)
  library(ggrepel)
  library(grid)
  library(gridExtra)
  library(Seurat)
  library(dplyr)
  library(monocle)
  library(reticulate)
})

setwd('/scRNASeq_workshop/arjun_updated_pdfs_06_17_2020//')


##### ASK 1
### 5C) Could we get a TSNE of the sample aggregates where the clusters are shown in 1 TSNE to 
### match other figures? (Heatmap is in supplemental)

# Old name (key) to new name (value)
name_translation = c(
  "Arg1+" = "Stress-responsive",
  "Ifn-responsive" = "IFN-responsive",
  "Ly6c+ H2-Ab+ Mono" = "Ly6c2+ H2-Ab1+ Mono-Int",
  "Ly6c+ Hp+ Mono" = "Ly6c2+ Hp+ Mono",
  "Ms4a7+ TAM" = "C1qa+ TAM",
  "UK" = "Mgl2+ TAM")

# New name(key) to color (value)
colors <- c(
  "Ly6c2+ H2-Ab1+ Mono-Int" = "deeppink2",
  "IFN-responsive" = "dark blue",
  "Stress-responsive" = "orange1",
  "C1qa+ TAM" = "deeppink4",
  "Ly6c2+ Hp+ Mono" = "lightpink2",
  "Mgl2+ TAM" = "turquoise3")

ordered_names <- c("Ly6c2+ Hp+ Mono",
                   "IFN-responsive",
                   "Stress-responsive",
                   "Ly6c2+ H2-Ab1+ Mono-Int",
                   "Mgl2+ TAM",
                   "C1qa+ TAM")

load('/scRNASeq_workshop/monocle/mouse/SC_integrated_scrubbed.Robj')

SC_integrated$seurat_idents_new_names <- factor(name_translation[as.vector(SC_integrated$seurat_idents)],
                                                levels=ordered_names)

pdf('figure_5C_lower.pdf', height=10, width=10)
TSNEPlot(SC_integrated, group.by="seurat_idents_new_names", cols = colors, pt.size=0.5)
dev.off()

#### Ask 2
### 5D) Is x-axis label OK? Wasn't 100% what scale this was
### RESPONSE: The Yaxis is -log10(pvalue)

#### Ask 3
### 5E) Would it be possible to get the Monocle reordered so that it matches Fig. 1 and goes 
### in descending order: Hp+ mono, IFN-resp, Stress resp, Mono-Int, Mgl2+ TAM (new to Fig 5), 
### C1q TAM?

HSMMs <- list()
for (p in c('P10', 'P11')){
  load(paste0('/scRNASeq_workshop/monocle/mouse/2019_11_15/', p, '_SCT_HSMM.Robj'))
  HSMMs[[p]] <- HSMM
  rm(HSMM)
  pData(HSMMs[[p]])$CellType2 <- factor(name_translation[as.vector(pData(HSMMs[[p]])$CellType)],
                                        levels=rev(ordered_names))
}



for (p in c('P10', 'P11')){
  pdf(paste0('figure_5E_', p, '_trajectory.pdf'), height=10, width=10)
  print(plot_cell_trajectory(HSMMs[[p]], color_by = "CellType2", show_branch_points = F) + 
          scale_color_manual(values = colors, name = "CellTypes"))
  dev.off()

  pdf(paste0('figure_5E_',p, '_violins.pdf'), height=10, width=10)
  print(ggplot(pData(HSMMs[[p]]), aes(x=CellType2, y=RelPseudotime, fill=CellType2)) + 
            geom_violin(trim = FALSE) + 
            scale_fill_manual(values = colors, name = "CellTypes", guide=guide_legend(reverse = TRUE)) + 
            coord_flip() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()))
  dev.off()
}

rm(HSMMs)
#### Ask 4
### For Fig. 3 Suppl, we have a heatmap of the M1/M2 marker genes on P1 Mono clusters 
### (Res06; 0,1,2,3,5). At the time of generating this we had been using a different set 
### up of clusters. My Seurat is acting up and I haven't had a chance to fix it. Since you 
### have been using current P1, would you be able to make a similar heatmap?

M1M2_genes <- read.table('/scRNASeq_workshop/cross_corr_mtx/M1M2_signatures.tsv', sep='\t', header=T, stringsAsFactors = F)
M1M2_genes <- M1M2_genes[order(M1M2_genes$cluster, M1M2_genes$gene), ]

load('/scRNASeq_workshop/monocle/mouse/P1/data.Robj')
names(data)
Idents(data$seurat_data) <- factor(data$seurat_data$res.0.6)

sobj <- subset(data$seurat_data, idents=c("0","1","2","3","5"))

pdf('Supplementary_figure_3A.pdf', height=20, width=20)
DoHeatmap(sobj, features=M1M2_genes$gene)
dev.off()

pdf('Supplementary_figure_3A.pdf', height=20, width=20)
DoHeatmap(sobj, features=M1M2_genes$gene)
dev.off()

genes = apply(data$seurat_data@assays$RNA@scale.data[M1M2_genes$gene[M1M2_genes$gene%in%rownames(data$seurat_data@assays$RNA@scale.data)],], 
            MARGIN = 1, function(x){
                  q1 = quantile(x, 0.1)
                  q3 = quantile(x, 0.9)
                  var(x[x >= q1 & x <= q3])
                  })
genes <- genes[genes > 0.01]
pdf('Supplementary_figure_3A_secondary.pdf', height=20, width=20)
DoHeatmap(sobj, features=names(genes))
dev.off()

rm(data, M1M2_genes, sobj, genes)
#### Ask 5
### And lastly, in Fig. 5 Suppl, I have volcano plots from our original peripheral blood 
### aggregates. It's a minor point so I think it would be fine to keep as is even if 
### generated differently, but was wondering if you would be able to make a volcano plot 
### with the same settings as you did for those in Fig. 5 so they match? If so, attaching 
### the Excel file with the DE analysis. P1 = DTR tumor-bearing blood, P9 = WT tumor-bearing 
### blood

source('volcanize_seurat_markers.R')

# The input data used FindAllMarkers with only.pos=FALSE so teh data is basically duplicated with swapped
# pct.1 and pct.2, and sign flipped logFC
clusters <- list(
    Cluster0 = read.table('WTvsDTRBlood_minpct0.1_logfc0_Cluster0.txt', header=T, sep='\t', stringsAsFactors = F),
    Cluster4 = read.table('WTvsDTRBlood_minpct0.1_logfc0_Cluster4.txt', header=T, sep='\t', stringsAsFactors = F)
    )

for (cc in names(clusters)){
  t1 <- clusters[[cc]][clusters[[cc]]$cluster=='P1', ]
  t2 <- clusters[[cc]][clusters[[cc]]$cluster=='P9', ]
  t2$temp  <- t2$pct.1
  t2$pct.1  <- t2$pct.2
  t2$pct.2  <- t2$temp
  t2$temp <- NULL
  t2$avg_logFC <- -t2$avg_logFC
  t2 <- t2[order(t2$avg_logFC),]
  t1 <- t1[order(t1$avg_logFC),]
  t2$cluster <- NULL
  t1$cluster <- NULL
  
  if (!all(t1 == t2)){
    stop()
  }
}

for (cc in names(clusters)){
  # Since P1 is the test, we'd want the volcano to be what's upregulated in the test.
  clusters[[cc]] <- clusters[[cc]][clusters[[cc]]$cluster=='P1',]
  clusters[[cc]]$cluster <- NULL
  rownames(clusters[[cc]]) <- clusters[[cc]]$gene
}

for (cc in names(clusters)){
  pdf(paste0('Supplementary_figure_5E_', cc, '.pdf'), height=10, width=10)
  print(volcanize_from_FindMarkers(clusters[[cc]], pos_group ="DTR tumor-bearing blood", neg_group = "WT tumor-bearing blood",
                                   background_color = 'black', text_color = 'red'))
  dev.off()
  pdf(paste0('Supplementary_figure_5E_', cc, '_secondary.pdf'), height=10, width=10)
  print(volcanize_from_FindMarkers(clusters[[cc]], pos_group ="DTR tumor-bearing blood", neg_group = "WT tumor-bearing blood",
                                   background_color = 'darkgrey'))
  dev.off()
}

rm(clusters)

### Not asks but redone volcanoes

P10_P11.sobj <- subset(SC_integrated, idents=c('P10', 'P11'))

P10_P11.sobj@active.assay <- 'SCT'

Idents(P10_P11.sobj) <- P10_P11.sobj@meta.data$seurat_idents_new_names

for (cluster in levels(P10_P11.sobj$seurat_idents_new_names)){
  print(paste0('processing ', cluster))
  temp <- subset(P10_P11.sobj, idents=c(cluster))
  Idents(temp) <- temp@meta.data$orig.ident
  markers <- FindMarkers(temp, ident.1='P11', ident.2='P10', assay = 'SCT', logfc.threshold = 0.001, min.pct = 0.1, test.use='MAST', only.pos=FALSE)
  
  pdf(paste0('Figure_5D_', gsub(" " , "_", gsub('+', 'pos', cluster, fixed = TRUE)), ".pdf"))
  print(volcanize_from_FindMarkers(markers, pos_group = "P11", neg_group = 'P10', background_color = 'black', 
                             text_color='red'))
  dev.off()

  pdf(paste0('Figure_5D_', gsub(" " , "_", gsub('+', 'pos', cluster, fixed = TRUE)), "_secondary.pdf"))
  print(volcanize_from_FindMarkers(markers, pos_group = "P11", neg_group = 'P10', background_color = 'darkgrey'))
  dev.off()
}



