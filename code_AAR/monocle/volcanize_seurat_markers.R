library(ggrepel)

# No guarantees this works... needs testing
volcanize_from_FindAllMarkers <- function(x, pos_group, neg_group, label_n=10, logFC_threshold=0.4, padj_threshold=0.05, background_color='lightgrey', sig_color='red', text_color='black', point_size=0.5) {
  x <- x[x$cluster %in% c(pos_group, neg_group), ]
  x <- transform(x, 
                 plmin_log_FC=ifelse(cluster==pos_group, 
                                     avg_logFC, 
                                     -avg_logFC)
  )
  top_genes <- (x %>% 
                  group_by(cluster) %>% 
                  filter(avg_logFC >= logFC_threshold & p_val_adj <= padj_threshold) %>% 
                  top_n(label_n, wt=avg_logFC))$gene
  x <- transform(x, col=factor(ifelse(gene %in% top_genes, 
                                      2, 
                                      ifelse(avg_logFC>logFC_threshold & p_val_adj <= padj_threshold, 
                                             1, 
                                             0)), 
                               levels=c(0,1,2)))  
  
  plot <- ggplot(x, aes(x=plmin_log_FC, y=-log10(p_val_adj))) + 
    geom_point(aes(col=col), size=point_size) + 
    scale_color_manual(values=c(background_color, 
                                sig_color, 
                                sig_color), drop=F) + 
    geom_text_repel(data=filter(x, col=='2'), 
                    aes(label=gene), 
                    size=5, 
                    color=text_color) +
    theme_minimal() +
    theme(legend.position = "none") +
    xlab('logFC') +
    xlim(-max(abs(x$avg_logFC)),
         max(abs(x$avg_logFC))) +
    ggtitle(label = paste0('Volcano of ', pos_group, ' vs. ', neg_group)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  if (any(x$p_val_adj == 0)) {
    pbuild <- ggplot_build(plot = plot)
    y.min <- pbuild$layout$panel_params[[1]]$y.range[1]
    y.max <- pbuild$layout$panel_params[[1]]$y.range[2]
    
    plot <- plot + scale_y_continuous(breaks=c(pbuild$layout$panel_params[[1]]$y$minor_breaks, y.max),
                                      labels=c(pbuild$layout$panel_params[[1]]$y$minor_breaks, expression(infinity))) +
      theme(axis.text.y = element_text(size = 15),
            axis.text.x = element_text(size = 15),
            plot.margin=unit(c(30, 30, 5.5, 5.5), "points"))
  }
  return(plot)
}

# This works great!
volcanize_from_FindMarkers <- function(x, pos_group='Positive_group', neg_group='Negative_group', 
                                       label_n=10, logFC_threshold=0.4, padj_threshold=0.05, 
                                       logpadj_visual_cutoff=NULL, background_color='lightgrey', 
                                       sig_color='red', text_color='black', point_size=0.5) {
  x <- transform(x, 
                 regulation=ifelse(avg_logFC<0, 'DOWN', 'UP'),
                 abs_log_FC = abs(avg_logFC)
  )
  x$gene <- rownames(x)
  
  top_genes <- (x %>% 
                  group_by(regulation) %>% 
                  filter(abs_log_FC >= logFC_threshold & p_val_adj <= padj_threshold) %>% 
                  top_n(label_n, wt=abs_log_FC))$gene
  
  if (!is.null(logpadj_visual_cutoff)){
    padj_visual_cutoff = 10**(logpadj_visual_cutoff)
    x$p_val_adj[x$p_val_adj<=padj_visual_cutoff] = 0
  }
  
  x <- transform(x, col=factor(ifelse(gene %in% top_genes, 
                                      2, 
                                      ifelse(abs_log_FC>logFC_threshold & p_val_adj <= padj_threshold, 
                                             1, 
                                             0)), 
                               levels=c(0,1,2)))
  
  plot <- ggplot(x, aes(avg_logFC, -log10(p_val_adj))) + 
    geom_point(aes(col=col), size=point_size) + 
    scale_color_manual(values=c(background_color, 
                                sig_color, 
                                sig_color), drop=F) + 
    geom_text_repel(data=filter(x, col=='2'), 
                    aes(label=gene), 
                    size=5, 
                    color=text_color) +
    theme_minimal() +
    theme(legend.position = "none") +
    xlab('logFC') +
    xlim(-max(x$abs_log_FC),
         max(x$abs_log_FC)) +
    ggtitle(label = paste0('Volcano of ', pos_group, ' vs. ', neg_group)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  if (any(x$p_val_adj == 0)) {
    pbuild <- ggplot_build(plot = plot)
    y.min <- pbuild$layout$panel_params[[1]]$y.range[1]
    y.max <- pbuild$layout$panel_params[[1]]$y.range[2]
    
    if (is.null(logpadj_visual_cutoff)){
      extreme_string <- expression(infinity)
    } else {
      extreme_string <- paste0("> ", logpadj_visual_cutoff)
    }
    
    # Ggplot 3.3.0 appears to have made an internal change to panel_params. Previously, the minor and makor breaks
    # on the y axis were stored as y.minor_source and y.major_source but now there's a ggproto object for the x and y
    # axes and the breaks are stored in y$minor_breaksbreaks and y$breaks
    if (is.null(pbuild$layout$panel_params[[1]]$y)){
      # GGplot 3.2.1
      breaks <- pbuild$layout$panel_params[[1]]$y.major_source
    } else {
      # GGplot 3.3.0
      breaks <- pbuild$layout$panel_params[[1]]$y$breaks
    }
    
    if (is.null(breaks)){
      stop('Couldn\'t resolve breaks in the plot. This is related to ggplot version')
    }
    
    plot <- plot + 
              scale_y_continuous(breaks=c(breaks, y.max),
                                 labels=c(breaks, extreme_string)) +
              theme(axis.text.y = element_text(size = 15),
                    axis.text.x = element_text(size = 15),
                    plot.margin=unit(c(30, 30, 5.5, 5.5), "points"))
  }
  return(plot)
}
