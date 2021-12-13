# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
library(ggplot2)
library(grid)
library(reshape2)

rotate <- function(x) t(apply(x, 2, rev))

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}


frac = 'no_dcs'

# Process Both
m1m2_genes <- m1m2$gene
gene_df <- get_genes(m1m2_genes, data[[frac]], drop=TRUE)

cormat = cor(t(gene_df))
upper_tri <- get_upper_tri(cormat)

# Get the values in each matrix for plotting the hist
m1_mat <- upper_tri[rownames(upper_tri) %in% m1_genes, colnames(upper_tri) %in% m1_genes]
m1_mat[lower.tri(m1_mat, diag=TRUE)]<- NA
m1_mat <- as.vector(m1_mat)
m1_mat <- m1_mat[!is.na(m1_mat)]

m2_mat <- upper_tri[rownames(upper_tri) %in% m2_genes, colnames(upper_tri) %in% m2_genes]
m2_mat[lower.tri(m2_mat, diag=TRUE)]<- NA
m2_mat <- as.vector(m2_mat)
m2_mat <- m2_mat[!is.na(m2_mat)]

m1m2_mat <- upper_tri[rownames(upper_tri) %in% m2_genes, colnames(upper_tri) %in% m1_genes]

# Now rotate it
upper_tri <- rotate(rotate(rotate(upper_tri)))

melted_cormat <- melt(upper_tri, na.rm = TRUE)


len = length(row.names(upper_tri))
last_m2 = sum(row.names(upper_tri) %in% m1m2[m1m2$cluster=='M2', 'gene'])
my.lines<-data.frame(x=c(0.5, last_m2+0.5),
                     y=c(len-last_m2+0.5, len-last_m2+0.5), 
                     xend=c(last_m2+0.55, last_m2+0.5),
                     yend=c(len-last_m2+0.55, 0.5))
my.labels <- data.frame(x=c(last_m2 + 0.7 * (len-last_m2), 0.65*last_m2),
                        y=c(0.7 * (len-last_m2), len-0.45*last_m2))


# Create a ggheatmap
p <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson r") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.95,0.4),
        legend.text = element_text(size=10),
        legend.title = element_text(size=15))+
  geom_segment(data=my.lines, aes(x,y,xend=xend, yend=yend), size=3, inherit.aes=F)+
  geom_text(data=my.labels, aes(x, y), label=c('M1 genes', 'M2 genes'), size=10, inherit.aes=F)+
  coord_fixed()

df <- rbind(melt(data.frame("M1/M1" = m1_mat)), melt(data.frame("M2/M2" = m2_mat)), melt(data.frame('M1/M2' = as.vector(m1m2_mat))))
df$variable <- sub("\\.", " vs. ", df$variable)

hh <- ggplot(df,aes(x=value, fill=variable)) + geom_density(alpha=0.35) + scale_fill_brewer(palette = 'Dark2') + labs(x='pearson r', y='density') + guides(fill=guide_legend(title="Genesets"), size=15) + 
  theme(axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15)
        )

vp <- viewport(width = 0.45, height = 0.45, x = 0.95, y = 0.55, just = c("right", "bottom"))

pdf(file = 'fig_3b.pdf', width = 20, height = 20)
print(p)
print(hh, vp=vp)
dev.off()
write.table(upper_tri[colnames(upper_tri), colnames(upper_tri)], file='fig_3b.tsv', sep = '\t', row.names = T, col.names=T, quote=F)
