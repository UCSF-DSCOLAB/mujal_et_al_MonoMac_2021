#P1 is bulk tumor myeloid cell sample. See Fig 1B
#Reclustered P1, Cluster 2. See Fig 2C
P1<-SetAllIdent(P1, id="res.0.6")
P1_Arg1_new<-names(P1@ident[which(P1@ident==2)])
P1_Arg1_new<-CreateSeuratObject(raw.data=P1@raw.data[,P1_Arg1_new],min.cells=0,min.genes=0,project="P1_Arg1_new")
P1_Arg1_new<-NormalizeData(object=P1_Arg1_new,normalization.method="LogNormalize",scale.factor=10000)
P1_Arg1_new<-FindVariableGenes(object=P1_Arg1_new,mean.function=ExpMean,dispersion.function=LogVMR,x.low.cutoff=0.0125,x.high.cutoff=5,y.cutoff=0.5)
P1_Arg1_new<-ScaleData(object=P1_Arg1_new)
P1_Arg1_new<-RunPCA(object=P1_Arg1_new,pc.genes=P1_Arg1_new@var.genes,do.print=TRUE,pcs.print=1:5,genes.print=5,pcs.compute=40)
PCElbowPlot(P1_Arg1_new,num.pc=40)
P1_Arg1_new<-FindClusters(object=P1_Arg1_new, reduction.type="pca",dims.use=1:15, resolution=c(0.2,0.4,0.5,0.6,0.7,0.8),print.output=0,save.SNN=TRUE)
P1_Arg1_new<-RunTSNE(object=P1_Arg1_new,dims.use=1:15,do.fast=TRUE)
