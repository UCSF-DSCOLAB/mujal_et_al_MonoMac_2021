#P1 is bulk B16 tumor myeloid cells. See Fig. 1B
P1<-SetAllIdent(P1, id="res.0.6")
#Subset Csf1r+ Mafb+ clusters. See. Fig. 1C
P1Mono<-SubsetData(P1, ident.use=c(0,1,2,3,5))

#Generate cell signature score
#markerfiles refers to gene signatures. See Suppl 1E

cells <- as.data.frame(GetCellEmbeddings(P1Mono, reduction.type="tsne"))
datafile <- P1Mono@data
threshold <- .5
experimentcode="P1Mono"


for (i in strsplit(markerfiles, "_", fixed=T)){
 	
 	i <- unlist(i)[1]
 	names <- c(names,i)}


mydata <- list()
count <- 0
 for (i in markerfiles){
 	if(endsWith(i, ".txt")==TRUE){
 	df <- cells
 	count <- count + 1
 	id <- names[count]
 	file <- read.table(paste("marker_files",i, sep="/"), sep="\t", stringsAsFactors=F, fill=T, header=T)
 	file <- file[order(file$adj.P.Val),][1:100,]
 	genes <- intersect(rownames(datafile), file[which(file$logFC>0 & file$Gene != "NA"), ]$Gene)
 	cat(paste(i,"\n", sep=""))
 	print(genes)
 	cat("\n")
 	
 	totals <-as.matrix(datafile[genes,])
 	#sctotals <- scale(t(totals))
 	
 	sig <- rowMeans(t(totals))
 	#sig <- sig + abs(min(sig))
 	sig <- sig / max (sig)
 	sig <- ifelse(sig > threshold, 1, 0)
 	
 	df$totals <- sig
 	mydata[[id]] <- df}}
      
#Generate overlay visualization on tSNE plot
 seethrough <- rgb(0,0,0, alpha=0.01)
 
 dir.create(paste(experimentcode,"_Threshold_",threshold*100, sep=""))
 
 for (i in c(1:length(mydata))) {
     ggplot(mydata[[i]], aes(x=tSNE_1, y=tSNE_2, color=as.factor(totals))) + geom_point()+scale_color_manual(values=c(seethrough,"red"))+theme(legend.position="none")
     ggsave(paste(experimentcode,"_Threshold_",threshold*100,"/",names[i],"_",threshold*100,"_saturated.pdf", sep=""), height=7, width=7)}
