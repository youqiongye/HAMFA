##1. Stratificating tumors based on hypoxic markers


##1. Stratificating tumors based on hypoxic markers.
Outpath <- "/extraspace/yye1/analysis/Hypoxia/New/DoubleCheck/Stratification.tumor/"
setwd(Outpath)
folder <- "Stratification.tumor"
if (!file.exists(folder)) { dir.create(folder) }
HypoxiaGenes <- read.delim("/extraspace/yye1/analysis/Hypoxia/Hypoxia_makers.txt",header=F)
###mRNA expression data from TCGA mRNA expression
exp.files.names <- list.files(path="/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",pattern="_20160513")
exp.files.Abs <- gsub("_mRNA_each_exp_20160513","",exp.files.names)

HypoxiaGenes.exp <- data.frame(HypoxiaGenes)
colnames(HypoxiaGenes.exp) <- "gene"
for(m in 1:length(exp.files.names)){
  sub.exp <- read.delim(paste("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",exp.files.names[m],sep=""),header=T)
  sub.exp$gene <- data.frame(do.call(rbind, strsplit(as.character(sub.exp$gene),'\\|')))$X1
  sub.exp <- sub.exp[which(sub.exp$gene %in% HypoxiaGenes$V1),] ##get target genes
  if(exp.files.Abs[m]=="SKCM"){
    sub.tumor.names <- colnames(sub.exp)[as.numeric(substr(colnames(sub.exp),14,15)) ==6]
  }else{
    sub.tumor.names <- colnames(sub.exp)[as.numeric(substr(colnames(sub.exp),14,15)) ==1]
  }
  sub.tumor.names <- sub.tumor.names[!is.na(sub.tumor.names)]
  print(length(sub.tumor.names))
  if(length(sub.tumor.names) >= 100){
    sub.exp <- sub.exp[,c("gene",sub.tumor.names)]
    sub.exp[sub.tumor.names] <- log10(sub.exp[sub.tumor.names]+1)
    sub_matrix <- as.matrix(sub.exp[sub.tumor.names])
    sub_matrix <- apply(sub_matrix,1,scale)
    colnames(sub_matrix) <- sub.exp$gene
    rownames(sub_matrix) <- sub.tumor.names
    distancem <- dist(sub_matrix)
    hclust_completem <- hclust(distancem, method = "ward.D")
    dendcompletem <- as.dendrogram(hclust_completem)
    rgb.palette <- colorRampPalette(c("blue","white","red"),space = "rgb")
    mycl <- cutree(hclust_completem,k=3)
    myclusters <- mycl[hclust_completem$order]
    SumZscore <- c()
    for(i in 1:3){
      a <- sum(apply(sub_matrix[which(row.names(sub_matrix) %in% names(myclusters)[myclusters == i]),],2,mean))
      SumZscore <- c(SumZscore,a)
    }
    index_hypoxic <- grep(max(SumZscore),SumZscore)
    index_normoxic <- grep(min(SumZscore),SumZscore)
    index_intermediate <- grep(SumZscore[-c(index_hypoxic,index_normoxic)],SumZscore)
    Stratification <- data.frame(myclusters)
    Stratification$SampleID <- rownames(Stratification)
    Stratification["myclusters"][Stratification["myclusters"]==index_hypoxic] <- "hypoxic"
    Stratification["myclusters"][Stratification["myclusters"]==index_normoxic] <- "normoxic"
    Stratification["myclusters"][Stratification["myclusters"]==index_intermediate] <- "intermediate"
    write.table(Stratification,file=paste("Stratification.tumor/",exp.files.Abs[m],".Hypoxia.stratification.txt",sep=""),quote = F,row.names = F,sep="\t")
    cols <- brewer.pal(max(myclusters), "Set1")
    pdf(paste(Outpath,exp.files.Abs[m],".Hypoxia.stratification.pdf",sep=""),width=14,height = 10)#,width=4000,height=1500,res=400)
    
    heatmap(t(sub_matrix), Colv=dendcompletem, Rowv = NA, scale="row",col=rgb.palette(100),
            breaks=c(seq(min(sub_matrix),0,length.out = 51),seq(0,max(sub_matrix),length.out = 50)),labCol = NA,cexRow = 2)
    fields::image.plot( legend.only=TRUE, col=rgb.palette(100),zlim=c(min(sub_matrix),max(sub_matrix)),legend.shrink = 0.5,horizontal = TRUE,smallplot = c(0.1,0.18, 0.9,0.95), legend.width = 0.2,legend.cex = 1.5,
                axis.args=list(at=pretty(c(min(sub_matrix),max(sub_matrix)),n=5),labels=pretty(c(min(sub_matrix),max(sub_matrix)),n=5),cex=1.7,lwd=0.1),legend.lab = "Row Z-score",legend.position="top")
    dev.off()
  }
}

