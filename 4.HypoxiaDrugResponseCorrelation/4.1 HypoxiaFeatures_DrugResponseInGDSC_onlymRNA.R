### Hypoxia associated molecular signature reponsed to drug sensitivity
source("~/code/geom_circle.R")

GDSC_Drug_Exp.spearman <- read.delim("/extraspace/yye1/share_data/DrugData/GDSC_Drug_Exp.spearman.csv",header=T)
GDSC_Drug_Exp.spearman.colnames <- read.delim("/extraspace/yye1/share_data/DrugData/GDSC_Drug_Exp.spearman.csv",header=F)[1,]
colnames(GDSC_Drug_Exp.spearman) <- GDSC_Drug_Exp.spearman.colnames[1,]
GDSC_screen_cpd <- readr::read_rds(path="/extraspace/yye1/share_data/DrugData/GDSC_compounds_target_GenePathway.rds")

HypoGDSC_Folder <- "/extraspace/yye1/analysis/Hypoxia/4.2HypoxiaFeatures_DrugResponseInGDSC/"

mRNA.PSM <- read.delim("/extraspace/yye1/analysis/Hypoxia/3.PSM/3.1mRNA/HypoxiaAssociated_mRNA.across.cancer.types.txt",header=T)
mRNA.PSM$feature.sig <- data.frame(do.call(rbind,strsplit(as.character(mRNA.PSM$feature.sig),"\\|")))$X1
mRNA.PSM.filter <- mRNA.PSM[which(mRNA.PSM$feature.sig != "?" & abs(mRNA.PSM$coef.sig) >= 1),]
mRNAChangeCount <- t(sapply(split(mRNA.PSM.filter[,"coef.sig"],mRNA.PSM.filter$feature.sig),function(x)c(length(x[x>0]),length(x[x<0])))) %>% data.frame
colnames(mRNAChangeCount) <- c("mRNA_Up","mRNA_Down")
mRNAChangeCount$GeneSymbol <- rownames(mRNAChangeCount)
mRNAChangeCount$Sum <- apply(mRNAChangeCount[,1:2],1,sum)


hypoxiaSigngene9ct <- mRNAChangeCount[mRNAChangeCount$Sum>=9,"GeneSymbol"]

HYassociatedThan9ctGDSC <- GDSC_Drug_Exp.spearman[which(GDSC_Drug_Exp.spearman$GDSC.exp.GENE_SYMBOLS %in% hypoxiaSigngene9ct),]
HYassociatedThan9ctGDSCspearman <- HYassociatedThan9ctGDSC[,gsub("_FDR","",grep("_FDR",colnames(HYassociatedThan9ctGDSC),value=T))]
HYassociatedThan9ctGDSCspearman[is.na(HYassociatedThan9ctGDSCspearman)] <- 0

HYassociatedThan9ctGDSCspearman[HYassociatedThan9ctGDSCspearman > 0.3] <- 1
HYassociatedThan9ctGDSCspearman[HYassociatedThan9ctGDSCspearman < -0.3] <- -1
HYassociatedThan9ctGDSCspearman[abs(HYassociatedThan9ctGDSCspearman) < 1] <- 0
HYassociatedThan9ctGDSCspearman.p <- HYassociatedThan9ctGDSC[,grep("_FDR",colnames(HYassociatedThan9ctGDSC),value=T)]
HYassociatedThan9ctGDSCspearman.p[HYassociatedThan9ctGDSCspearman.p < 0.05] <- 0
HYassociatedThan9ctGDSCspearman <- HYassociatedThan9ctGDSCspearman+HYassociatedThan9ctGDSCspearman.p
HYassociatedThan9ctGDSCspearman[abs(HYassociatedThan9ctGDSCspearman) != 1] <- 0
HYassociatedThan9ctGDSCspearman$gene <- HYassociatedThan9ctGDSC$GDSC.exp.GENE_SYMBOLS
HYassociatedThan9ctGDSCspearman.m <- reshape2::melt(HYassociatedThan9ctGDSCspearman,id.vars = "gene",measure.vars = colnames(HYassociatedThan9ctGDSCspearman)[1:(ncol(HYassociatedThan9ctGDSCspearman)-1)])
colnames(HYassociatedThan9ctGDSCspearman.m) <- c("gene","drug","correlation")
HYassociatedThan9ctGDSCspearman.m <- HYassociatedThan9ctGDSCspearman.m[which(abs(HYassociatedThan9ctGDSCspearman.m$correlation)>0),]
HYassociatedThan9ctGDSCspearman.m <- HYassociatedThan9ctGDSCspearman.m[HYassociatedThan9ctGDSCspearman.m$gene %in% names(table(HYassociatedThan9ctGDSCspearman.m$gene)[table(HYassociatedThan9ctGDSCspearman.m$gene)>=3]),]
HYassociatedThan9ctGDSCspearman.mm <- merge(HYassociatedThan9ctGDSCspearman.m, GDSC_screen_cpd,by.x="drug",by.y="DRUG.NAME")
HYassociatedThan9ctGDSCspearman.mm[,"TARGET.PATHWAY"][HYassociatedThan9ctGDSCspearman.mm$TARGET.PATHWAY %in% c("chromain  histone acetylation","chromatin  histone methylation","chromatin  other")] <- "chromatin signature"

pathway_loc <- data.frame(pathway=names(table(unique(HYassociatedThan9ctGDSCspearman.mm[,c("drug","TARGET.PATHWAY")])$TARGET.PATHWAY)),
                          radius=table(unique(HYassociatedThan9ctGDSCspearman.mm[,c("drug","TARGET.PATHWAY")])$TARGET.PATHWAY))[c(1,3)]
set.seed(1)
random_radius <- rep(c(1,2,3))
pathway_loc$pathway_x <- c(sin(seq(-1.56,1.4,length.out=9)),2+cos(seq(-3.14,0,length.out=8))) #cos(seq(-3.14,0,length.out=17))#c(sin(seq(-1,1,length.out=9)),cos(seq(-1,1,length.out=8)))
pathway_loc$pathway_y <-  c(0.8*cos(seq(-1.45,1.4,length.out=9)),-0.1+0.8*sin(seq(-3.14,0,length.out=8)))#sin(seq(-3.14,0,length.out=17)) #c(cos(seq(-1,1,length.out=9)),sin(seq(-1,1,length.out=8)))
drug_loc <- merge(unique(HYassociatedThan9ctGDSCspearman.mm[,c("drug","TARGET.PATHWAY")]),pathway_loc,by.x="TARGET.PATHWAY",by.y="pathway")
drugx<- c()
for( i in pathway_loc$pathway){
  sub <- drug_loc[which(drug_loc$TARGET.PATHWAY == i),]
  locx <- unique(sub[,"pathway_x"]) + 0.1 * cos(seq(-3.14,3.14,length.out = nrow(sub)+1)[1:nrow(sub)])
  drugx<- c(drugx,locx)
}
drug_loc$drugx<- drugx

drugy<- c()
for( i in pathway_loc$pathway){
  sub <- drug_loc[which(drug_loc$TARGET.PATHWAY == i),]
  locy <- unique(sub[,"pathway_y"]) + 0.1 * sin(seq(-3.14,3.14,length.out = nrow(sub)+1)[1:nrow(sub)])
  drugy<- c(drugy,locy)
}
drug_loc$drugy<- drugy
drug_num <- table(HYassociatedThan9ctGDSCspearman.mm$drug)[table(HYassociatedThan9ctGDSCspearman.mm$drug)>0]
drug_num <- data.frame(drug_num)
drug_num$drug <- rownames(drug_num)
drug_loc <- merge(drug_loc,drug_num,by="drug")
correlationgenes <- sapply(split(HYassociatedThan9ctGDSCspearman.mm[,"correlation"],HYassociatedThan9ctGDSCspearman.mm$gene),function(x)table(factor(x,levels=c(-1,1))))
correlationgenes <- data.frame(t(correlationgenes))
correlationgenes$gene <- rownames(correlationgenes)
correlationgenes$cor <- correlationgenes$X1 - correlationgenes$X.1
correlationgenes$corenrich <- rep("pos",nrow(correlationgenes))
correlationgenes["corenrich"][correlationgenes["cor"] < 0] <- "neg"
for(i in c("pos","neg")){
  if(i=="pos"){
    sub <- correlationgenes[which(correlationgenes$corenrich=="pos"),]
    sub$genex <-  0.2 * cos(seq(-3.14,3.14,length.out = nrow(sub)))
    sub$geney <- -0.5 +0.18*sin(seq(-3.14,3.14,length.out = nrow(sub)))
    correlationgenes.m <- sub
  }else{
    sub <- correlationgenes[which(correlationgenes$corenrich=="pos"),]
    sub$genex <-  2+0.2 * cos(seq(-3.14,3.14,length.out = nrow(sub)))
    sub$geney <- 0.5 +0.18*sin(seq(-3.14,3.14,length.out = nrow(sub)))
    correlationgenes.m <- rbind(correlationgenes.m,sub)
  }
}

HYassociatedThan9ctGDSCspearman.mmm <- merge(HYassociatedThan9ctGDSCspearman.mm,correlationgenes.m[,c("gene","corenrich","genex","geney")],by="gene")
drug_gene <- merge(HYassociatedThan9ctGDSCspearman.mmm,drug_loc,by="drug")
correlationgenes <- correlationgenes[order(correlationgenes$cor),]
correlationgenes.l <- correlationgenes
correlationgenes.l$genex <- seq(-0.6,2.6,length.out = nrow(correlationgenes.l))
correlationgenes.l$geney <- rep(-1.5,nrow(correlationgenes.l))

HYassociatedThan9ctGDSCspearman.mmm <- merge(HYassociatedThan9ctGDSCspearman.mm,correlationgenes.l[,c("gene","corenrich","genex","geney")],by="gene")
drug_gene <- merge(HYassociatedThan9ctGDSCspearman.mmm,drug_loc,by="drug")

pdf(paste(HypoGDSC_Folder,"GDSC Hypoxia-associated mRNAs and drug correlation at least 9 cancer types.pdf",sep=""),width=10,heigh=4.2)
ggplot(pathway_loc,aes(x=pathway_x,y=pathway_y))+
  geom_circle(radius = 0.04)+
  geom_segment(data=drug_gene[which(drug_gene$correlation==1),],aes(x =genex, y = geney, xend = drugx, yend = drugy), color = "lightpink",alpha=0.2)+
  geom_segment(data=drug_gene[which(drug_gene$correlation==-1),],aes(x =genex, y = geney, xend = drugx, yend = drugy), color = "cyan",alpha=0.2)+
  geom_point(data=drug_loc,aes(x=drugx,y=drugy,size=drug_num),alpha=0.9,color="darkorange")+
  scale_size_continuous(range = c(0.2,2),breaks = c(1,25,50))+
  #scale_x_continuous(limits = c(-1.6,3.2),expand = c(0,0))+
  #scale_y_continuous(limits = c(-1.6,1.6),expand = c(0,0))+
  geom_point(data=correlationgenes.l,aes(x=genex,y=geney),color="darkgreen",size=0.2)+
  geom_text(data=pathway_loc,aes(x=pathway_x,y=pathway_y,label=pathway),size=3)+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks = element_blank())
dev.off()

correlationgenes.lm <- reshape2::melt(correlationgenes.l,id.vars = "gene",measure.vars = c("X.1","X1"))
correlationgenes.lm$value <- -correlationgenes.lm$value
pdf(paste(HypoGDSC_Folder,"GDSC Hypoxia-associated mRNAs and drug correlation pairs at least 9 cancer types.pdf",sep=""),width=10,heigh=2.5)
ggplot(correlationgenes.lm,aes(x=gene,y=value,fill=factor(variable)))+
  geom_bar(color=NA,width = 0.4,stat = "identity")+
  scale_x_discrete(limits=correlationgenes.l$gene)+
  scale_fill_manual(limits=c("X.1","X1"),values = c("cyan","lightpink"),guide=F)+
  ylab("Compounds Count")+
  theme(axis.text.x=element_text(size=5,angle=90,hjust=1,vjust=0.5,color="black"),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_line(colour=NA),
        axis.line.y = element_line(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=10,color="black",angle=90,vjust=0.5,hjust=0.5))#+
dev.off()