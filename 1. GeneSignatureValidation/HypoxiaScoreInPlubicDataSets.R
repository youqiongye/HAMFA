#Hypoxia-score comparison of cancer cell lines and tumor fragments under hypoxic and normoxic condition in 10 datasets
library(magrittr)
require("biomaRt")
my.t.test <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj[c("p.value")])
}
###Hypoxia feature loading
HypoxiaMarkers <- read.delim("HypoxiaMarkerSignature.txt",header=F)$V1
HypoxiaMarkersL <- list(HypoxiaMarkers=HypoxiaMarkers)
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
RefseqSymbolList <- getBM(attributes=c("hgnc_symbol","refseq_mrna","ensembl_transcript_id"),mart=grch37)
RefseqSymbolList <- RefseqSymbolList[RefseqSymbolList$refseq_mrna !="",]
HypoxiaMarkersRefseq <- RefseqSymbolList[RefseqSymbolList$hgnc_symbol %in% HypoxiaMarkers, ]

#1. GSE18494
#2. GSE55935
#3. GSE3188
#4. GSE70051
#5. GSE75034
#6. GSE3051
#7. GSE33521
#8. GSE79069
#9. GSE73556
#10. GSE30979

##1. GSE18494
##Data loading and processing
GPL9419_GSE18494 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GPL9419_GSE18494",comment.char = "#")
GSE18494 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GSE18494_series_matrix.txt",comment.char = "!")
GSE18494$ID_REF <- gsub("_at","",GSE18494$ID_REF)
GSE18494 <- merge(RefseqSymbolList,GSE18494,by.x="refseq_mrna",by.y="ID_REF")
GSE18494_Samples <- read.csv("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GSE18494_SampleAnn.tab",header = F) %>% tidyr::spread(V2,V3)
GSE18494.m <- as.matrix(GSE18494[grep("GSM",colnames(GSE18494),value=T)])
rownames(GSE18494.m) <- GSE18494$hgnc_symbol
##Hypoxia score calculation
GSE18494_ES <- GSVA::gsva(GSE18494.m,HypoxiaMarkersL)
GSE18494_ES <- data.frame(t(GSE18494_ES ))
GSE18494_ES$GSM_ID <-rownames(GSE18494_ES)
GSE18494_ES <- merge(GSE18494_ES,GSE18494_Samples,by.x="GSM_ID","V1")

GSE18494_ES_Diff <- lapply(split(GSE18494_ES[,c("HypoxiaMarkers","time")],GSE18494_ES$`cell line`),function(x){
  b <- x[x$time != "control",]
  ctrl <- x[x$time == "control","HypoxiaMarkers"]
  Diff <- t(sapply(split(b[,"HypoxiaMarkers"],b$time),function(y){c(Dif = mean(y)-mean(ctrl),unlist(my.t.test(y,ctrl)))})) %>% data.frame()
  Diff$time <- rownames(Diff)
  return(Diff)
}) %>% dplyr::bind_rows() %>% dplyr::mutate(`cell line` = rep(names(split(GSE18494_ES[,c("HypoxiaMarkers","time")],GSE18494_ES$`cell line`)),each = 3))

pdf("GSE18494_ES expression.pdf",height = 3,width = 6)
ggplot(GSE18494_ES, aes(x=time,y=HypoxiaMarkers))+
  geom_boxplot()+
  ggbeeswarm::geom_quasirandom(aes(color=time))+
  facet_wrap(.~`cell line`,scales="free_y")+
  scale_x_discrete(limit=c("control","4h","8h","12h"))+
  scale_color_manual(limit=c("control","4h","8h","12h"),values=c("black","green","blue","red"),name="")+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"),
          axis.text.y = element_text(color="black",size=10),
          axis.title.y = element_text(color="black",size=12),
          axis.ticks.x = element_blank(),legend.position = "bottom",legend.direction = "horizontal",
          axis.line.y = element_line(colour = "black"),strip.background = element_blank())
dev.off()

#2. GSE55935
##Data loading and processing
GPL10558_GSE55935 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GPL10558_GSE55935.txt",comment.char = "#")[,c("ID","ILMN_Gene")]
colnames(GPL10558_GSE55935) <- c("ID","Gene.Symbol")
GSE55935_SampleAnn <- read.table("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GSE55935_SampleAnn.tab",header=T)
GPL10558_GSE55935$Gene.Symbol %>% unique
GSE55935 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GSE55935_series_matrix.txt",comment.char = "!")#[,c("ID",GSE55935_SampleAnn$Sample_title)]
GSE55935 <- merge(GPL10558_GSE55935,GSE55935,by.x="ID",by.y = "ID_REF")
GSE55935 <- GSE55935[,c("ID","Gene.Symbol",GSE55935_SampleAnn$Sample_title)]
GSE55935.m <- as.matrix(GSE55935[grep("GSM",colnames(GSE55935),value=T)])
GSE55935.m <- t(apply(GSE55935.m, 1,function(x){(x-mean(x))/sd(x)}))
rownames(GSE55935.m) <- GSE55935$Gene.Symbol

##Hypoxia score calculation
GSE55935_ES <- GSVA::gsva(GSE55935.m,HypoxiaMarkersL)$es.obs
GSE55935_ES <- data.frame(t(GSE55935_ES ))
GSE55935_ES$GSM_ID <-rownames(GSE55935_ES)
GSE55935_ES <- merge(GSE55935_ES,GSE55935_SampleAnn,by.x="GSM_ID","Sample_title")
GSE55935_ES_p <- unlist(my.t.test(GSE55935_ES[GSE55935_ES$Clusters=="Hypoxic",]$HypoxiaMarkers,GSE55935_ES[GSE55935_ES$Clusters=="Normoxic",]$HypoxiaMarkers))

pdf("GSE55935_ES expression.pdf",height = 3,width = 3)
ggplot(GSE55935_ES, aes(x=Clusters,y=HypoxiaMarkers))+
  geom_boxplot(width=0.7)+
  ggbeeswarm::geom_quasirandom(aes(color=Clusters,fill=Clusters,shape=CellLine),width
                               = 0.2)+
  scale_x_discrete(limit=c("Normoxic","Hypoxic"))+
  scale_shape_manual(limit = unique(GSE55935_ES$CellLine),values = 21:24)+
  #geom_line(aes(group=CellLine))+
  scale_fill_manual(limit=c("Normoxic","Hypoxic"),values=c("black","red"),guide=F)+
  scale_color_manual(limit=c("Normoxic","Hypoxic"),values=c("black","red"),name="")+
  theme( panel.background=element_rect(colour=NA,fill="white"),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         axis.title.x=element_blank(),
         axis.text.x=element_text(vjust=0.5,hjust=1,color="black"),
         axis.text.y = element_text(color="black",size=10),
         axis.title.y = element_text(color="black",size=12),
         axis.ticks.x = element_blank(),legend.position = "bottom",legend.direction =
           "horizontal",
         axis.line.y = element_line(colour = "black"),strip.background = element_blank())
dev.off()


#3. GSE3188
GSE3188_NormoxiaSam <- c("GSM71498","GSM71499","GSM71500")
GSE3188_HypoxiaSam <- c("GSM71501","GSM71502","GSM71503")
GPL96_GSE3188 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GPL96_GSE3188.txt",comment.char = "#")[,c("ID","Gene.Symbol")]
GPL96_GSE3188 <- GPL96_GSE3188[GPL96_GSE3188$Gene.Symbol != "",]
GPL96_GSE3188$Gene.Symbol <- data.frame(do.call(rbind,strsplit(toupper(GPL96_GSE3188$Gene.Symbol)," /// ")))$X1


GSE3188 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GSE3188-GPL96_series_matrix.txt",comment.char = "!")#[,c("ID",GSE3188_SampleAnn$Sample_title)]
GSE3188 <- merge(GPL96_GSE3188,GSE3188,by.x="ID",by.y = "ID_REF")
GSE3188 <- GSE3188[,c("ID","Gene.Symbol",GSE3188_HypoxiaSam,GSE3188_NormoxiaSam)]

GSE3188.m <- as.matrix(GSE3188[grep("GSM",colnames(GSE3188),value=T)])
GSE3188.m <- apply(GSE3188.m,2,function(x)log2(x+1))
GSE3188.m <- t(apply(GSE3188.m, 1,function(x){(x-mean(x))/sd(x)}))
rownames(GSE3188.m) <- GSE3188$Gene.Symbol
GSE3188_ES <- GSVA::gsva(GSE3188.m,HypoxiaMarkersL)$es.obs

GSE3188_ES <- data.frame(t(GSE3188_ES ))
GSE3188_ES$GSM_ID <-rownames(GSE3188_ES)
GSE3188_ES$Clusters <-  ifelse(GSE3188_ES$GSM_ID %in% GSE3188_HypoxiaSam,"Hypoxic","Normoxic")


GSE3188_GPL2507_NormoxiaSam <- c("GSM71605","GSM71607")
GSE3188_GPL2507_HypoxiaSam <- c("GSM71609","GSM71611")

GPL2507_GSE3188 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GPL2507_GSE3188.txt",comment.char = "#")[,c("ID","Gene.symbol")]
colnames(GPL2507_GSE3188) <- c("ID","Gene.Symbol")
GPL2507_GSE3188 <- GPL2507_GSE3188[GPL2507_GSE3188$Gene.Symbol != "",]
#GPL2507_GSE3188$Gene.Symbol <- data.frame(do.call(rbind,strsplit(toupper(GPL2507_GSE3188$Gene.Symbol)," /// ")))$X1

GSE3188_GPL2507 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GSE3188-GPL2507_series_matrix.txt",comment.char = "!")#[,c("ID",GSE3188_GPL2507_SampleAnn$Sample_title)]
GSE3188_GPL2507 <- merge(GPL2507_GSE3188,GSE3188_GPL2507,by.x="ID",by.y = "ID_REF")
GSE3188_GPL2507 <- GSE3188_GPL2507[,c("ID","Gene.Symbol",GSE3188_GPL2507_NormoxiaSam,GSE3188_GPL2507_HypoxiaSam)]

GSE3188_GPL2507.m <- as.matrix(GSE3188_GPL2507[grep("GSM",colnames(GSE3188_GPL2507),value=T)])
#GSE3188_GPL2507.m <- apply(GSE3188_GPL2507.m,2,function(x)log2(x+1))
GSE3188_GPL2507.m <- t(apply(GSE3188_GPL2507.m, 1,function(x){(x-mean(x))/sd(x)}))
rownames(GSE3188_GPL2507.m) <- GSE3188_GPL2507$Gene.Symbol
GSE3188_GPL2507_ES <- GSVA::gsva(GSE3188_GPL2507.m,HypoxiaMarkersL)$es.obs

GSE3188_GPL2507_ES <- data.frame(t(GSE3188_GPL2507_ES ))
GSE3188_GPL2507_ES$GSM_ID <-rownames(GSE3188_GPL2507_ES)
GSE3188_GPL2507_ES$Clusters <-  ifelse(GSE3188_GPL2507_ES$GSM_ID %in% GSE3188_GPL2507_HypoxiaSam,"Hypoxic","Normoxic")

GSE3188_ES <- rbind(GSE3188_ES,GSE3188_GPL2507_ES)
GSE3188_ES_p <- unlist(my.t.test(GSE3188_ES[GSE3188_ES$Clusters=="Hypoxic",]$HypoxiaMarkers,GSE3188_ES[GSE3188_ES$Clusters=="Normoxic",]$HypoxiaMarkers))

pdf("GSE3188_ES expression.pdf",height = 3,width = 3)
ggplot(GSE3188_ES, aes(x=Clusters,y=HypoxiaMarkers))+
  geom_boxplot(width=0.7)+
  ggbeeswarm::geom_quasirandom(aes(color=Clusters),width = 0.2)+
  scale_x_discrete(limit=c("Normoxic","Hypoxic"))+
  scale_color_manual(limit=c("Normoxic","Hypoxic"),values=c("black","red"),name="")+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_text(vjust=0.5,hjust=1,color="black"),
          axis.text.y = element_text(color="black",size=10),
          axis.title.y = element_text(color="black",size=12),
          axis.ticks.x = element_blank(),legend.position = "bottom",legend.direction = "horizontal",
          axis.line.y = element_line(colour = "black"),strip.background = element_blank())
dev.off()

#4. GSE70051
GSE70051_SampleAnn <- read.csv("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GSE70051_SampleAnn.tab",header=F)
GSE70051_SampleAnn <- GSE70051_SampleAnn %>% tidyr::spread(V2,V3)
GPL570_GSE70051 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GPL570_GSE9649.txt",comment.char = "#")[,c("ID","Gene.Symbol")]
GPL570_GSE70051 <- GPL570_GSE70051[GPL570_GSE70051$Gene.Symbol != "",]
GPL570_GSE70051$Gene.Symbol <- data.frame(do.call(rbind,strsplit(toupper(GPL570_GSE70051$Gene.Symbol)," /// ")))$X1
GSE70051 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GSE70051_series_matrix.txt",comment.char = "!")#[,c("ID",GSE70051_SampleAnn$Sample_title)]
GSE70051 <- merge(GPL570_GSE70051,GSE70051,by.x="ID",by.y = "ID_REF")
GSE70051 <- GSE70051[,c("ID","Gene.Symbol",GSE70051_SampleAnn$V1)]

GSE70051.m <- as.matrix(GSE70051[grep("GSM",colnames(GSE70051),value=T)])
GSE70051.m <- t(apply(GSE70051.m, 1,function(x){(x-mean(x))/sd(x)}))
rownames(GSE70051.m) <- GSE70051$Gene.Symbol
GSE70051_ES <- GSVA::gsva(GSE70051.m,HypoxiaMarkersL)$es.obs
GSE70051_ES <- data.frame(t(GSE70051_ES ))
GSE70051_ES$GSM_ID <-rownames(GSE70051_ES)
GSE70051_ES <- merge(GSE70051_ES,GSE70051_SampleAnn,by.x="GSM_ID","V1")
GSE70051_ES <- GSE70051_ES[GSE70051_ES$`oxygen level` %in% c("21% Oxygen","1% Oxygen","0.1% Oxygen") & GSE70051_ES$ph == "pH 7.5",]#& GSE70051_ES$`cell line` %in% c("SiHa","FaDuDD"),]

GSE70051_ES_p_0.1O2 <- unlist(my.t.test(GSE70051_ES[GSE70051_ES$`oxygen level`=="0.1% Oxygen",]$HypoxiaMarkers,GSE70051_ES[GSE70051_ES$`oxygen level`=="21% Oxygen",]$HypoxiaMarkers))

GSE70051_ES_p_1O2 <- unlist(my.t.test(GSE70051_ES[GSE70051_ES$`oxygen level`=="1% Oxygen",]$HypoxiaMarkers,GSE70051_ES[GSE70051_ES$`oxygen level`=="21% Oxygen",]$HypoxiaMarkers))

D_GSE70051 <- data.frame(CellLine = paste(unique(GSE70051_ES$`cell line`),collapse = ","),
                         HypoxiaScoreDiff = mean(GSE70051_ES[GSE70051_ES$`oxygen level`=="0.1% Oxygen",]$HypoxiaMarkers)-mean(GSE70051_ES[GSE70051_ES$`oxygen level`=="21% Oxygen",]$HypoxiaMarkers),
                         p.value = GSE70051_ES_p_0.1O2,DataSet = "GSE70051",Time="24h")

pdf("GSE70051_ES expression.pdf",height = 3,width = 3)
ggplot(GSE70051_ES, aes(x=`oxygen level`,y=HypoxiaMarkers))+
  geom_boxplot()+
  ggbeeswarm::geom_quasirandom(aes(color=`oxygen level`))+
  #  facet_wrap(.~`cell line`,scales="free_y")+
  scale_x_discrete(limit=c("21% Oxygen","0.1% Oxygen"))+
  scale_color_manual(limit=c("21% Oxygen","0.1% Oxygen"),values=c("black","red","red","red"),name="")+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_text(vjust=0.5,hjust=1,color="black"),
          axis.text.y = element_text(color="black",size=10),
          axis.title.y = element_text(color="black",size=12),
          axis.ticks.x = element_blank(),legend.position = "bottom",legend.direction = "horizontal",
          axis.line.y = element_line(colour = "black"),strip.background = element_blank())
dev.off()

#5. GSE75034
GPL10558_GSE75034 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GPL10558_GSE55935.txt",comment.char = "#")[,c("ID","ILMN_Gene")]
colnames(GPL10558_GSE75034) <- c("ID","Gene.Symbol")

GPL10558_GSE75034 <- GPL10558_GSE75034[GPL10558_GSE75034$Gene.Symbol != "",]
GSE75034_SampleAnn <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GSE75034_SampleAnn.tab",header=T)

GSE75034 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GSE75034-GPL10558_series_matrix.txt",comment.char = "!")#[,c("ID",GSE75034_SampleAnn$Sample_title)]
GSE75034 <- merge(GPL10558_GSE75034,GSE75034,by.x="ID",by.y="ID_REF")
GSE75034 <- GSE75034[,c("ID","Gene.Symbol",intersect(GSE75034_SampleAnn$Sample_title,colnames(GSE75034)))]

GSE75034.m <- as.matrix(GSE75034[grep("GSM",colnames(GSE75034),value=T)])
GSE75034.m <- apply(GSE75034.m,2,function(x){log2(x+1)})
GSE75034.m <- t(apply(GSE75034.m, 1,function(x){(x-mean(x))/sd(x)}))
rownames(GSE75034.m) <- GSE75034$Gene.Symbol

GSE75034_ES <- GSVA::gsva(GSE75034.m,HypoxiaMarkersL)$es.obs
GSE75034_ES <- data.frame(t(GSE75034_ES ))
GSE75034_ES$GSM_ID <-rownames(GSE75034_ES)
GSE75034_ES <- merge(GSE75034_ES,GSE75034_SampleAnn,by.x="GSM_ID","Sample_title")
GSE75034_ES_p <- unlist(my.t.test(GSE75034_ES[GSE75034_ES$Clusters=="Hypoxia",]$HypoxiaMarkers,GSE75034_ES[GSE75034_ES$Clusters=="Normoxia",]$HypoxiaMarkers,paired=T))
F_GSE75034 <- data.frame(CellLine = paste(unique(GSE75034_ES$CellLine),collapse = ","),
                         HypoxiaScoreDiff = mean(GSE75034_ES[GSE75034_ES$Clusters=="Hypoxia",]$HypoxiaMarkers)-mean(GSE75034_ES[GSE75034_ES$Clusters=="Normoxia",]$HypoxiaMarkers),
                         p.value=GSE75034_ES_p,DataSet = "GSE75034",Time="24h")

pdf("GSE75034_ES expression.pdf",height = 3,width = 3)
ggplot(GSE75034_ES, aes(x=Clusters,y=HypoxiaMarkers))+
  geom_boxplot(width=0.7)+
  ggbeeswarm::geom_quasirandom(aes(color=Clusters,fill=Clusters,shape=CellLine),width= 0.2)+
  scale_x_discrete(limit=c("Normoxia","Hypoxia"))+
  scale_shape_manual(limit = unique(GSE75034_ES$CellLine),values = c(21:25,8))+
  geom_line(aes(group=CellLine),linetype="dashed")+
  scale_fill_manual(limit=c("Normoxia","Hypoxia"),values=c("black","red"),guide=F)+
  scale_color_manual(limit=c("Normoxia","Hypoxia"),values=c("black","red"),guide=F)+
  theme( panel.background=element_rect(colour=NA,fill="white"),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y = element_text(color="black",size=10),
         axis.title.y = element_text(color="black",size=12),
         axis.ticks.x = element_blank(),legend.position = "bottom",legend.direction =
           "horizontal",
         axis.line.y = element_line(colour = "black"),strip.background = element_blank())
dev.off()

#6.GSE3051	HeLa	Cervical cancer

GPL570_GSE3051<- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GPL570_GSE9649.txt",comment.char = "#")[,c("ID","Gene.Symbol")]
GPL570_GSE3051<- GPL570_GSE3051[GPL570_GSE3051$Gene.Symbol != "",]
GPL570_GSE3051$Gene.Symbol <- data.frame(do.call(rbind,strsplit(toupper(GPL570_GSE3051$Gene.Symbol)," /// ")))$X1
GSE3051<- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GSE3051_series_matrix.txt",comment.char = "!")#[,c("ID",GSE3051_SampleAnn$Sample_title)]
GSE3051 <- merge(GPL570_GSE3051,GSE3051,by.x="ID",by.y = "ID_REF")

GSE3051.m <- as.matrix(GSE3051[grep("GSM",colnames(GSE3051),value=T)])
GSE3051.m <- t(apply(GSE3051.m, 1,function(x){(x-mean(x))/sd(x)}))
rownames(GSE3051.m) <- GSE3051$Gene.Symbol
GSE3051_ES <- GSVA::gsva(GSE3051.m,HypoxiaMarkersL)$es.obs

GSE3051_ES <- data.frame(t(GSE3051_ES ))
GSE3051_ES$GSM_ID <-rownames(GSE3051_ES)
GSE3051_ES$Clusters <- rep(c("Hypoxic","Normoxic"),each=3)
GSE3051_ES_p <- unlist(my.t.test(GSE3051_ES[GSE3051_ES$Clusters=="Hypoxic",]$HypoxiaMarkers,GSE3051_ES[GSE3051_ES$Clusters=="Normoxic",]$HypoxiaMarkers,paired=T))

I_GSE3051 <- data.frame(CellLine = "Hela",HypoxiaScoreDiff = mean(GSE3051_ES[GSE3051_ES$Clusters=="Hypoxic",]$HypoxiaMarkers)-mean(GSE3051_ES[GSE3051_ES$Clusters=="Normoxic",]$HypoxiaMarkers),
                        p.value=GSE3051_ES_p,DataSet="GSE3051",Time = "24h" )

pdf("GSE3051_ES expression.pdf",height = 3,width = 3)
ggplot(GSE3051_ES, aes(x=Clusters,y=HypoxiaMarkers))+
  geom_boxplot(width=0.7)+
  ggbeeswarm::geom_quasirandom(aes(color=Clusters),width = 0.2)+
  scale_x_discrete(limit=c("Normoxic","Hypoxic"))+
  scale_color_manual(limit=c("Normoxic","Hypoxic"),values=c("black","red"),name="")+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_text(vjust=0.5,hjust=1,color="black"),
          axis.text.y = element_text(color="black",size=10),
          axis.title.y = element_text(color="black",size=12),
          axis.ticks.x = element_blank(),legend.position = "bottom",legend.direction = "horizontal",
          axis.line.y = element_line(colour = "black"),strip.background = element_blank())+
  annotate("text",x=1,y=0,label=paste("p = ",signif(GSE3051_ES_p,digits = 2),sep=""))
dev.off()

##7	GSE33521	HeLa	Cervical cancer
GSE33521_SampleAnn <- data.frame(GSM_ID = c("GSM829385","GSM829386","GSM829387","GSM829388"),clusters=c("Normoxia","Normoxia","Hypoxia","Hypoxia"))

GPL6480_GSE33521 <- GPL6480_GSE31079 
GSE33521 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GSE33521_series_matrix.txt",comment.char = "!")#[,c("ID",GSE33521_SampleAnn$Sample_title)]
GSE33521 <- merge(GPL6480_GSE33521,GSE33521,by.x="ID",by.y = "ID_REF")
GSE33521 <- GSE33521[,c("ID","Gene.Symbol",GSE33521_SampleAnn$GSM_ID)]

GSE33521.m <- as.matrix(GSE33521[grep("GSM",colnames(GSE33521),value=T)])
GSE33521.m <- t(apply(GSE33521.m, 1,function(x){(x-mean(x))/sd(x)}))
rownames(GSE33521.m) <- GSE33521$Gene.Symbol
GSE33521_ES <- GSVA::gsva(GSE33521.m,HypoxiaMarkersL)$es.obs

GSE33521_ES <- data.frame(t(GSE33521_ES ))
GSE33521_ES$GSM_ID <-rownames(GSE33521_ES)

GSE33521_ES <- merge(GSE33521_ES,GSE33521_SampleAnn,by="GSM_ID")
GSE33521_ES_p <- my.t.test(c(-0.506,-0.202),c(0.008,0.274),paired=T)#unlist(my.t.test(GSE33521_ES[GSE33521_ES$clusters=="Hypoxia",]$HypoxiaMarkers,GSE33521_ES[GSE33521_ES$clusters=="Normoxia",]$HypoxiaMarkers))
J_GSE33521 <- data.frame(CellLine = "Hela",HypoxiaScoreDiff = mean(GSE33521_ES[GSE33521_ES$clusters=="Hypoxia",]$HypoxiaMarkers) - mean(GSE33521_ES[GSE33521_ES$clusters=="Normoxia",]$HypoxiaMarkers),
                         p.value=GSE33521_ES_p,DataSet= "GSE33521",Time="6h")

pdf("GSE33521_ES expression.pdf",height = 3,width = 3)
ggplot(GSE33521_ES, aes(x=clusters,y=HypoxiaMarkers))+
  geom_boxplot()+
  ggbeeswarm::geom_quasirandom(aes(color=clusters))+
  #  facet_wrap(.~`cell line`,scales="free_y")+
  scale_x_discrete(limit=c("Normoxia","Hypoxia"))+
  scale_color_manual(limit=c("Normoxia","Hypoxia"),values=c("black","red","red","red"),name="")+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_text(vjust=0.5,hjust=1,color="black"),
          axis.text.y = element_text(color="black",size=10),
          axis.title.y = element_text(color="black",size=12),
          axis.ticks.x = element_blank(),legend.position = "bottom",legend.direction = "horizontal",
          axis.line.y = element_line(colour = "black"),strip.background = element_blank())+
  annotate("text",x=1,y=0,label=paste("p = ",signif(GSE33521_ES_p,digits = 2),sep=""))
dev.off()
#8. GSE79069	HeLa	Cervical cancer
GPL10558_GSE79069 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GPL10558_GSE55935.txt",comment.char = "#")[,c("ID","ILMN_Gene")]
colnames(GPL10558_GSE79069) <- c("ID","Gene.Symbol")
GSE79069_SampleAnn <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GSE79069_SampleAnn.tab")

GSE79069 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GSE79069_series_matrix.txt",comment.char = "!")#[,c("ID",GSE79069_SampleAnn$Sample_title)]
GSE79069 <- merge(GPL10558_GSE79069,GSE79069,by.x="ID",by.y = "ID_REF")
GSE79069 <- GSE79069[,c("ID","Gene.Symbol",GSE79069_SampleAnn$GSM_ID)]

GSE79069.m <- as.matrix(GSE79069[grep("GSM",colnames(GSE79069),value=T)])
GSE79069.m <- apply(GSE79069.m ,2,function(x)log2(x+1))
GSE79069.m <- t(apply(GSE79069.m, 1,function(x){(x-mean(x))/sd(x)}))
rownames(GSE79069.m) <- GSE79069$Gene.Symbol
GSE79069_ES <- GSVA::gsva(GSE79069.m,HypoxiaMarkersL)$es.obs

GSE79069_ES <- data.frame(t(GSE79069_ES ))
GSE79069_ES$GSM_ID <-rownames(GSE79069_ES)

GSE79069_ES <- merge(GSE79069_ES,GSE79069_SampleAnn,by="GSM_ID")
GSE79069_ES <- GSE79069_ES[GSE79069_ES$Time=="20h",]

GSE79069_ES_p <-unlist(my.t.test(GSE79069_ES[GSE79069_ES$Clusters=="Hypoxia",]$HypoxiaMarkers,GSE79069_ES[GSE79069_ES$Clusters=="Normoxia",]$HypoxiaMarkers))
K_GSE79069 <- data.frame(CellLine = "Hela",HypoxiaScoreDiff = mean(GSE79069_ES[GSE79069_ES$Clusters=="Hypoxia",]$HypoxiaMarkers) - mean(GSE79069_ES[GSE79069_ES$Clusters=="Normoxia",]$HypoxiaMarkers),
                         p.value = GSE79069_ES_p,DataSet = "GSE79069",Time="20h")


pdf("GSE79069_ES expression.pdf",height = 3,width = 3)
ggplot(GSE79069_ES, aes(x=Clusters,y=HypoxiaMarkers))+
  geom_boxplot()+
  ggbeeswarm::geom_quasirandom(aes(color=Clusters))+
  #  facet_wrap(.~`cell line`,scales="free_y")+
  scale_x_discrete(limit=c("Normoxia","Hypoxia"))+
  scale_color_manual(limit=c("Normoxia","Hypoxia"),values=c("black","red","red","red"),name="")+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_text(vjust=0.5,hjust=1,color="black"),
          axis.text.y = element_text(color="black",size=10),
          axis.title.y = element_text(color="black",size=12),
          axis.ticks.x = element_blank(),legend.position = "bottom",legend.direction = "horizontal",
          axis.line.y = element_line(colour = "black"),strip.background = element_blank())+
  annotate("text",x=1,y=0,label=paste("p = ",signif(GSE79069_ES_p,digits = 2),sep=""))
dev.off()

#9. GSE73556	D456MG GICs 	 Glioblastoma 
GSE73556_SampleAnn <- data.frame(GSM_ID = c("GSM1897949","GSM1897950","GSM1897951","GSM1897946","GSM1897947","GSM1897948"),clusters=c("Normoxia","Normoxia","Normoxia","Hypoxia","Hypoxia","Hypoxia"))

GPL4133_GSE73556 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GPL4133_GSE73556.txt",comment.char = "#")[,c("ID","GENE_SYMBOL")]
colnames(GPL4133_GSE73556) <- c("ID","Gene.Symbol")
GPL4133_GSE73556 <- GPL4133_GSE73556[GPL4133_GSE73556$Gene.Symbol != "",]

GSE73556 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GSE73556_series_matrix.txt",comment.char = "!")#[,c("ID",GSE73556_SampleAnn$Sample_title)]
GSE73556 <- merge(GPL4133_GSE73556,GSE73556,by.x="ID",by.y = "ID_REF")
GSE73556 <- GSE73556[,c("ID","Gene.Symbol",GSE73556_SampleAnn$GSM_ID)]

GSE73556.m <- as.matrix(GSE73556[grep("GSM",colnames(GSE73556),value=T)])
GSE73556.m <- apply(GSE73556.m,2,function(x){log2(x+1)})
GSE73556.m <- t(apply(GSE73556.m, 1,function(x){(x-mean(x))/sd(x)}))
rownames(GSE73556.m) <- GSE73556$Gene.Symbol
GSE73556_ES <- GSVA::gsva(GSE73556.m,HypoxiaMarkersL)$es.obs

GSE73556_ES <- data.frame(t(GSE73556_ES ))
GSE73556_ES$GSM_ID <-rownames(GSE73556_ES)

GSE73556_ES <- merge(GSE73556_ES,GSE73556_SampleAnn,by="GSM_ID")
GSE73556_ES_p <-unlist(my.t.test(GSE73556_ES[GSE73556_ES$clusters=="Hypoxia",]$HypoxiaMarkers,GSE73556_ES[GSE73556_ES$clusters=="Normoxia",]$HypoxiaMarkers))
L_GSE73556 <- data.frame(CellLine="D456MG GICs",HypoxiaScoreDiff = mean(GSE73556_ES[GSE73556_ES$clusters=="Hypoxia",]$HypoxiaMarkers)- mean(GSE73556_ES[GSE73556_ES$clusters=="Normoxia",]$HypoxiaMarkers),
                         p.value = GSE73556_ES_p,DataSet ="GSE73556",Time = "12h")
pdf("GSE73556_ES expression.pdf",height = 3,width = 3)
ggplot(GSE73556_ES, aes(x=clusters,y=HypoxiaMarkers))+
  geom_boxplot()+
  ggbeeswarm::geom_quasirandom(aes(color=clusters))+
  #  facet_wrap(.~`cell line`,scales="free_y")+
  scale_x_discrete(limit=c("Normoxia","Hypoxia"))+
  scale_color_manual(limit=c("Normoxia","Hypoxia"),values=c("black","red","red","red"),name="")+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_text(vjust=0.5,hjust=1,color="black"),
          axis.text.y = element_text(color="black",size=10),
          axis.title.y = element_text(color="black",size=12),
          axis.ticks.x = element_blank(),legend.position = "bottom",legend.direction = "horizontal",
          axis.line.y = element_line(colour = "black"),strip.background = element_blank())+
  annotate("text",x=1,y=0,label=paste("p = ",signif(GSE73556_ES_p,digits = 2),sep=""))
dev.off()

#10. GSE30979
GSE30979_SampleAnn <- read.table("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GSE30979_SampleAnn.tab",header=T) %>% t() %>% data.frame
GSE30979_SampleAnn$cluster <- data.frame(do.call(rbind,strsplit(rownames(GSE30979_SampleAnn),"\\.")))$X1      #ifelse(substr(rownames(GSE30979_SampleAnn),1,6)=="Normal","Normoxic","Hypoxia")
GSE30979_SampleAnn <- GSE30979_SampleAnn[2:nrow(GSE30979_SampleAnn),]
colnames(GSE30979_SampleAnn) <- c("Sample_title","Clusters")

GPL6244_GSE30979 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GPL6244_GSE30979.txt",comment.char = "#")
GPL6244_GSE30979 <- GPL6244_GSE30979[GPL6244_GSE30979$gene_assignment != "---",]
GPL6244_GSE30979$hgnc_symbol <- data.frame(do.call(rbind,strsplit(GPL6244_GSE30979$gene_assignment," // ")))$X2
GPL6244_GSE30979 <- GPL6244_GSE30979[,c("ID","hgnc_symbol")]
GSE30979 <- read.delim("/extraspace/yye1/analysis/data/GEO_data/Hypoxia/GSE30979_series_matrix.txt",comment.char = "!")
GSE30979 <- merge(GPL6244_GSE30979,GSE30979,by.x="ID",by.y = "ID_REF")

GSE30979.m <- as.matrix(GSE30979[grep("GSM",colnames(GSE30979),value=T)])
GSE30979.m <- t(apply(GSE30979.m, 1,function(x){(x-mean(x))/sd(x)}))
rownames(GSE30979.m) <- GSE30979$hgnc_symbol
GSE30979_ES <- GSVA::gsva(GSE30979.m,HypoxiaMarkersL)$es.obs

GSE30979_ES <- data.frame(t(GSE30979_ES ))
GSE30979_ES$GSM_ID <-rownames(GSE30979_ES)

GSE30979_ES <- merge(GSE30979_ES,GSE30979_SampleAnn,by.x="GSM_ID","Sample_title")
GSE30979_ES_p <- unlist(my.t.test(GSE30979_ES[GSE30979_ES$Clusters=="Hypoxic",]$HypoxiaMarkers,GSE30979_ES[GSE30979_ES$Clusters=="Normoxic",]$HypoxiaMarkers))

pdf("GSE30979_ES expression.pdf",height = 3,width = 3)
ggplot(GSE30979_ES, aes(x=Clusters,y=HypoxiaMarkers))+
  geom_boxplot(width=0.7)+
  ggbeeswarm::geom_quasirandom(aes(color=Clusters),width = 0.2)+
  scale_x_discrete(limit=c("Normoxic","Hypoxic"))+
  scale_color_manual(limit=c("Normoxic","Hypoxic"),values=c("black","red"),name="")+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_text(vjust=0.5,hjust=1,color="black"),
          axis.text.y = element_text(color="black",size=10),
          axis.title.y = element_text(color="black",size=12),
          axis.ticks.x = element_blank(),legend.position = "bottom",legend.direction = "horizontal",
          axis.line.y = element_line(colour = "black"),strip.background = element_blank())
dev.off()


