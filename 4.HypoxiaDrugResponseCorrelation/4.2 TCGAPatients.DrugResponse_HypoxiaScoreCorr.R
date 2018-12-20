##The correlation of imputred drug response and hypoxia in TCGA patients
my.cor.test<- function(...) {
  obj<-try(cor.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(c(NA,NA)) else return(obj[c("estimate","p.value")])
}
TCGA_DS_HS_CorFolder <- "/extraspace/yye1/analysis/Hypoxia/TCGAPatientsDR_HypoxiaCor/"
GDSC_screen_cpd <- read.delim("/extraspace/yye1/analysis/Hypoxia/New/DrugResponse/GDSC_screen_cpd.tab")
DrugResponse <- read.csv("/extraspace/yye1/share_data/DrugData/imputed_drug_response_corrected_GLDC.csv",header=T)
colnames(DrugResponse) <- sapply(colnames(DrugResponse),function(x){ifelse(stringr::str_length(x)>4,substr(x,2,5),x)})
DrugResponse.m <- reshape2::melt(DrugResponse,id.vars="Drug",measure.vars=colnames(DrugResponse)[2:ncol(DrugResponse)])

HypoxiaStr <- read.csv("/extraspace/yye1/analysis/Hypoxia/Stratification.tumor/Stratification.DataAll.csv",header=T)
HypoxiaStr$variable <-substr(HypoxiaStr$SampleID,9,12)
DrugResponse.Hypoxia <- merge(HypoxiaStr,DrugResponse.m,by="variable")
DrugResponse.Hypoxia <- DrugResponse.Hypoxia[which(DrugResponse.Hypoxia$type %in% names(table(unique(DrugResponse.Hypoxia[,1:4])$type)[table(unique(DrugResponse.Hypoxia[,1:4])$type) > 30])),]
HypoxiaScore <- readr::read_rds("TCGA_HypoxiaScore.rds")
HypoxiaScore %>% dplyr::mutate(barcode=substr(barcode,1,16)) %>% dplyr::inner_join(DrugResponse.Hypoxia,by=c("barcode"="SampleID")) -> DrugResponse.Hypoxia

readr::write_rds( DrugResponse.Hypoxia,path=paste(TCGA_DS_HS_CorFolder,"DrugResponse.Hypoxia.rds",sep=""))
DrugResponse.Hypoxia <- readr::read_rds(paste(TCGA_DS_HS_CorFolder,"DrugResponse.Hypoxia.rds",sep=""))

DrugResponse.HypoxiaCorr <- sapply(split(DrugResponse.Hypoxia[,c("HypoxiaScore","value")],list(DrugResponse.Hypoxia$type,DrugResponse.Hypoxia$Drug)),function(x)my.cor.test(as.numeric(x$HypoxiaScore),as.numeric(x$value),method="spearman"))
DrugResponse.HypoxiaCorr <- data.frame(t(DrugResponse.HypoxiaCorr ))
DrugResponse.HypoxiaCorr$cancer_types <- data.frame(do.call(rbind,strsplit(rownames(DrugResponse.HypoxiaCorr),"\\.")))$X1
DrugResponse.HypoxiaCorr$Drug <- gsub("BLCA.|BRCA.|CESC.|ESCA.|GBM.|HNSC.|KIRP.|LGG.|LIHC.|LUAD.|LUSC.|OV.|PAAD.|PRAD.|SARC.|TGCT.|THCA.|THYM.","",row.names(DrugResponse.HypoxiaCorr))
DrugResponse.HypoxiaCorr$FDR <- p.adjust(DrugResponse.HypoxiaCorr$p.value,method="fdr")
DrugResponse.HypoxiaCorr$Sign <- ifelse(DrugResponse.HypoxiaCorr$FDR < 0.05 & DrugResponse.HypoxiaCorr$estimate > 0.2,"Pos",ifelse(DrugResponse.HypoxiaCorr$FDR < 0.05 & DrugResponse.HypoxiaCorr$estimate < -0.2,"Neg",""))
DrugResponse.HypoxiaCorr$FDRlog10 <- -log10(DrugResponse.HypoxiaCorr$FDR+10^-10)
DrugResponse.HypoxiaCorr$estimate <- as.numeric(DrugResponse.HypoxiaCorr$estimate)
DrugResponse.HypoxiaCorr <- DrugResponse.HypoxiaCorr[order(DrugResponse.HypoxiaCorr$estimate),]
DrugResponse.HypoxiaCorr$p.value <- unlist(DrugResponse.HypoxiaCorr$p.value)

DrugResponse.HypoxiaCorrSign <- DrugResponse.HypoxiaCorr[DrugResponse.HypoxiaCorr$Sign !="",]
DrugResponse.HypoxiaCorrSign1 <- merge(DrugResponse.HypoxiaCorrSign,GDSC_screen_cpd,by="Drug")
write.table(DrugResponse.HypoxiaCorrSign1,file=paste(TCGA_DS_HS_CorFolder,"DrugResponse.HypoxiaCorrSign.tab",sep=""),quote = F,row.names = F,sep="\t")

DrugResponse.HypoxiaCorrSignCount <- sapply(split(DrugResponse.HypoxiaCorrSign[,"Sign"],DrugResponse.HypoxiaCorrSign$cancer_types),function(x)c(length(x[x=="Neg"]),length(x[x=="Pos"])))
DrugResponse.HypoxiaCorrSignCount <- data.frame(t(DrugResponse.HypoxiaCorrSignCount))
DrugResponse.HypoxiaCorrSignCount$cancer_types <- rownames(DrugResponse.HypoxiaCorrSignCount)
colnames(DrugResponse.HypoxiaCorrSignCount) <- c("Neg","Pos","cancer_types")
DrugResponse.HypoxiaCorrSignCountm <- reshape2::melt(DrugResponse.HypoxiaCorrSignCount,id.vars="cancer_types",measure.vars=c("Neg","Pos"))

pdf(paste(TCGA_DS_HS_CorFolder,"DrugResponse.HypoxiaCorrCount.pdf",sep=""),width=5,height=3)
ggplot(DrugResponse.HypoxiaCorrSignCountm)+
  geom_bar(aes(x=cancer_types,y=value,fill=variable),color=NA,stat = "identity",width = 0.6)+
  scale_x_discrete(limit=DrugResponse.HypoxiaCorrSignCount[order(DrugResponse.HypoxiaCorrSignCount$Neg,DrugResponse.HypoxiaCorrSignCount$Pos),"cancer_types"])+
  scale_fill_manual(limit=c("Pos","Neg"),values=c("magenta","green"),guide=F)+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_text(color="black",angle=90,size=8),
        axis.ticks=element_line(color="black"),
        axis.line = element_line(color = "black"))
dev.off()


DrugResponse.HypoxiaCorrCount <- sapply(split(DrugResponse.HypoxiaCorrSign[,"Sign"],DrugResponse.HypoxiaCorrSign$Drug), function(x)c(length(x[x=="Neg"]),length(x[x=="Pos"])))
DrugResponse.HypoxiaCorrCount <- data.frame(t(DrugResponse.HypoxiaCorrCount))
colnames(DrugResponse.HypoxiaCorrCount) <- c("Resistance","Sensitivity")
DrugResponse.HypoxiaCorrCount$Drug <- rownames(DrugResponse.HypoxiaCorrCount)
DrugResponse.HypoxiaCorrCountF <- DrugResponse.HypoxiaCorrCount[apply(DrugResponse.HypoxiaCorrCount[,1:2],1,sum)>2,]

pdf(paste(TCGA_DS_HS_CorFolder,"DrugResponse.HypoxiaCorr ordered by alphabet.pdf",sep=""),width=10,height = 4)#,width=4000,height=1500,res=400)
ggplot(DrugResponse.HypoxiaCorrSign1)+
  geom_point(aes(x=Drug,y=cancer_types,size=FDRlog10,color=estimate))+
  scale_color_gradientn(limit=c(-0.8,0.8),colours=colorRampPalette(c("green","white","magenta"),space="rgb")(100),name="Drug sensitivity")+#,
  #  scale_color_gradient2(low = "blue",mid="white", high = "red",midpoint = 0,na.value="white",breaks=seq(-2,2,length.out = 3),labels=c(-2,0,2),name="Difference")+
  scale_size_continuous(limit=c(-log10(0.05),10.1),range = c(0.3, 2.8),breaks=c(-log10(0.05),5,10),labels=c("0.05","1e-5","<1e-10"))+
  scale_y_discrete(limit= DrugResponse.HypoxiaCorrSignCount[rev(order(DrugResponse.HypoxiaCorrSignCount$Neg,DrugResponse.HypoxiaCorrSignCount$Pos)),"cancer_types"])+
  #scale_x_discrete(limit= DrugResponse.HypoxiaCorrCountF[order(DrugResponse.HypoxiaCorrCountF$Resistance-DrugResponse.HypoxiaCorrCountF$Sensitivity),"Drug"])+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid=element_line(colour="grey",linetype="dashed"),
        panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.2),
        axis.title=element_blank(),
        axis.text.y=element_text(size=10,colour = "black"),
        axis.text.x=element_text(size=8,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.direction="horizontal",
        legend.position="bottom",
        legend.key = element_rect(fill="white",colour = "black"))+
  geom_point(data = DrugResponse.HypoxiaCorrSign1[grep("Yes",DrugResponse.HypoxiaCorrSign1$Consistence),],aes(x=Drug,y=cancer_types,size=FDRlog10),shape=1)+
  geom_point(data = DrugResponse.HypoxiaCorrSign1[grep("No",DrugResponse.HypoxiaCorrSign1$Consistence),],aes(x=Drug,y=cancer_types,size=FDRlog10),shape=4)+
  geom_point(data = Report_HypoxiaDrugCor[Report_HypoxiaDrugCor$Consistence=="Re" & Report_HypoxiaDrugCor$Sensitivity=="sensitive",],aes(x=Drug,y=cancer_types),shape=24,color="green")+
  geom_point(data = Report_HypoxiaDrugCor[Report_HypoxiaDrugCor$Consistence=="Re" & Report_HypoxiaDrugCor$Sensitivity=="resistance",],aes(x=Drug,y=cancer_types),shape=24,color="magenta")
dev.off()

