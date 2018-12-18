###Hypoxia effect tumor survival.
#!usr/bin/Rscript
Outpath <- "/extraspace/yye1/analysis/Hypoxia/Hypoxia_survival/"
setwd(Outpath)
library(methods)
library(survival)

FileLocation <- "/extraspace/yye1/analysis/Hypoxia/New/DoubleCheck/Stratification.tumor/"
cancerNames <- c("BLCA","BRCA","CESC","ESCA","GBM","HNSC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PAAD","SKCM","PRAD","SARC","STAD","THCA","TGCT","THYM","UCEC")
folder <- "survival_months"
if (!file.exists(folder)) { dir.create(folder) }
nn1 <- 20

for(cancer in cancerNames) {
  ## read apa result file one by one for each cancer type
  Stratification.Data <- read.delim(paste(FileLocation,cancer,".Hypoxia.stratification.txt",sep=""),header=T)
  Stratification.Data <- Stratification.Data[Stratification.Data$myclusters != "intermediate",]
  Stratification.Data$ID <- substr(Stratification.Data$SampleID,1,12)
  Stratification.Data$ID <- gsub("\\.","\\-",Stratification.Data$ID)

  if (nrow(Stratification.Data) >= 20) {# at lest 20 samples required
    nameOfClinicFile <- paste( "/extraspace/TCGA/TCGA_clinical/",cancer, "_clinical_clean.txt", sep = "")
    if (file.exists(nameOfClinicFile)) {
      clinic_data <- read.table(nameOfClinicFile, sep = "\t", header = TRUE, comment.char = "", quote = "", fill = TRUE)
      clinic_ID <- as.vector(clinic_data[,1])
      commonsample <- intersect(Stratification.Data$ID, clinic_ID) ## get the sample with clinc info
      clinic_data <- clinic_data[match(commonsample, clinic_ID),]
      clinic_data <- merge(clinic_data,Stratification.Data,by.x="barcode",by.y="ID")
      
      clinic_data$os_months <- clinic_data$os_days/30
      
      Result <- matrix(NA, 1, 7, byrow=TRUE)
      colnames(Result) <- c("coef", "Exp(coef)", "Coxp", "KMp", "Hypoxia_L","Hypoxia_H",'CI_95%_for_HR')
      keeplink <- which(!is.na(clinic_data[,"os_status"]) & !is.na(clinic_data[,"os_days"]) & !is.na(clinic_data[,"os_days"])  & clinic_data[,"os_days"] >= 0 )
      Time.dfs <- as.numeric(as.vector(clinic_data[keeplink,"os_months"]))
      cen.status <- ifelse(as.vector(clinic_data[keeplink,"os_status"]) == "Dead", 1,0)
      group_refined <-  clinic_data[keeplink,"myclusters"]
      group_refined <- factor(group_refined,levels=c("normoxic","hypoxic"))
       if(length(group_refined) >= nn1){
       
        test.data1 <- list(time     = Time.dfs,
                           status   = cen.status,
                           group    = as.factor(group_refined))
        tempresult<-try(model1 <- coxph(Surv(time, status) ~ group, data=test.data1, na.action=na.exclude),silent=TRUE)
        if(!is(tempresult, "try-error")){
          Result[1, c("coef", "Exp(coef)", "Coxp")] <- summary(model1)$coefficients[1,c("coef", "exp(coef)", "Pr(>|z|)" )]
          HR_detail = summary(model1)
          Result[1, c('CI_95%_for_HR')] = paste(signif(HR_detail$conf.int[,"lower .95"],2),'-',signif(HR_detail$conf.int[,"upper .95"],2),sep=" ")
          Result[1, c( "Hypoxia_L","Hypoxia_H")] <- table(as.factor(group_refined))
          test.data1 <- list(time     = Time.dfs,
                             status   = cen.status,
                             group    = as.factor(group_refined))
          model1 <- survdiff(Surv(time, status) ~ group, data=test.data1, na.action=na.exclude)
          Result[1, c("KMp")] <- 1-pchisq(model1$chisq, df=length(levels(factor(group_refined)))-1)
          # FDRAdjustPvalue <- p.adjust(Result[1, c("KMp")], method="fdr")
          
          #  }
        }
        
      }######### end of if each marker has at least 20 no NA samples
      
      if(as.numeric(Result[1, "KMp"]) < 0.05){
        fit <- survfit(Surv(time, status) ~ group, data=test.data1, na.action=na.exclude)
        pdf(paste(folder,"/",cancer,".survival_only_oxygen.pdf",sep=""),width=6,height = 6)
        
        plot(fit,col=c("blue","red"),lty=1,lwd=2,mark.time=TRUE,main=paste(cancer,"Kaplan-Meier Curves ",sep=" "),xlab = "Survival in days",cex.lab=1.5,cex.axis=1.2)
        legend("topright", attributes(as.factor(test.data1$group))$levels, col=c("blue","red","gray"),lty=1,box.lwd=2,lwd=2)
        text(max(Time.dfs)/6,0.1,paste("p = ",signif(as.numeric(Result[1, c("KMp")]),digits = 2),sep=""),cex=1.3)
        dev.off()
        #print(summary(fit)$table[, "0.95LCL"])
      }
      if(cancer == "BLCA"){
        Result.All <- Result
      }else{
        Result.All <- rbind(Result.All, Result)
      }
    } ##########  end of survival

  }}
rownames(Result.All) <- cancerNames


Result.All <- data.frame(Result.All)
tmp <- Result.All$Coxp
CoxFDRAdjustPvalue <- p.adjust(tmp, method="fdr")
CoxOtherAdjustPvalue <- p.adjust(tmp, method="bonferroni")
tmp <- Result.All$KMp
KMpFDRAdjustPvalue <- p.adjust(tmp, method="fdr")
KMpOtherAdjustPvalue <- p.adjust(tmp, method="bonferroni")
Result.All <- cbind( cancerNames,Result.All,  FDR = as.numeric(CoxFDRAdjustPvalue), Bonferroni = CoxOtherAdjustPvalue,
                     KMFDR = KMpFDRAdjustPvalue, KMBonferroni = KMpOtherAdjustPvalue)
write.table(Result.All,file="survival.analysis.results.txt",quote = F,row.names = F,sep="\t")



Result.All <- read.delim("survival.analysis.results.txt")
Result.All$Exp.coef. <- as.numeric(Result.All$Exp.coef.)
Result.All1 <- Result.All
Result.All1[,"Exp.coef."][Result.All1[,"Exp.coef."] > 3] <- 3
Result.All1[,"Exp.coef."][Result.All1[,"Exp.coef."] >= 1.5 & Result.All1[,"Exp.coef."] < 3] <- 2
Result.All1[,"Exp.coef."][Result.All1[,"Exp.coef."] > 1 & Result.All1[,"Exp.coef."] < 1.5] <- 1
Result.All1[,"Exp.coef."][Result.All1[,"Exp.coef."] < 0.1] <- -3
Result.All1[,"Exp.coef."][Result.All1[,"Exp.coef."] >= 0.1 & Result.All1[,"Exp.coef."] < 0.15] <- -2

Result.All1[,"Exp.coef."][Result.All1[,"Exp.coef."] > 0.15 & Result.All1[,"Exp.coef."] < 1] <- -1
Result.All1$sizev <- ifelse(Result.All1$FDR < 0.15,1,0)
pdf("Hazard ratio summary.pdf",width=8,height = 2)
ggplot(Result.All1,aes(y= 1,x=cancerNames))+
  geom_point(aes(color=factor(Exp.coef.),size=sizev),shape=15)+
  scale_color_manual(limit=c(-3,-2,-1,1,2,3),values=colorRampPalette(c("blue","white","red"),space="rgb")(20)[c(1,4,8,12,16,20)],labels=c("<0.1","0.1-0.15","0.15-1.0","1.0-1.5","1.5-3",">3"),name="HR")+
  scale_size(limit=c(0,1),range=c(2,5),breaks = c(0,1),labels=c("FDR>1.5","FDR<=0.15"),name="")+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        panel.grid.major=element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(size=12,colour = "black",angle=90,vjust=0.5,hjust=1),
        axis.ticks=element_blank(),legend.text=element_text(size=10),legend.background = element_blank(),
        axis.line=element_blank(),legend.position = "bottom",legend.direction = "horizontal",
        legend.title=element_text(size=12,vjust=0))
dev.off()

