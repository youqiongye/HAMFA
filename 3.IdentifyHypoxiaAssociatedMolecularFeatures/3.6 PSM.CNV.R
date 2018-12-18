###Identify hypoxia associated DNA copy number
Outpath <- "/extraspace/yye1/analysis/Hypoxia/New/DoubleCheck/HypoxiaAssociatedFeatures/"
setwd(Outpath)
library(dummies)
cancerNames <- c("BLCA","BRCA","CESC","ESCA","GBM","HNSC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PAAD","SKCM","PRAD","SARC","STAD","THCA","TGCT","THYM","UCEC")
folder <- "3.6cnv.lesion"
setwd("/extraspace/yye1/analysis/Hypoxia/New/DoubleCheck/3.PSM/")
if (!file.exists(folder)) { dir.create(folder) }
scripts.dir <- paste("/extraspace/yye1/analysis/Hypoxia/New/DoubleCheck/3.PSM/",folder,sep="")
setwd(scripts.dir)
analysis="myclusters" #Oxygen_Content
sum.cnvAll <- data.frame()
####must exist stum  clinical, stratification, cnv files.
for(cancer in cancerNames){
  data <- read.csv(paste("/extraspace/yye1/analysis/Hypoxia/New/ClinicalData/ProcessedClinicalData/",cancer,"_TumorPurity_ClinicalData.csv",sep = ""),header=T)
  rownames(data) <- data[,1]
  data <- data[,-1]
  
  analysis <- "myclusters"
  if(analysis=="myclusters"){
    # convert hypoxic and normoxic to numeric 1,0 to suppress the warning message in lm
    data$myclusters <- ifelse(data$myclusters=="hypoxic",1,0)
    colnames(data)[which(colnames(data)=="myclusters")] <- "Z"
  }
  
  # convert to dummy
  
  dummy.feature <- setdiff(colnames(data),c("Z","age_at_initial_pathologic_diagnosis","Purity"))#,"pathologic_stage"))
  if(length(dummy.feature)>0){
    data.dum <- dummy.data.frame(data, names=dummy.feature)
    dummy.list <- attr(data.dum,"dummies")
    rm.col <- c()
    for (i in 1:length(dummy.list))
    {
      rm.col <- c(rm.col, dummy.list[[i]][length(dummy.list[[i]])])
    }
    data.dum <- data.dum[,-rm.col]
    data.dum$X0 <- rep(1, nrow(data.dum))
    #form <- as.formula("Z~.") # should exclude X0
    exclude.col <- match(c("Z","X0"), colnames(data.dum))
    colnames(data.dum) <- gsub(" ", ".", colnames(data.dum))
    form <- as.formula(paste("Z~",paste(colnames(data.dum)[-exclude.col],collapse="+"),sep=""))
  }else{
    data.dum <- data
    data.dum$X0 <- rep(1, nrow(data.dum))
    #form <- as.formula("Z~.") # should exclude X0
    exclude.col <- match(c("Z","X0"), colnames(data.dum))
    colnames(data.dum) <- gsub(" ", ".", colnames(data.dum))
    form <- as.formula(paste("Z~",paste(colnames(data.dum)[-exclude.col],collapse="+"),sep=""))
  }
  
  # perform calculation
  source("~/code/cal.R")
  library(doMC)
  library(foreach)
  registerDoMC(15)
  # cnv.exp  /extraspace/yye1/share_data/TCGA_CNV_GISTIC/
  cnv <- read.delim(paste("/extraspace/yye1/share_data/TCGA_CNV_GISTIC/GDAC_",cancer,"/all_lesions.conf_99.txt",sep=""),header=T)
  cnv <- cnv[grep("CN values",cnv$Unique.Name),]
  cnv$Unique.Name <- gsub("Peak|- CN values|Peak","",cnv$Unique.Name)
  cnv$Unique.Name <- gsub(" ","",cnv$Unique.Name)
  cnv <- cnv[,c("Unique.Name",grep("TCGA",colnames(cnv),value=T))]
  
  if(length(duplicated(substr(colnames(cnv),1,12))[duplicated(substr(colnames(cnv),1,12))==TRUE]) >0){
    cnv <- cnv[,colnames(cnv)[-which(colnames(cnv) %in% colnames(cnv)[duplicated(substr(colnames(cnv),1,12))])]]
  }
  cnv[,2:ncol(cnv)] <- apply(cnv[,2:ncol(cnv)],2,function(x){2^(1+x)})
  cnv.pri <- data.frame(t(cnv[,2:ncol(cnv)]))
  colnames(cnv.pri) <- as.vector(cnv[,1])
  rownames(cnv.pri) <- substr(rownames(cnv.pri),1,12)
  cnv.pri <- rm.zero.col(cnv.pri)
  
  folder <- paste(cancer,"_",analysis,sep="")
  if (!file.exists(folder)) { dir.create(folder) }
  cnv.result <- weight.test(data.dum, form, cnv.pri, is.continuous=TRUE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=FALSE, cancer, data.type= "cnv", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=FALSE)
  sum.cnv <- summarize.fdr(cnv.pri, cnv.result, print=TRUE)
  summarize.fdr(cnv.pri, cnv.result, print=TRUE, cutoff=0.05)
  write.summary(sum.cnv, cancer, analysis,"cnv")
  write.result(cnv.result, cancer, analysis,"cnv")
  save(cnv.result, file=paste(cancer,"_", analysis,"_result.RData",sep=""))
  if(length(which(cnv.result$fdr < 0.05)) > 0){
    sum.cnv <- data.frame(sum.cnv)
    sum.cnv$class <- rep(cancer,times=nrow(sum.cnv))
    if(nrow(sum.cnvAll) == 0){
      sum.cnvAll <- sum.cnv
    }else{
      sum.cnvAll <- rbind(sum.cnvAll,sum.cnv)
    }
    perm.result <- foreach (seed =1:100, .combine='rbind') %dopar%
    {
      ## cnv, only for KIRC        
      perm.cnv.result <- weight.test(data.dum, form, cnv.pri, is.continuous=TRUE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=FALSE, cancer, "cnv", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=TRUE, seed=seed)
      perm.sum.cnv <- summarize.fdr(cnv.pri, perm.cnv.result)
      
      write(c(seed, perm.sum.cnv$n.sig), file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_count_",seed,".txt", sep=""), ncolumns=7)
      save(seed,perm.cnv.result, file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_",seed,".RData", sep=""))
      
    }
    
    
    cutoff <- 0.05
    seedV <- 1:100
    perm.cal(cancer, analysis, "cnv", cnv.pri, cutoff=cutoff, seedV=seedV)
    
    if(FALSE)
    {
      cnv.ttest <- myttest(data.dum, cnv.pri, cancer,"cnv")
      sum.cnv <- summarize.fdr(cnv.pri, cnv.ttest)
      save(cnv.ttest, file=paste(cancer,"_cnv_ttest.RData", sep=""))
    }
  }
}

write.table(sum.cnvAll,file="cnv.genes.across.cancer.types.txt",quote = F,sep="\t",row.names = F)