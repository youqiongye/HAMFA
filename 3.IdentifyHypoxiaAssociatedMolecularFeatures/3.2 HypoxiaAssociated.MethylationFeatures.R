##Identify hypoxia associated DNA methylation
Outpath <- "/extraspace/yye1/analysis/Hypoxia/New/DoubleCheck/HypoxiaAssociatedFeatures/"
setwd(Outpath)
library(dummies)
folder <- "3.2methy"
if (!file.exists(folder)) { dir.create(folder) }
scripts.dir <- paste("/extraspace/yye1/analysis/Hypoxia/New/DoubleCheck/3.PSM/",folder,sep="")
setwd(scripts.dir)
analysis="myclusters" #Oxygen_Content
sum.methyAll <- data.frame()
####must exist stum  clinical, stratification, methy files.
cancerNames <- c("BLCA","BRCA","CESC","ESCA","GBM","HNSC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PAAD","SKCM","PRAD","SARC","STAD","THCA","TGCT","THYM","UCEC")
for(cancer in cancerNames){
  methy <- readRDS(paste("/extraspace/yye1/share_data/TCGA_methy_mostNegatively/",cancer,"_methy_probe_mostNegatively_to_GeneExp.rds",sep=""))
  if(ncol(methy) > 50){
    data <- read.csv(paste("/extraspace/yye1/analysis/Hypoxia/New/ClinicalData/ProcessedClinicalData/",cancer,"_TumorPurity_ClinicalData.csv",sep = ""),header=T)
    
    data <- unique(data)
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
    registerDoMC(20)
    # methy.exp
    #  methy <- readRDS(paste("/extraspace/yye1/share_data/TCGA_methy_mostNegatively/",cancer,"_methy_probe_mostNegatively_to_GeneExp.rds",sep=""))
    methy <- methy[,c("gene",colnames( methy )[as.numeric(substr(colnames( methy ),14,15)) %in% c(1,6)])]
    methy <- methy[,colnames(methy)[!duplicated(substr(colnames(methy),1,12))]]
    methy.pri <- data.frame(t(methy[,2:ncol(methy)]))
    colnames(methy.pri) <- as.vector(methy[,1])
    rownames(methy.pri) <- substr(rownames(methy.pri),1,12)
    methy.pri <- rm.zero.col(methy.pri)
    library(randomForest)
    methy.pri <- na.roughfix(methy.pri)
    
    folder <- paste(cancer,"_",analysis,sep="")
    if (!file.exists(folder)) { dir.create(folder) }
    methy.result <- weight.test(data.dum, form, methy.pri, is.continuous=TRUE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=FALSE, cancer, data.type= "methy", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=FALSE)
    sum.methy <- summarize.fdr(methy.pri, methy.result, print=TRUE)
    summarize.fdr(methy.pri, methy.result, print=TRUE, cutoff=0.05)
    write.summary(sum.methy, cancer, analysis,"methy")
    write.result(methy.result, cancer, analysis,"methy")
    save(methy.result, file=paste(cancer,"_", analysis,"_result.RData",sep=""))
    if(length(which(methy.result$fdr < 0.05)) > 0){
      sum.methy <- data.frame(sum.methy)
      sum.methy$class <- rep(cancer,times=nrow(sum.methy))
      if(nrow(sum.methyAll) == 0){
        sum.methyAll <- sum.methy
      }else{
        sum.methyAll <- rbind(sum.methyAll,sum.methy)
      }
      perm.result <- foreach (seed =1:100, .combine='rbind') %dopar%
      {
        ## methy, only for KIRC        
        perm.methy.result <- weight.test(data.dum, form, methy.pri, is.continuous=TRUE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=FALSE, cancer, "methy", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=TRUE, seed=seed)
        perm.sum.methy <- summarize.fdr(methy.pri, perm.methy.result)
        
        write(c(seed, perm.sum.methy$n.sig), file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_count_",seed,".txt", sep=""), ncolumns=7)
        save(seed,perm.methy.result, file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_",seed,".RData", sep=""))
        
      }
      
      
      cutoff <- 0.05
      seedV <- 1:100
      perm.cal(cancer, analysis, "methy", methy.pri, cutoff=cutoff, seedV=seedV)
      
      if(FALSE)
      {
        methy.ttest <- myttest(data.dum, methy.pri, cancer,"methy")
        sum.methy <- summarize.fdr(methy.pri, methy.ttest)
        save(methy.ttest, file=paste(cancer,"_methy_ttest.RData", sep=""))
      }
    }
  }
  
}


write.table(sum.methyAll,file="methy.genes.across.cancer.typesAll.txt",quote = F,sep="\t",row.names = F)




