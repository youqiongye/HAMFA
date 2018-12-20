#Identify hypoxia associated mRNA
library(dummies)
Outpath <- "/extraspace/yye1/analysis/Hypoxia/New/DoubleCheck/HypoxiaAssociatedFeatures/"
setwd(Outpath)
folder <- "3.1mRNA"
if (!file.exists(folder)) { dir.create(folder) }
scripts.dir <- paste(Outpath,folder,sep="")
setwd(scripts.dir)
analysis="myclusters" #Oxygen_Content
sum.mRNAAll <- data.frame()
####must exist stum  clinical, stratification, mRNA files.
cancerNames <- c("BLCA","BRCA","CESC","ESCA","GBM","HNSC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PAAD","SKCM","PRAD","SARC","STAD","THCA","TGCT","THYM","UCEC")
for(cancer in cancerNames){
  data <- read.csv(paste("/extraspace/yye1/analysis/Hypoxia/New/DoubleCheck/2.ClinicalData/ProcessedClinicalData/",cancer,"_TumorPurity_ClinicalData.csv",sep = ""),header=T)
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
  registerDoMC(21)
  # mRNA.exp
  mRNAseq <- read.delim(paste("/extraspace/yye1/share_data/TCGA_mRNAlog2/",cancer,"_mRNAlog2.tab",sep=""),header=T)
  mRNAseq <-  mRNAseq[,c("gene",colnames( mRNAseq )[as.numeric(substr(colnames( mRNAseq ),14,15)) %in% c(1,6)])]
  
  mRNAseq <- mRNAseq[,colnames(mRNAseq)[!duplicated(substr(colnames(mRNAseq),1,12))]]
  mRNAseq.pri <- mRNAseq[,2:ncol(mRNAseq)]
  colnames(mRNAseq.pri) <- substr(colnames(mRNAseq.pri),1,12)
  rownames(mRNAseq.pri) <- mRNAseq$gene
  mRNAseq.pri <- t(mRNAseq.pri)
  mRNAseq.pri <- rm.zero.col(mRNAseq.pri)
  folder <- paste(cancer,"_",analysis,sep="")
  if (!file.exists(folder)) { dir.create(folder) }
  
  mRNAseq.result <- weight.test(data.dum, form, mRNAseq.pri, is.continuous=TRUE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=FALSE, cancer, data.type= "mRNAseq", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=FALSE)
  sum.mRNA <- summarize.fdr(mRNAseq.pri, mRNAseq.result, print=TRUE)
  summarize.fdr(mRNAseq.pri, mRNAseq.result, print=TRUE, cutoff=0.05)
  write.summary(sum.mRNA, cancer, analysis,"mRNA")
  write.result(mRNAseq.result, cancer, analysis,"mRNA")
  save(mRNAseq.result, file=paste(cancer,"_", analysis,"_result.RData",sep=""))
  if(length(which(mRNAseq.result$fdr < 0.05)) > 0){
    sum.mRNA <- data.frame(sum.mRNA)
    sum.mRNA$class <- rep(cancer,times=nrow(sum.mRNA))
    if(nrow(sum.mRNAAll) == 0){
      sum.mRNAAll <- sum.mRNA
    }else{
      sum.mRNAAll <- rbind(sum.mRNAAll,sum.mRNA)
    }
    perm.result <- foreach (seed =1:100, .combine='rbind') %dopar%
    {
      ## mRNAseq, only for KIRC        
      perm.mRNAseq.result <- weight.test(data.dum, form, mRNAseq.pri, is.continuous=TRUE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=FALSE, cancer, "mRNA", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=TRUE, seed=seed)
      perm.sum.mRNA <- summarize.fdr(mRNAseq.pri, perm.mRNAseq.result)
      
      write(c(seed, perm.sum.mRNA$n.sig), file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_count_",seed,".txt", sep=""), ncolumns=7)
      save(seed,perm.mRNAseq.result, file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_",seed,".RData", sep=""))
      
    }
    
    
    cutoff <- 0.05
    seedV <- 1:100
    perm.cal(cancer, analysis, "mRNAseq", mRNAseq.pri, cutoff=cutoff, seedV=seedV)
    
    if(FALSE)
    {
      mRNA.ttest <- myttest(data.dum, mRNAseq.pri, cancer,"mRNA")
      sum.mRNA <- summarize.fdr(mRNAseq.pri, mRNA.ttest)
      save(mRNAseq.ttest, file=paste(cancer,"_mRNA_ttest.RData", sep=""))
    }
  }
}

write.table(sum.mRNAAll,file="mRNAseqlog2.genes.across.cancer.typesAll.txt",quote = F,sep="\t",row.names = F)

