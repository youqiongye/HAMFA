###Identify hypoxia associated miRNA
Outpath <- "/extraspace/yye1/analysis/Hypoxia/New/DoubleCheck/HypoxiaAssociatedFeatures/"
setwd(Outpath)

library(dummies)
folder <- "3.5miRNAseq"
if (!file.exists(folder)) { dir.create(folder) }
scripts.dir <- paste("/extraspace/yye1/analysis/Hypoxia/New/DoubleCheck/3.PSM/",folder,sep="")
setwd(scripts.dir)
analysis="myclusters" #Oxygen_Content
sum.miRNAseqAll <- data.frame()
####must exist stum  clinical, stratification, miRNAseq files.
cancerNames <- c("BLCA","BRCA","CESC","ESCA","GBM","HNSC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PAAD","SKCM","PRAD","SARC","STAD","THCA","TGCT","THYM","UCEC")

#stratification.files.Abs <- stratification.files.Abs[which(stratification.files.Abs != "OV" & stratification.files.Abs !="UCEC")]
for(cancer in cancerNames){
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
  registerDoMC(15)
  # miRNAseq.exp #
  miRNAseq <- read.delim(paste("/extraspace/TCGA/TCGA_exp_DataPortal/miRNA_exp/",cancer,"_miRNA_each_exp_20160513",sep=""),header=T)
  if(ncol(miRNAseq) >= 50){
  if(cancer=="SKCM"){
    selectCols <- c("gene",colnames( miRNAseq )[as.numeric(substr(colnames( miRNAseq ),14,15)) ==6])
    selectCols <- selectCols[!is.na(selectCols)]
    miRNAseq <-  miRNAseq[,selectCols]
  }else{
    selectCols <- c("gene",colnames( miRNAseq )[as.numeric(substr(colnames( miRNAseq ),14,15)) ==1])
    selectCols <- selectCols[!is.na(selectCols)]
    miRNAseq <-  miRNAseq[,selectCols]
  }
  
    if(length(duplicated(substr(colnames(miRNAseq),1,12))[duplicated(substr(colnames(miRNAseq),1,12))==TRUE]) >0){
      miRNAseq <- miRNAseq[,colnames(miRNAseq)[-which(colnames(miRNAseq) %in% colnames(miRNAseq)[duplicated(substr(colnames(miRNAseq),1,12))])]]
    }
    miRNAseq[,2:ncol(miRNAseq)] <- miRNAseq[,2:ncol(miRNAseq)]
    miRNAseq.pri <- data.frame(t(miRNAseq[,2:ncol(miRNAseq)]))
    miRNAseq.pri <- apply(miRNAseq.pri,2,function(x)log2(x+1))
    colnames(miRNAseq.pri) <- as.vector(miRNAseq[,1])
    rownames(miRNAseq.pri) <- substr(rownames(miRNAseq.pri),1,12)
    miRNAseq.pri <- rm.zero.col(miRNAseq.pri)
    
    folder <- paste(cancer,"_",analysis,sep="")
    if (!file.exists(folder)) { dir.create(folder) }
    miRNAseq.result <- weight.test(data.dum, form, miRNAseq.pri, is.continuous=TRUE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=FALSE, cancer, data.type= "miRNAseq", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=FALSE)
    sum.miRNAseq <- summarize.fdr(miRNAseq.pri, miRNAseq.result, print=TRUE)
    summarize.fdr(miRNAseq.pri, miRNAseq.result, print=TRUE, cutoff=0.05)
    write.summary(sum.miRNAseq, cancer, analysis,"miRNAseq")
    write.result(miRNAseq.result, cancer, analysis,"miRNAseq")
    save(miRNAseq.result, file=paste(cancer,"_", analysis,"_result.RData",sep=""))
    if(length(which(miRNAseq.result$fdr < 0.05)) > 0){
      sum.miRNAseq <- data.frame(sum.miRNAseq)
      sum.miRNAseq$class <- rep(cancer,times=nrow(sum.miRNAseq))
      if(nrow(sum.miRNAseqAll) == 0){
        sum.miRNAseqAll <- sum.miRNAseq
      }else{
        sum.miRNAseqAll <- rbind(sum.miRNAseqAll,sum.miRNAseq)
      }
      perm.result <- foreach (seed =1:100, .combine='rbind') %dopar%
      {
        ## miRNAseq, only for KIRC        
        perm.miRNAseq.result <- weight.test(data.dum, form, miRNAseq.pri, is.continuous=TRUE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=FALSE, cancer, "miRNAseq", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=TRUE, seed=seed)
        perm.sum.miRNAseq <- summarize.fdr(miRNAseq.pri, perm.miRNAseq.result)
        
        write(c(seed, perm.sum.miRNAseq$n.sig), file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_count_",seed,".txt", sep=""), ncolumns=7)
        save(seed,perm.miRNAseq.result, file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_",seed,".RData", sep=""))
        
      }
      
      
      cutoff <- 0.05
      seedV <- 1:100
      perm.cal(cancer, analysis, "miRNAseq", miRNAseq.pri, cutoff=cutoff, seedV=seedV)
      
      if(FALSE)
      {
        miRNAseq.ttest <- myttest(data.dum, miRNAseq.pri, cancer,"miRNAseq")
        sum.miRNAseq <- summarize.fdr(miRNAseq.pri, miRNAseq.ttest)
        save(miRNAseq.ttest, file=paste(cancer,"_miRNAseq_ttest.RData", sep=""))
      }
    }
  }
  
}



write.table(sum.miRNAseqAll,file="miRNAseq.genes.across.cancer.typesAll.txt",quote = F,sep="\t",row.names = F)



