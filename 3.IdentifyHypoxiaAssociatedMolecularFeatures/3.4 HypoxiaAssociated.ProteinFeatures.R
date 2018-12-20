#Identify hypoxia associated protein
library(dummies)
folder <- "3.4rppa"
if (!file.exists(folder)) { dir.create(folder) }
scripts.dir <- paste("/extraspace/yye1/analysis/Hypoxia/New/DoubleCheck/3.PSM/",folder,sep="")
setwd(scripts.dir)
analysis="myclusters" #Oxygen_Content
sum.rppaAll <- data.frame()
####must exist stum  clinical, stratification, rppa files.
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
  # rppa.exp #
  rppa <- read.delim(paste("/extraspace/shared_data/TCGA_protein/",cancer,"_protein.re_20160627",sep=""),header=T)
  if(length(duplicated(substr(colnames(rppa),1,12))[duplicated(substr(colnames(rppa),1,12))==TRUE]) >0){
    rppa <- rppa[,colnames(rppa)[-which(colnames(rppa) %in% colnames(rppa)[duplicated(substr(colnames(rppa),1,12))])]]
  }
  rppa.pri <- data.frame(t(rppa[,2:ncol(rppa)]))
  colnames(rppa.pri) <- as.vector(rppa[,1])
  rownames(rppa.pri) <- substr(rownames(rppa.pri),1,12)
  rppa.pri <- rm.zero.col(rppa.pri)
  
  folder <- paste(cancer,"_",analysis,sep="")
  if (!file.exists(folder)) { dir.create(folder) }
  rppa.result <- weight.test(data.dum, form, rppa.pri, is.continuous=TRUE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=FALSE, cancer, data.type= "rppa", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=FALSE)
  sum.rppa <- summarize.fdr(rppa.pri, rppa.result, print=TRUE)
  summarize.fdr(rppa.pri, rppa.result, print=TRUE, cutoff=0.05)
  write.summary(sum.rppa, cancer, analysis,"rppa")
  write.result(rppa.result, cancer, analysis,"rppa")
  save(rppa.result, file=paste(cancer,"_", analysis,"_result.RData",sep=""))
  if(length(which(rppa.result$fdr < 0.05)) > 0){
    sum.rppa <- data.frame(sum.rppa)
    sum.rppa$class <- rep(cancer,times=nrow(sum.rppa))
    if(nrow(sum.rppaAll) == 0){
      sum.rppaAll <- sum.rppa
    }else{
      sum.rppaAll <- rbind(sum.rppaAll,sum.rppa)
    }
    perm.result <- foreach (seed =1:100, .combine='rbind') %dopar%
    {
      ## rppa, only for KIRC        
      perm.rppa.result <- weight.test(data.dum, form, rppa.pri, is.continuous=TRUE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=FALSE, cancer, "rppa", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=TRUE, seed=seed)
      perm.sum.rppa <- summarize.fdr(rppa.pri, perm.rppa.result)
      
      write(c(seed, perm.sum.rppa$n.sig), file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_count_",seed,".txt", sep=""), ncolumns=7)
      save(seed,perm.rppa.result, file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_",seed,".RData", sep=""))
      
    }
    
    
    cutoff <- 0.05
    seedV <- 1:100
    perm.cal(cancer, analysis, "rppa", rppa.pri, cutoff=cutoff, seedV=seedV)
    
    if(FALSE)
    {
      rppa.ttest <- myttest(data.dum, rppa.pri, cancer,"rppa")
      sum.rppa <- summarize.fdr(rppa.pri, rppa.ttest)
      save(rppa.ttest, file=paste(cancer,"_rppa_ttest.RData", sep=""))
    }
  }
}


write.table(sum.rppaAll,file="rppa.genes.across.cancer.typesAll.txt",quote = F,sep="\t",row.names = F)


