###Identify hypoxia associated mutation
Outpath <- "/extraspace/yye1/analysis/Hypoxia/New/DoubleCheck/HypoxiaAssociatedFeatures/"
setwd(Outpath)
library(dummies)
folder <- "3.3Mutation"
if (!file.exists(folder)) { dir.create(folder) }
scripts.dir <- paste("/extraspace/yye1/analysis/Hypoxia/New/DoubleCheck/3.PSM/",folder,sep="")
setwd(scripts.dir)
analysis="myclusters" #Oxygen_Content
sum.mutAll <- data.frame()
cancerNames <- c("BLCA","BRCA","CESC","ESCA","GBM","HNSC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PAAD","SKCM","PRAD","SARC","STAD","THCA","TGCT","THYM","UCEC")

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
  # mutation
  mutation <- read.delim(paste("/extraspace/TCGA/Mutation/mutation_",cancer,".txt",sep=""),header=T)
  mutation.pri <- data.frame(t(mutation[,2:ncol(mutation)]))
  colnames(mutation.pri) <- as.vector(mutation[,1])
  # only focuse on the highly mutated genes
  # if the top 100 mutated gene mutation frequency less than 0.05, keep top 100 mutated genes, or keep all genes which matutation frequency large than 0.05
  top <- 100
  mut.cutoff <- 0.05
  mut.rate <- colSums(mutation.pri)/nrow(mutation.pri)
  mut.rate.sort <- sort(mut.rate, decreasing=T, index.return=T)
  if(length(which(mut.rate > mut.cutoff)) > 100){
    keep.index <- mut.rate.sort$ix[match( names(which(mut.rate.sort$x > mut.cutoff)),names(mut.rate.sort$x))]
  }else{
    top.index <-mut.rate.sort$ix[1:top]
    top.cutoff <- mut.rate.sort$x[top] # keep ties
    keep.index <- which(mut.rate>= top.cutoff)
  }
  #folder <- "sample_list"
  #if (!file.exists(folder)) { dir.create(folder) }
  # For LUSC: gender, ATT/IPW could not get balanced weight. OVERLAP and MW get balanced results. Double check!
  folder <- paste(cancer,"_",analysis,sep="")
  if (!file.exists(folder)) { dir.create(folder) }
  mut.result <- weight.test(data.dum, form, mutation.pri[,keep.index], is.continuous=FALSE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=FALSE, cancer, "mut", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=FALSE)
  sum.mut <- summarize.fdr(mutation.pri[,keep.index], mut.result, print=TRUE)
  summarize.fdr(mutation.pri[,keep.index], mut.result, print=TRUE, cutoff=0.05)
  write.summary(sum.mut, cancer, analysis,"mut")
  write.result(mut.result, cancer, analysis,"mut")
  if(length(which(mut.result$fdr < 0.05)) > 0){
    sum.mut <- data.frame(sum.mut)
    sum.mut$class <- rep(cancer,times=nrow(sum.mut))
    if(nrow(sum.mutAll) == 0){
      sum.mutAll <- sum.mut
    }else{
      sum.mutAll <- rbind(sum.mutAll,sum.mut)
    }
    perm.result <- foreach (seed =1:100, .combine='rbind') %dopar%
    {      
      perm.mut.result <- weight.test(data.dum, form, mutation.pri[,keep.index], is.continuous=FALSE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=TRUE, cancer, "mut", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=TRUE, seed=seed)
      perm.sum.mut <- summarize.fdr(mutation.pri[,keep.index], perm.mut.result)
      
      write(c(seed, perm.sum.mut$n.sig), file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_count_",seed,".txt", sep=""), ncolumns=7)
      save(seed,perm.mut.result, file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_",seed,".RData", sep=""))
      
    }
    
    
    cutoff <- 0.05
    seedV <- 1:100
    perm.cal(cancer, analysis, "mut", mutation.pri[,keep.index], cutoff=cutoff, seedV=seedV)
    
    if(FALSE)
    {
      mut.ttest <- myttest(data.dum, mutation.pri, cancer,"mut")
      sum.mut <- summarize.fdr(mutation.pri, mut.ttest)
      save(mRNAseq.ttest, file=paste(cancer,"_mut_ttest.RData", sep=""))
    }
  }
}


write.table(sum.mutAll,file="Mutation.genes.across.cancer.types.txt",quote = F,sep="\t",row.names = F)





