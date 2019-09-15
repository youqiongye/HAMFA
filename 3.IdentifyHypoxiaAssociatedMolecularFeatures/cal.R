# weighted test agains continuous response variables, such as gene expression


# drawback: the ylabel of the inversed plot is negative, can be edited in AI
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
hist.mirror2 <- function(data.plot, label)
    {
        myplot <- ggplot(data.plot, aes(ps, fill=status))+geom_histogram(data=subset(data.plot, status=="Yes"), aes(ps, y=..count..))+
            geom_histogram(data=subset(data.plot,status=="No"), aes(ps, y= -..count..))+
                scale_fill_hue(label)+
                    labs(x="Propensity Score")
        print(myplot)
    }


weight.test <- function(data,form, molecular.pri, is.continuous=TRUE, weight="MW", mirror.plot=FALSE,cancer, data.type, outdir=".", perm=FALSE, seed=seed)
    {  #molecular.pri <- mutation.pri[,keep.index]
        rownames(molecular.pri) <- gsub("\\.","\\-",rownames(molecular.pri))    
        common <- intersect(rownames(data), rownames(molecular.pri))
        print(paste("Number of samples:", length(common)))
        
        molecular.common <- molecular.pri[match(common, rownames(molecular.pri)),]
        clinical.common <- data[match(common, rownames(data)),]
        print(summary(factor(clinical.common$Z)))
        # write sample info to file
        write.table(t(rbind(rownames(clinical.common),ifelse(clinical.common$Z==1, "hypoxic","normoxic"))), file=paste(cancer,"_",data.type,"_samples.txt", sep=""),quote = F,row.names = F,sep="\t")
        
        print("----------------")
        if("age_at_initial_pathologic_diagnosis" %in% colnames(clinical.common)){
          age.cutoff <- 65
          print(paste("age <", age.cutoff,":", length(which(clinical.common$age_at_initial_pathologic_diagnosis< age.cutoff))))
          print(paste("age >=", age.cutoff,":", length(which(clinical.common$age_at_initial_pathologic_diagnosis>= age.cutoff))))
          print("----------------")
          
        }
  

        if(perm)
            {
                n <- nrow(clinical.common)
                set.seed(seed)
                perm <- sample(1:n,n)
                clinical.common$Z <- clinical.common$Z[perm]
            }
 
        print(paste("Weighting scheme:", weight))
        source("GIPW_function_omega.R")
        ans <- GIPW.std.omega(dat=clinical.common, form.ps=form, weight=weight,trt.est=F)
        
        # draw mirror plot

        if(mirror.plot)
            {
                ps <- ans$ps
                status <- ifelse(clinical.common$Z==1, "Yes","No")
                data.plot <- data.frame(ps, status )
                pdf(paste(outdir,"/",cancer,"_", data.type, ifelse(perm,paste("_perm_",seed, sep=""), ""),"_mirror_raw.pdf",sep=""))
                if(analysis=="omt")
                    {
                        label="OMT"
                    }
                if(analysis=="gender")
                    {
                        label="Female"
                    }
                if(analysis=="race")
                    {
                        label="Non-White"
                }
                if(analysis=="myclusters"){
                  label = "Hypoxic"
                }
                print(hist.mirror2(data.plot, label))
                dev.off()
            }
        wt <- ans$W
        source("check_balance.R")
        index.tr <- which(clinical.common$Z==1)
        index.ctr <- which(clinical.common$Z==0)

        cutoff=0.1 # 10%
        for (i in 2:(ncol(clinical.common)-1))
            {
                 print(paste("check", colnames(clinical.common)[i]))
                if(nrow(clinical.common[which(clinical.common[,i]==1),]) >=3){
                  #print(summary(clinical.common[,i]))
                  std.diff <- check.balance(index.tr, index.ctr, clinical.common[,i], wt, printout=TRUE)
                  if(std.diff>cutoff)
                  {
                    print(paste("Fail: standardized difference= ", std.diff))
                    print("=======================================")
                    #  stop()
                  }else{
                    print("Pass!")
                  }
                }else{
                  print("too small sample size")
                }
                
               
            }

        
        molecular.pvalues <- c()
        molecular.coefs <- c()
        molecular.0 <- c()
        molecular.1 <- c()
        molecular.0.w <- c()
        molecular.1.w <- c()
        for (i in 1:ncol(molecular.common))
            {
                tmp <- sapply(split(as.numeric(molecular.common[,i]), clinical.common$Z), mean,na.rm=T)
                molecular.0 <- c(molecular.0, tmp[1])
                molecular.1 <- c(molecular.1, tmp[2])
                # get weighted mean
                molecular.0.w <- c(molecular.0.w, sum(as.numeric(molecular.common[index.ctr,i])*wt[index.ctr],na.rm=T)/sum(wt[index.ctr],na.rm=T))
                molecular.1.w <- c(molecular.1.w, sum(as.numeric(molecular.common[index.tr,i])*wt[index.tr],na.rm=T)/sum(wt[index.tr],na.rm=T))
                #print(i)
                if (i%%1000==0)
                    {
                        print(i)
                    }
                
                if(is.continuous)
                    {
                        #pvalue <- try(summary(lm(clinical.common$Z~molecular.common[,i], weights=wt))$coef[2,4])
                        pvalue <- try(summary(lm(molecular.common[,i]~clinical.common$Z, weights=wt))$coef[2,4])
                        coef <- try(summary(lm(molecular.common[,i]~clinical.common$Z, weights=wt))$coef[2,1])
                        
                        if (class(pvalue)=="try-error")
                            {
                                pvalue <- NA
                                coef <- NA
                            }
                    }else{
                        
                        #pvalue <- summary(glm(clinical.common$Z~molecular.common[,i], family=binomial, weights=wt))$coef[2,4]
                        pvalue <- summary(glm(as.numeric(molecular.common[,i])~clinical.common$Z, family=binomial, weights=wt),na.rm=T)$coef[2,4]
                        coef <- summary(glm(as.numeric(molecular.common[,i])~clinical.common$Z, family=binomial, weights=wt),na.rm=T)$coef[2,1]
                        if (class(pvalue)=="try-error")
                            {
                                pvalue <- NA
                                coef <- NA
                            }
                    }
                molecular.pvalues <- c(molecular.pvalues, pvalue)
                molecular.coefs <- c(molecular.coefs, coef)
            }
        molecular.fdr <- p.adjust(molecular.pvalues,"fdr")
       
         return(list(feature=colnames(molecular.common),coef=molecular.coefs, pvalue=molecular.pvalues,fdr=molecular.fdr , mean.0= molecular.0, mean.1= molecular.1, mean.0.w=molecular.0.w, mean.1.w=molecular.1.w))
    }

summarize.fdr <- function(molecular.data, molecular.result, cutoff=0.05, print=FALSE )
    {
        pvalue <- molecular.result$pvalue
        fdr <- molecular.result$fdr
        coef <- molecular.result$coef
        mean.0 <- molecular.result$mean.0
        mean.1 <- molecular.result$mean.1
        mean.0.w <- molecular.result$mean.0.w
        mean.1.w <- molecular.result$mean.1.w
        signif.index <- which(fdr<cutoff)
        print(paste("Features with FDR <", cutoff, "=", length(signif.index)))
        if(print)
            {
                print(cbind(colnames(molecular.data)[signif.index],signif(pvalue[signif.index],3), signif(fdr[signif.index],3), signif(coef[signif.index],3), signif(mean.0[signif.index],3), signif(mean.1[signif.index],3)))
            }
        return(list(feature.sig=colnames(molecular.data)[signif.index], pvalue.sig=pvalue[signif.index], fdr.sig= fdr[signif.index], n.sig=length(signif.index), coef.sig=coef[signif.index], mean0.sig=mean.0[signif.index], mean1.sig=mean.1[signif.index], mean0.sig.w=mean.0.w[signif.index], mean1.sig.w=mean.1.w[signif.index]))#, n.sig=length(which(molecular.result$pvalue<cutoff))) )
        
    }


rm.zero.col <- function(molecular.data)
    {       
        ## remove non-coding genes
        genes.raw <- colnames(molecular.data)
        genes <- unlist(lapply(genes.raw, function(x) unlist(strsplit(x,"[|]"))[1]))
        ids <- as.numeric(unlist(lapply(genes.raw, function(x) unlist(strsplit(x,"[|]"))[2])))
        if(FALSE) # do not remove any gene
            {
        rm1 <- grep("?", genes, fixed=TRUE) # question mark
        #rm2 <- grep("^MIR", genes) # MIR
        #rm3 <- grep("^SNOR", genes)
        #rm <- c(rm1, rm2, rm3)
        molecular.data <- molecular.data[,-rm1]
        print(paste(length(rm1),"genes were removed."))
    }
        
        
        # remove non-variable genes, i.e., sd=0
        zero.col <- which(apply(molecular.data, 2, sd)==0)
        if (length(zero.col)>0)
            {
                molecular.data <- molecular.data[,-zero.col]
                print(paste(length(zero.col),"columns were removed because of no variation."))
            }
        
        return(molecular.data)
    }

write.summary <- function(summary, cancer, analysis, type)
    {
        file <- paste(cancer,"_",analysis,"_summary_",type,".txt", sep="")
        data <- data.frame(summary$feature.sig, summary$fdr.sig, summary$coef.sig, summary$mean0.sig, summary$mean1.sig, summary$mean0.sig.w, summary$mean1.sig.w)
        if (analysis=="gender")
            {
                colnames(data) <- c("feature","fdr","coef","mean_MALE", "mean_FEMALE", "mean_MALE_weighted", "mean_FEMALE_weighted")
            }
        if(analysis=="race")
            {
                colnames(data) <- c("feature","fdr","coef","mean_nonWHITE", "mean_WHITE","mean_nonWHITE_weighted", "mean_WHITE_weighted")
            }
        if (analysis=="omt")
            {
                colnames(data) <- c("feature","fdr","coef","mean_nonOMT", "mean_OMT","mean_nonOMT_weighted", "mean_OMT_weighted")
            }
        write.table(data,file=file, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
    }

write.result <- function(result, cancer, analysis, type, fdr.cutoff=0.05)
    {
        file <- paste(cancer,"_",analysis,"_result_",type,".txt", sep="")
        data <- data.frame(result$feature, result$pvalue, result$fdr, result$coef, result$mean.0, result$mean.1, result$mean.0.w, result$mean.1.w)
        if (analysis=="gender")
            {
                colnames(data) <- c("feature","pvalue","fdr","coef","mean_MALE", "mean_FEMALE", "mean_MALE_weighted", "mean_FEMALE_weighted")
            }
        if(analysis=="race")
            {
                colnames(data) <- c("feature","fdr","coef","mean_nonWHITE", "mean_WHITE","mean_nonWHITE_weighted", "mean_WHITE_weighted")
            }
        if (analysis=="omt")
            {
                colnames(data) <- c("feature","fdr","coef","mean_nonOMT", "mean_OMT","mean_nonOMT_weighted", "mean_OMT_weighted")
            }
        write.table(data,file=file, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)      
    }


plot.perm <- function(perm, n, cancer, analysis, type, cutoff)
    {        
        p=length(which(perm>= n))/length(perm)
        print(paste("Permutation P-value =", p))
        print(paste(type,": n.sig =", n))
        print(paste("Median perm n.sig =", median(perm)))
        file <- paste(cancer,"_",analysis,"_perm_",type,"_",cutoff,".pdf", sep="")
        pdf(file, width=3, height=3)
        par(mgp=c(2,1,0))
        if(type=="mut") {main="Mutation"}
        if(type=="cnv") {main="SCNA"}
        if(type=="methy") {main="Methy"}
        if(type=="mRNAseq") {main="mRNA"}
        if(type=="miRNA") {main="miRNA"}
        if(type=="rppa") {main="Protein"}
        if (n > max(perm))
            {
                myhist <- hist(perm, xlim=c(0, max(max(perm),n)), main="", xlab="# Genes" )
                
            }else{
                myhist <- hist(perm, main="", xlab="# Genes")
            }
        text(n, max(myhist$counts),paste("p-value =", p), pos=ifelse(n> 0.5* max(perm),2,4), col=ifelse(p<=0.05, "red","blue"), cex=1.1, xpd=T)
        abline(v=n, col=ifelse(p<=0.05, "red","blue"), lty=2, lwd=2)
        dev.off()
    }

perm.cal <- function(cancer, analysis, type, molecular.pri, cutoff=0.05, seedV=1:100)
    {
        print(paste("Calculate permutation for", cancer,":",type,":"))
      #  load(file=paste(cancer,"_", analysis,"_result.RData",sep=""))
        result <- get(paste(type,".result", sep=""))
        summary <- summarize.fdr(molecular.pri, result, cutoff=cutoff )
        n <- summary$n.sig
        perm <- c()
        print("For permutation:")
        for (seed in seedV)
            {
                ##print(seed)
                perm.result <- perm.summary <- c()
                 if(type=="mut")
                     {
                         load(paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_",seed,".RData", sep=""))
                     }else if(type=="pre")
                          {
                              load(paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_pre_",seed,".RData", sep=""))
                          }else{
                 
                              load(paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_",seed,".RData", sep=""))
                          }
                perm.result <- get(paste("perm.",type,".result", sep=""))
                perm.summary <- summarize.fdr(molecular.pri, perm.result, cutoff=cutoff )
                perm <- c(perm, perm.summary$n.sig)
                do.call(rm, list(paste("perm.",type,".result", sep="")))
                               
            }
        ##print (perm)
        ##print (n)
        plot.perm(perm, n, cancer, analysis, type, cutoff)
    }

# only for mRNAseq, only write the protein coding genes, adjust fdr accordingly
keep.coding <- function(mRNAseq.result,cancer)
    {
        genes.raw <- mRNAseq.result$feature
        genes <- unlist(lapply(genes.raw, function(x) unlist(strsplit(x,"[|]"))[1]))
        ids <- as.numeric(unlist(lapply(genes.raw, function(x) unlist(strsplit(x,"[|]"))[2])))
        gene.data <- read.table("/rsrch1/bcb/yyuan3/extraspace/yyuan3/Lab_share/protein_coding_gene_annotation.tsv",header=TRUE, sep="\t", quote="")
        name.id <- paste(gene.data$Gene, gene.data$ID, sep="|")
        keep.index <- match(gene.data$ID, ids)
        print(paste("After remove non-coding, total genes:", length(keep.index)))
        out <- data.frame(mRNAseq.result)[keep.index,]
        #out$feature <- as.vector(out$feature)
        out$fdr <- p.adjust(out$pvalue, "fdr")
        colnames(out) <- c("feature","coef","pvalue","fdr","mean_MALE", "mean_FEMALE", "mean_MALE_weighted", "mean_FEMALE_weighted")
        write.table(out, file=paste(cancer,"_gender_result_mRNAseq_coding.txt", sep=""), sep="\t", quote=FALSE,row.names=FALSE)
        
        
    }

## t-test
myttest <- function(data, molecular.pri, cancer, data.type, perm=FALSE, seed=seed)
    {
        
      common <- intersect(rownames(data), rownames(molecular.pri))
        print(paste("Number of samples:", length(common)))
        molecular.common <- molecular.pri[match(common, rownames(molecular.pri)),]
        clinical.common <- data[match(common, rownames(data)),]
        print(summary(factor(clinical.common$Z)))
        print("----------------")

        if(perm)
            {
                n <- nrow(clinical.common)
                set.seed(seed)
                perm <- sample(1:n,n)
                clinical.common$Z <- clinical.common$Z[perm]
            }
        
        molecular.pvalues <- c()
        molecular.coefs <- c()
        molecular.0 <- c()
        molecular.1 <- c()
         for (i in 1:ncol(molecular.common))
             {
                 tmp <- by(molecular.common[,i], clinical.common$Z, mean)
                 molecular.0 <- c(molecular.0, tmp[1])
                 molecular.1 <- c(molecular.1, tmp[2])
                 if (i%%1000==0)
                    {
                        print(i)
                    }
                 result <- try(t.test(molecular.common[,i]~clinical.common$Z))
                 if (class(result)=="try-error")
                     {
                         pvalue <- coef <- NA
                     }else{
                         
                         pvalue <- result$p.value
                         coef <- result$statistic
                     }
                 molecular.pvalues <- c(molecular.pvalues, pvalue)
                 molecular.coefs <- c(molecular.coefs, coef)
             }
        molecular.fdr <- p.adjust(molecular.pvalues,"fdr")
        return(list(feature=colnames(molecular.common),coef=molecular.coefs, pvalue=molecular.pvalues,fdr=molecular.fdr , mean.0= molecular.0, mean.1= molecular.1))
        
    }

format.barcode <- function(samples)
    {
            
        samples.new <- unlist(lapply(samples, function(x) paste(unlist(strsplit(x, "-"))[1:3], collapse="-")))
        types <- unlist(lapply(samples, function(x) paste(unlist(strsplit(x, "-"))[4], collapse="-")))
        print(summary(as.factor(types)))
        types[substr(types, 1,2)=="01" | substr(types, 1,2)=="03"] <- "primary"
        types[substr(types, 1,2)=="02"] <- "recurrent"
        types[substr(types, 1,2)=="06"] <- "metastatic"
        types[substr(types, 1,2)=="11"| substr(types, 1,2)=="10"] <- "normal"
        return(list(barcodes=samples.new, types=types))
        
    }

format.mature <- function(mature.raw)
    {
        samples <- format.barcode(rownames(mature.raw))
        ## only primary cancer will be used for downstream analysis
        ## KIRC has two 05: additinoal- new primary
        mature.pri <- mature.raw[samples$types=="primary",,drop=FALSE]
        mature.rec <- mature.raw[samples$types=="recurrent",,drop=FALSE]
        mature.met <- mature.raw[samples$types=="metastatic",,drop=FALSE]
        mature.normal <- mature.raw[samples$types=="normal",,drop=FALSE]
        rownames(mature.pri) <- samples$barcodes[samples$type=="primary"]
        rownames(mature.rec) <- samples$barcodes[samples$type=="recurrent"]
        rownames(mature.met) <- samples$barcodes[samples$type=="metastatic"]
        rownames(mature.normal) <- samples$barcodes[samples$type=="normal"]
        print(paste("Samples with mature: ", nrow(mature.pri)))
        print(paste("Samples with normal mature: ", nrow(mature.normal)))
        ## fix NA
        min.offset <- 2^min(mature.pri[-which(is.na(mature.pri))])/10
        #min.offset <- 10^-6
        mature.pri <- 2^mature.pri
        mature.pri[is.na(mature.pri)] <- 0
        mature.pri <- log2(mature.pri+min.offset)
        # remove features with no variation
        zero.col <- which(apply(mature.pri, 2, sd)==0)
        if(length(zero.col)>0)
            {
                mature.pri <- mature.pri[,-zero.col]
                print(paste("features removed because of no variation:", length(zero.col)))
            }
        return(mature.pri)
    }


format.pre <- function(pre.raw)
    {
        samples <- format.barcode(rownames(pre.raw))
        ## only primary cancer will be used for downstream analysis
        ## KIRC has two 05: additinoal- new primary
        pre.pri <- pre.raw[samples$types=="primary",,drop=FALSE]
        pre.rec <- pre.raw[samples$types=="recurrent",,drop=FALSE]
        pre.met <- pre.raw[samples$types=="metastatic",,drop=FALSE]
        pre.normal <- pre.raw[samples$types=="normal",,drop=FALSE]
        rownames(pre.pri) <- samples$barcodes[samples$type=="primary"]
        rownames(pre.rec) <- samples$barcodes[samples$type=="recurrent"]
        rownames(pre.met) <- samples$barcodes[samples$type=="metastatic"]
        rownames(pre.normal) <- samples$barcodes[samples$type=="normal"]
        print(paste("Samples with pre: ", nrow(pre.pri)))
        print(paste("Samples with normal pre: ", nrow(pre.normal)))
        ## fix NA
        min.offset <- min(pre.pri[pre.pri!=0])/10
        #min.offset <- 2^min(pre.pri[-which(is.na(pre.pri))])/10
        #min.offset <- 10^-6
        #pre.pri <- 2^pre.pri
        #pre.pri[is.na(pre.pri)] <- 0
        pre.pri <- log2(pre.pri+min.offset)
        # remove features with no variation
        zero.col <- which(apply(pre.pri, 2, sd)==0)
        if(length(zero.col)>0)
            {
                pre.pri <- pre.pri[,-zero.col]
                print(paste("features removed because of no variation:", length(zero.col)))
            }
        return(pre.pri)
    }
