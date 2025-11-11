TOTAL_MARKERS <- 	176505
EFFECTIVE_MARKERS <- 82836.66
effective_ratio <- EFFECTIVE_MARKERS/TOTAL_MARKERS
library(rMVP)
library(readr)
library(data.table)
library(dplyr)
library(tidyverse)
MVP.Data(fileVCF="SAP.320_PI_fil.vcf",
         filePhe="trait_plasticity.csv",
         sep.phe=",",
         fileKin=TRUE,
         filePC=TRUE,
         priority="memory",
         maxLine=10000,
         out="mvp.RNASEQ.matched320"
)
genotype <- attach.big.matrix("mvp.RNASEQ.matched320.geno.desc")
phenotype <- read.table("mvp.RNASEQ.matched320.phe",head=TRUE)
#colnames(phenotype) <- ifelse(grepl("^PC", colnames(phenotype)),
#                              paste0(colnames(phenotype), "_"),
 #                             colnames(phenotype))
map <- read.table("mvp.RNASEQ.matched320.geno.map" , head = TRUE)
Kinship <- attach.big.matrix("mvp.RNASEQ.matched320.kin.desc")
Covariates_PC <- bigmemory::as.matrix(attach.big.matrix("mvp.RNASEQ.matched320.pc.desc"))
for(x in 1:100){
  phe1=phenotype # make a copy of phenotype
  nline = nrow(phe1)
  phe1[sample(c(1:nline), as.integer(nline*0.1)),2:ncol(phe1)]=NA  # randomly choose 10% phenotype to be NA
  colnames(phe1)=paste0(colnames(phenotype),x)  # rename the phenotype by attaching bootstrap number
  for(i in 2:ncol(phe1)){
    imMVP<-MVP(phe = phe1[,c(1,i)], geno = genotype, map = map, K=Kinship, CV.FarmCPU=Covariates_PC,
               maxLoop = 10, method = "FarmCPU", p.threshold = (0.05/EFFECTIVE_MARKERS),
               threshold = (0.05/effective_ratio),
               file.output = 'pmap.signal')
  }
}
traits=colnames(phenotype[-1])
get.support=function(trait){ # write a function to summarise the occurrence of signals, trait is what i have in the rmvp output filenames, disregarding the number of bootstrap
  files = list.files(pattern = paste0(trait, ".*FarmCPU_signals.csv"))
  print(files)
  if (length(files)>=1){
    signals <-
      files %>%
      map_df(~read.csv(.,skip=1,header=F,colClasses = c("factor","factor","integer","factor","factor","numeric","numeric","numeric", "numeric")))
    header = read.csv(paste0(trait,"1.FarmCPU_signals.csv"))
    colnames(signals)=colnames(header)
    colnames(signals)[9]<-"pvalue"
    signals=signals %>%
      group_by(SNP,CHROM,POS) %>%
      summarise(P=mean(pvalue), RMIP = n()/100) #%>% ## if {{trait}} does not work otherwise change this name to something else
    #separate(SNP, c("CHR","BP"),remove=F)
    write.table(signals, file=paste0("Z", trait, "signals.csv"), quote = F,row.names = F,sep=",")
  }
  else{
    print(paste0("file not found", trait))
  }
}
for(x in traits){get.support(x)}


