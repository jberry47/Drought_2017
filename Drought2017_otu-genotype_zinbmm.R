library(data.table)
library(foreach)
library(doParallel)
library(NBZIMM)
library(stringr)

cl <- makeCluster(20)
registerDoParallel(cl,cores = 20)

merged1 <- data.frame(data.table::fread("Drought2017-rhz_root_soil-merged1.csv"),stringsAsFactors = F)

sub <- merged1
otu <- sub[,str_detect(colnames(sub),"OTU")]
Genotype <- factor(sub$geno,levels=c("G1","E17","E14","S2","S10","E6","E9","E11","S5","E4"))
PC1 <- sub$PC1
PC2 <- sub$PC2
unique_id <- sub$unique_id

test_otus <- foreach(i = colnames(otu),.combine = "rbind",.packages = "NBZIMM") %dopar% {
  otu1 <- otu[,i]
  if(mean(otu1==0,na.rm=T)<0.75){
    tryCatch({
      if(all(otu1>0)){
        f_microbe_abb = glmm.nb(otu1 ~ Genotype + PC1 + PC2, random = ~ 1 | unique_id)
        out <- anova(f_microbe_abb)
        out$Componant <- rownames(out)
        rownames(out) <- NULL
        out$OTU <- i
        out$model <- "NB"
        out
      }else{
        f_microbe_abb = glmm.zinb(otu1 ~ Genotype + PC1 + PC2, random = ~ 1 | unique_id)
        out <- anova(f_microbe_abb)
        out$Componant <- rownames(out)
        rownames(out) <- NULL
        out$OTU <- i
        out$model <- "ZINB"
        out 
      }
    },warning=function(war){},error=function(err){})
  }
}
stopCluster(cl)
write.csv(test_otus,"manual_testing_rhz_ws_genotype_effect_norm.csv",quote = F,row.names = F)