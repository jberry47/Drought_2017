library(raster)
library(ggplot2)
library(deldir)
library(vegan)
library(plyr)
library(gridExtra)
library(reshape2)
library(FactoMineR)
library(factoextra)
library(chngpt)
library(stringr)
library(gstat)
library(sp)

#*************************************************************************************************
# Reading and formatting data
#*************************************************************************************************
setwd("/media/jberry/Extra Drive 1/Danforth/16S_Bacteria/Field/2017_Drought")

## Reading, cleaning, and writing the OTU table
# dat <- read.table("Drought2017.nonchimeras_OTU_table_ee05-100-095.txt",sep = "\t", comment.char = "", header = T)
# dat <- t(dat)
# dat[1:5,1:5]
# dat_colnames <- dat[1,]
# dat_rownames <- rownames(dat)[-1]
# dat <- dat[-1,]
# colnames(dat) <- dat_colnames
# dat <- apply(dat,2, as.numeric)
# rownames(dat) <- dat_rownames
# R <- rowSums(dat)
# dat <- dat[R > 10000,] #removes bad samples
# C <- colSums(dat)
# dat <- dat[,C < 50000 &  C > 100] #removes bad OTUs
# mR <- max(R)
# dat <- data.frame(t(apply(dat,1,function(i)floor(mR*(i/sum(i))))))
# dat[1:5,1:5]
# write.csv(dat,"Drought2017_rz-root-soil_OTUtab_cleaned.csv",quote = F)

## Reading other field data, joining, and writing
# dat <- read.table("d2819df383166d1751b1cb504ac64590-schachtman_2017_soil.tsv",sep="\t",header = T,stringsAsFactors = F)
# dat <- dat[dat$seq_path != "" & dat$location=="Scottsbluff NE",]
# dat <- dat[dat$treatment %in% c("WW","WS"),]
# dat1 <- dat2 <- dat3 <- dat
# dat1$my_id <- unlist(lapply(strsplit(dat1$barcode,""),function(i){i[length(i)-6] <- "R"; paste(i,collapse = "")}))
# dat2$my_id <- unlist(lapply(strsplit(dat2$barcode,""),function(i){i[length(i)-6] <- "Z"; paste(i,collapse = "")}))
# dat3$my_id <- unlist(lapply(strsplit(dat2$barcode,""),function(i){i[length(i)-6] <- "S"; paste(i,collapse = "")}))
# dat <- rbind(dat1,dat2,dat3)
# jgi <- read.table("d2819df383166d1751b1cb504ac64590-JGI_seq_2017.tsv",sep = "\t",header = T,stringsAsFactors = F)
# jgi <- jgi[jgi$seq_path != "",]
# jgi$my_id <- jgi$barcode
# jgi <- jgi[jgi$treatment %in% c("WW","WS"),]
# otu_tab <- read.csv("Drought2017_rz-root-soil_OTUtab_cleaned.csv",stringsAsFactors = F)
# colnames(otu_tab)[1] <- "my_id"
# biom <- openxlsx::read.xlsx("sorghum_systems_download.xlsx",2)
# df <- data.frame("shannon" = diversity(otu_tab[,-1]))
# df$my_id <- otu_tab$my_id
# merged <- join(join(join(join(dat,df,by="my_id"),jgi[,c("tissuetype","my_id")],by="my_id"),biom[,-c(1,2,4:8,ncol(biom))],by="unique_id"),otu_tab,by="my_id")
# data.table::fwrite(merged,"Drought2017_rz-root-soil_joined_cleaned.csv")

merged <- data.frame(data.table::fread("Drought2017_rz-root-soil_joined_cleaned.csv"),stringsAsFactors = F)
merged <- merged[!is.na(merged$tissuetype == ""),]
merged[1:5,1:20]

#*************************************************************************************************
# Alpha diversity spatial effects
#*************************************************************************************************
thresh_test <- read.csv("manual_testing_wsOnly_sweetenergyOnly_dryStalk-otu_hinge_norm.csv",stringsAsFactors = F)
thresh_list <- split(thresh_test,thresh_test$OTU)

dat <- merged[merged$sample_date == "8/1/17",]
dat <- dat[dat$treatment %in% c("WW","WS"),]
dat$coord_x <- ceiling(as.numeric(substr(dat$unique_id,5,6))/2)
dict <- data.frame("treatment"=sapply(c(0,0,1,1,1,1,0,0,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0,1,1,0,0,0,0,1,1),function(i) if(i==1){"WS"}else{"WW"}),"rep"=rep(1:8,each=4),"row"=rep(c(0,1),16),"coord_y"=(1:43)[-c(3,8,13,16,19,24,27,30,33,36,41)],stringsAsFactors = F)
dat$row <- as.numeric(substr(dat$unique_id,5,6))%%2
dat <- join(dict,dat,by=c("treatment","rep","row"),type="left")
head(dat)

trait <- "OTU5109"
sub <- dat[!is.na(dat[,trait]),]
design_effects <- aggregate(data=sub,as.formula(paste0(trait,"~treatment+geno")),function(i) mean(as.numeric(unlist(i))))
design_sd <- aggregate(data=sub,as.formula(paste0(trait,"~treatment+geno")),function(i) sd(as.numeric(unlist(i))))
sub$trait_n <- apply(sub[,c("treatment","geno",trait)],1,function(i){(as.numeric(i[3])-design_effects[design_effects$treatment == i[1] & design_effects$geno ==i[2],trait])/design_sd[design_sd$treatment == i[1] & design_sd$geno ==i[2],trait]})
coordinates(sub) <- ~coord_x+coord_y
vgm1 <- variogram(trait_n~1,sub)
fit <- fit.variogram(vgm1,model=vgm(min(as.numeric(unlist(sub[,trait]@data))),"Sph"))
sub.grid <- data.frame("x"=rep(1:12,each=43),"y"=rep(1:43,12))
coordinates(sub.grid) <- ~x+y
kriged <- gstat::krige(trait_n~1,sub,sub.grid,fit,maxdist=3)
out <- setNames(data.frame(kriged)[,1:3],c("x","y","pred"))
out <- join(out,setNames(dict,c("treatment","rep","row","y")),by="y")
out[is.na(out$treatment),"treatment"] <- "Maize"
out[is.na(out$rep),"rep"] <- out[which(is.na(out$rep))-1,"rep"]
out$pred_n <- with(out,(pred-mean(pred,na.rm=T))/sd(pred,na.rm=T))
out$rep <- ordered(out$rep,levels=c(8:1))
sub <- dat[!is.na(dat[,trait]),]

my_sub <- out[interaction(out$x,out$y) %in% interaction(sub$coord_x,sub$coord_y)[sub[sub$treatment=="WS" & sub$type %in% c("Energy","Sweet"),trait]>thresh_list[[trait]]$thresh[1]],]
my_sub <- my_sub[my_sub$treatment =="WS",]

p <- ggplot(out,aes(x,y))+
  facet_grid(rep~.,scales = "free_y")+
  scale_y_continuous(expand = c(0, 0),breaks=1:43)+
  scale_x_continuous(breaks=1:12,limits = c(0.5,14),expand = c(0, 0))+
  geom_tile(aes(fill=pred))+
  geom_tile(data=my_sub,aes(fill=pred),color="black",size=1)+
  geom_text(aes(x=13.2,y=y,label=treatment,color=treatment),show.legend = F)+
  scale_color_manual(values=c("gray20","orange","blue"))+
  scale_fill_gradient2(mid = "gray95",high="darkgreen",limits=c(-2.5,2.5),na.value = "gray60")+
  labs(fill="Z-scale",color="Treatment")+
  ggtitle(trait)+
  theme_light()+
  theme(panel.spacing = unit(0.02, "lines"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=12,color="white"),
    strip.text.y=element_text(size=12,color="white"))
p 
ggsave(paste0("2017_drought_8-1-17only_",trait,"_spatial.png"),width=5.9,height=7.04,plot = p, dpi = 300)


#*************************************************************************************************
# Rarefaction curve
#*************************************************************************************************
sub <- merged[merged$treatment %in% c("WS") & merged$sample_date == "8/1/17" & merged$tissuetype == "Root",]
sub <- sub[!duplicated(sub$my_id),]
y <- sub[,str_detect(colnames(sub),"OTU")]
bad <- which(is.na(rowSums(y)))
y <- y[-bad,]
rownames(y) <- sub$my_id[-bad]
test <- vegan::rarecurve(y,step = 50000)

N <- do.call(rbind,lapply(seq_along(test),function(i){
  n <- names(attr(test[[i]],"Subsample"))
  nn <- length(n)
  n <- n[n != ""]
  data.frame("Sample.Size"=attr(test[[i]], "Subsample"),"Diversity"=test[[i]],"Sample"=rep(n,nn),stringsAsFactors = F)
}))

p <- ggplot(N,aes(Sample.Size,Diversity,group=Sample))+
  geom_line()+
  xlab("Sample Size")+
  ylab("Unique Species")+
  theme_light()+
  #scale_y_continuous(limits=c(0,10000),breaks = seq(0,10000,2500))+
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(axis.text = element_text(size = 14))+
  theme(axis.ticks.length=unit(0.2,"cm"))
p 
ggsave("2017_drought_WSonly_rarefaction-curve.png",p,width = 5.57,height=4.92,dpi=300)


#*************************************************************************************************
# Sample site CAPS
#*************************************************************************************************
sub <- merged[merged$treatment %in% c("WS") & merged$sample_date == "8/1/17",]
sub <- sub[!duplicated(sub$my_id),]
y <- sub[,str_detect(colnames(sub),"OTU")]
bad <- which(is.na(rowSums(y)))
if(length(bad) == 0){y <- y; x <- sub[,!str_detect(colnames(sub),"OTU")]}else{y <- y[-bad,]; x <- sub[-bad,!str_detect(colnames(sub),"OTU")]}

cap_uc <- vegan::capscale(y~tissuetype+Condition(geno),x,distance = "canberra")
cap_uc <- data.frame(summary(cap_uc)$sites)
cap_uc <- cbind(cap_uc,x)
p <- ggplot(cap_uc,aes(eval(parse(text=colnames(cap_uc)[1])),eval(parse(text=colnames(cap_uc)[2]))))+
  geom_point(aes(color=tissuetype))+
  xlab(colnames(cap_uc)[1])+
  ylab(colnames(cap_uc)[2])+
  geom_hline(yintercept=0,linetype="dashed",color="gray20",size=1)+
  geom_vline(xintercept=0,linetype="dashed",color="gray20",size=1)+
  theme_light()+
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(axis.text = element_text(size = 14))+
  theme(axis.ticks.length=unit(0.2,"cm"))
p
ggsave("2017_drought_WSonly_080117only_sample-site_caps.png",p,width = 6.2,height=5.21,dpi=300)

#*************************************************************************************************
# Single soil nutrient CAPS
#*************************************************************************************************
sub <- merged[merged$treatment %in% c("WW","WS") & merged$sample_date == "8/1/17" & merged$tissuetype == "RHZ",]
sub <- sub[!duplicated(sub$my_id),]
y <- sub[,str_detect(colnames(sub),"OTU")]
bad <- which(is.na(rowSums(y)))
if(length(bad) == 0){y <- y; x <- sub[,!str_detect(colnames(sub),"OTU")]}else{y <- y[-bad,]; x <- sub[-bad,!str_detect(colnames(sub),"OTU")]}

cap_uc <- vegan::capscale(y~Potassium._ppm_K+Condition(geno*treatment),x,distance = "canberra")
anova(cap_uc)
cap_uc <- data.frame(summary(cap_uc)$sites)
cap_uc <- cbind(cap_uc,x)
ggplot(cap_uc,aes(eval(parse(text=colnames(cap_uc)[1])),eval(parse(text=colnames(cap_uc)[2]))))+
  geom_point(aes(color=Potassium._ppm_K))+
#  stat_ellipse(aes(fill=Potassium._ppm_K,color=Potassium._ppm_K),geom = "polygon",alpha=0.25)+
  xlab(colnames(cap_uc)[1])+
  ylab(colnames(cap_uc)[2])+
  geom_hline(yintercept=0,linetype="dashed",color="gray20",size=1)+
  geom_vline(xintercept=0,linetype="dashed",color="gray20",size=1)+
  theme_light()+
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(axis.text = element_text(size = 14))+
  theme(axis.ticks.length=unit(0.2,"cm"))
# +
#   guides(color = guide_legend(title = "Potassium._ppm_K"))
#   guides(fill = guide_legend(title = "Potassium._ppm_K"))



#*************************************************************************************************
# All soil nutrients microbiome effect
#*************************************************************************************************
sub <- merged[merged$treatment %in% c("WW","WS") & merged$sample_date == "8/1/17" & merged$tissuetype == "Soil",]
sub <- sub[!duplicated(sub$my_id),]
y <- sub[,str_detect(colnames(sub),"OTU")]
bad <- which(is.na(rowSums(y)))
if(length(bad) == 0){y <- y; x <- sub[,!str_detect(colnames(sub),"OTU")]}else{y <- y[-bad,]; x <- sub[-bad,!str_detect(colnames(sub),"OTU")]}

library(foreach)
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl,cores = 8)
measures <- colnames(x)[-c(1:6,7,8,9,10,11,12,14,16,27:40)]
soil_effects <- foreach(i = measures,.combine = "rbind",.packages=c("vegan","stringr"),.export = c("y")) %dopar% {
  ind_fmla <- as.formula(paste0("y ~ ",i,"+Condition(geno*treatment)",collapse = ""))
  cap_uc <- vegan::capscale(ind_fmla,x,distance = "canberra")
  my_test <- anova(cap_uc)
  c(i,my_test$`Pr(>F)`[1]) 
}
stopCluster(cl)

soil_effects1 <- data.frame(soil_effects,stringsAsFactors = F)
rownames(soil_effects1) <- NULL
soil_effects1$X2 <- as.numeric(soil_effects1$X2)
colnames(soil_effects1) <- c("Source","Effect")
head(soil_effects1)

p <- corrplot::corrplot(cor(sapply(x[,measures[measures != "sample_date"]],as.numeric)),type = "upper",order = "hclust",tl.col = "black",cl.pos = "b")
p[upper.tri(p)] <- NA
p_melt <- na.omit(melt(p))
my_labs <- data.frame("Var1"=unique(p_melt$Var2),"Var2"=unique(p_melt$Var2)) 
p1 <- ggplot(p_melt,aes(Var1,Var2))+
  geom_point(aes(color=value),size=10)+
  geom_text(data=my_labs,aes(label=Var2),check_overlap = T,nudge_x = -0.75,hjust=1,size=4.25)+
  expand_limits(x=-4)+
  scale_color_gradient2(limits=c(-1,1))+
  theme_void()+
  theme(plot.margin=unit(c(0,0,0,0.48), "cm"))+
  theme(legend.position="top")+
  guides(color = guide_colourbar(barwidth = 15,title = NULL))+
  theme(axis.text.x = element_text(size = 12,color="black"),
    axis.text.y = element_blank(),
    axis.title= element_blank())+ 
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
  

soil_effects1$Source <- ordered(soil_effects1$Source,levels=rownames(p))
p2 <- ggplot(na.omit(soil_effects1),aes(Source,-log10(Effect)))+
  geom_bar(stat = "identity")+
  geom_hline(yintercept = -log10(0.05),linetype="dashed",color="gray20",size=1.5)+
#  ggtitle("Microbiome Effect")+
  theme_light()+
  theme(axis.text = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title= element_text(size = 18))+ 
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+
  theme(plot.margin=unit(c(1.25,0.25,4.8,0.48), "cm"))+
  coord_flip()

a <- grid.arrange(p1,p2,nrow=1,widths=c(1,0.75))
ggsave("2017_drought_soil-nutrient-microbiome_effect_08012017_norm.png",a,width = 11.33,height=7.50,dpi=300)


#*************************************************************************************************
# In WS, OTU association to genotype accounting for soil effects
#*************************************************************************************************
sub <- merged[merged$treatment %in% c("WS") & merged$sample_date == "8/1/17" & merged$tissuetype == "RHZ",]
sub <- sub[!duplicated(sub$my_id),]
y <- sub[,str_detect(colnames(sub),"OTU")]
bad <- which(is.na(rowSums(y)))
if(length(bad) == 0){y <- y; x <- sub[,!str_detect(colnames(sub),"OTU")]}else{y <- y[-bad,]; x <- sub[-bad,!str_detect(colnames(sub),"OTU")]}

sub1 <- apply(x[,-c(1:6,7,8,9,10,11,12,14,16,27:40)],MARGIN = 2,as.numeric)
sub1[1:5,1:5]
pca = PCA(sub1, graph = F)
pca_df = data.frame("my_id" = sub$my_id[-bad],
  "PC1" = pca$ind$coord[, 1],
  "PC2" = pca$ind$coord[, 2])
head(pca_df)
varexp = signif(c(pca$eig[1,2], pca$eig[2,2]), 4)
merged1 <- merged[merged$my_id %in% pca_df$my_id,]
merged1 <- join(merged1,pca_df,by="my_id")

#plot PCA
p1 = ggplot(merged1, aes(PC1, PC2)) +
  geom_point(aes(color=tissuetype)) + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_light() +
  xlab(paste("PC1 (", varexp[1], "%)", sep = "")) +
  ylab(paste("PC2 (", varexp[2], "%)", sep = ""))
p1
ggsave("2017_drought_rhz-nutrient-microbiome_effect_08012017_WSonly_PCA.png",p1,width = 5.7,height=4.7,dpi=300)

write.csv(signif(pca$var$contrib,4),"ws_pca_loadings.csv",quote = F)
skree <- data.frame(pca$eig)
skree$PC <- paste0("PC",1:12)
skree$PC <- ordered(skree$PC, levels=paste0("PC",1:12))
p <- ggplot(skree,aes(PC,percentage.of.variance))+
  geom_bar(stat = "identity")+
  theme_light()+
  ylab("Percentage of Variance")+
  theme(axis.text = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title= element_text(size = 18))+ 
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
p 
ggsave("2017_drought_soil-rhz-microbiome_effect_08012017_WSonly_PCA_skree.png",p,width=9,height=4.5,dpi=300)

library(NBZIMM)
data.table::fwrite(merged1,"Drought2017-rhz_root_soil-merged1.csv")
sub <- merged1
otu <- sub[,str_detect(colnames(sub),"OTU")]
N = rowSums(otu)
Genotype <- factor(sub$geno,levels=c("G1","E17","E14","S2","S10","E6","E9","E11","S5","E4"))
Barcodes <- sub$my_id
PC1 <- sub$PC1
PC2 <- sub$PC2
unique_id <- sub$unique_id
otu1 <- otu[,3]
f_microbe_abb = glmm.nb(otu1 ~ Genotype + PC1, random = ~ 1 | unique_id)
out <- anova(f_microbe_abb)

mean(otu[,"OTU102063"],na.rm=T)

library(foreach)
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl,cores = 8)
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
head(test_otus)


f_geno <- read.csv("manual_testing_ws_genotype_effect.csv",stringsAsFactors = F)
head(f_geno)
write.csv(with(f_geno[f_geno$Componant == "PC2",],OTU[p.value < 0.05] ),"pc2_otus.csv",row.names = F,quote = F)

aggregate(data=merged1,unique_id~geno+treatment+sample_date,function(i)length(unique(i)))



#*************************************************************************************************
# In WS, OTU association to biomass with threshold regression
#*************************************************************************************************
sub <- merged[merged$treatment %in% c("WS") & merged$sample_date == "8/1/17" & merged$type %in% c("Sweet","Energy"),]
y <- sub[,str_detect(colnames(sub),"OTU")]
bad <- which(is.na(rowSums(y)))
sub <- sub[-bad,]
#sub <- sub[!is.na(sub$fresh_weight_stalk_g),]
#sub$fresh_dry <- sub$fresh_weight_stalk_g - sub$dry_weight_stalk_g
y <- sub[,str_detect(colnames(sub),"OTU")]

my_otus <- f_geno$OTU[f_geno$Componant == "Genotype" & f_geno$p.value > 0.05]

library(foreach)
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl,cores = 8)

thresh_test <- foreach(i = my_otus,.combine = "rbind",.packages = "chngpt") %do% {
  if(i %in% colnames(sub)){
    fit <- chngptm(formula.1 = dry_weight_stalk_g~1,formula.2 = as.formula(paste0("~",i)), data = sub, type = "hinge",family = "gaussian", est.method = "fastgrid",var.type = "bootstrap",save.boot = T)
    out <- data.frame(coef(summary(fit)))
    out$Source <- rownames(out)
    rownames(out) <- NULL
    out$thresh <- as.numeric(fit$chngpt)
    out$OTU <- i
    out 
  }
}
stopCluster(cl)
write.csv(thresh_test,"manual_testing_wsOnly_sweetenergyOnly_dryStalk-otu_hinge_norm.csv",quote = F,row.names = F)

thresh_test <- read.csv("manual_testing_wsOnly_sweetenergyOnly_dryStalk-otu_hinge_norm.csv",stringsAsFactors = F)
thresh_list <- split(thresh_test,thresh_test$OTU)
thresh_sub <- thresh_list[unlist(lapply(thresh_list,function(i) i$p.value.[2]<0.05))]
thresh_sub[[3]]
as.character(na.omit(as.character(names(thresh_sub))))

#for method = "segmented"
# which_otu <- "OTU1085"
# ggplot(sub[,!duplicated(colnames(sub))],aes_string(which_otu,"fresh_weight_stalk_g"))+
#   geom_point(aes(color=geno),size=2.5)+
#   geom_vline(data=thresh_list[[which_otu]],aes(xintercept = thresh[1]),linetype="dashed")+
#   geom_segment(data=thresh_list[[which_otu]],aes(x=min(sub[,which_otu]),xend=thresh[1],y=est[Source=="(Intercept)"],yend=est[Source=="(Intercept)"]+est[Source==which_otu]*thresh[1]),color="red")+
#   geom_segment(data=thresh_list[[which_otu]],aes(x=thresh[1],xend=max(sub[,which_otu]),y=est[Source=="(Intercept)"]+est[Source==which_otu]*thresh[1],yend=est[Source=="(Intercept)"]+est[Source==which_otu]*max(sub[,which_otu])+est[Source==paste0("(",which_otu,"-chngpt)+")]*(max(sub[,which_otu])-thresh[1])),color="red")+
#   theme_light()+
#   theme(axis.text = element_text(size = 14),
#     axis.title= element_text(size = 18))+
#   theme(aspect.ratio=1)

#for method = "hinge"
which_otu <- "OTU282"
p <- ggplot(sub[,!duplicated(colnames(sub))],aes_string(which_otu,"dry_weight_stalk_g"))+
 geom_point(aes(color=geno),size=2.5)+
 geom_vline(data=thresh_list[[which_otu]],aes(xintercept = thresh[1]),linetype="dashed")+
 geom_segment(data=thresh_list[[which_otu]],aes(x=min(sub[,which_otu]),xend=thresh[1],y=est[Source=="(Intercept)"],yend=est[Source=="(Intercept)"]),color="red")+
 geom_segment(data=thresh_list[[which_otu]],aes(x=thresh[1],xend=max(sub[,which_otu]),y=est[Source=="(Intercept)"],yend=est[Source=="(Intercept)"]+est[Source==paste0("(",which_otu,"-chngpt)+")]*(max(sub[,which_otu])-thresh[1])),color="red")+
  theme_light()+
  theme(axis.text = element_text(size = 14),
    axis.title= element_text(size = 18))
ggsave("2017_drought_wsOnly_energyOnly_OTU866-height_hinge.png",width=6.55,height=5.38,plot = p, dpi = 300)

pdf("2017_drought_wsOnly_energyOnly_height-otu_hinge_hits.pdf")
for(i in as.character(na.omit(as.character(names(thresh_sub))))){
  which_otu <- i
  print(
    ggplot(sub[,!duplicated(colnames(sub))],aes_string(which_otu,"plant_height_cm"))+
      geom_point(aes(color=geno),size=2.5)+
      geom_vline(data=thresh_list[[which_otu]],aes(xintercept = thresh[1]),linetype="dashed")+
      geom_segment(data=thresh_list[[which_otu]],aes(x=min(sub[,which_otu]),xend=thresh[1],y=est[Source=="(Intercept)"],yend=est[Source=="(Intercept)"]),color="red")+
      geom_segment(data=thresh_list[[which_otu]],aes(x=thresh[1],xend=max(sub[,which_otu]),y=est[Source=="(Intercept)"],yend=est[Source=="(Intercept)"]+est[Source==paste0("(",which_otu,"-chngpt)+")]*(max(sub[,which_otu])-thresh[1])),color="red")+
      theme_light()+
      theme(axis.text = element_text(size = 14),
        axis.title= element_text(size = 18))+
      theme(aspect.ratio=1)
  )
}
dev.off()  


#*************************************************************************************************
# In WS, OTU association to biomass with zinbmm
#*************************************************************************************************
sub <- merged[merged$treatment %in% c("WS") & merged$sample_date == "8/1/17" & merged$type %in% c("Energy","Sweet"),]
y <- sub[,str_detect(colnames(sub),"OTU")]
bad <- which(is.na(rowSums(y)))
sub <- sub[-bad,]
y <- sub[,str_detect(colnames(sub),"OTU")]

Height <- sub$plant_height_cm
Biomass <- sub$dry_weight_stalk_g
Barcodes <- sub$my_id
my_otus <- f_geno$OTU[f_geno$Componant == "Genotype" & f_geno$p.value > 0.05]

library(foreach)
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl,cores = 8)
test_otus_height <- foreach(i = my_otus,.combine = "rbind",.packages = "NBZIMM") %dopar% {
  if(i %in% colnames(y)){
    otu1 <- y[,i]
    if(mean(otu1==0)<0.85){
      tryCatch({
        if(all(na.omit(otu1)>0)){
          f_microbe_abb = glmm.nb(otu1 ~ Biomass, random = ~ 1 | Barcodes)
          out <- anova(f_microbe_abb)
          out$Componant <- rownames(out)
          rownames(out) <- NULL
          out$OTU <- i
          out$model <- "NB"
          out
        }else{
          f_microbe_abb = glmm.zinb(otu1 ~ Biomass, random = ~ 1 | Barcodes)
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
}
stopCluster(cl)
write.csv(test_otus_height,"manual_testing_wsOnly_sweetenergyOnly_dryStalk-otu_zinbmm_norm.csv",quote = F,row.names = F)


#*************************************************************************************************
# Intersection of model hits and characterization
#*************************************************************************************************
nb_mod <- read.csv("manual_testing_wsOnly_sweetenergyOnly_dryStalk-otu_zinbmm_norm.csv",stringsAsFactors = F)
head(nb_mod)
write.csv(with(nb_mod[nb_mod$Componant == "Height",],OTU[p.value < 0.05] ),"nb_height_otus.csv",row.names = F,quote = F)

thresh_mod <- read.csv("manual_testing_wsOnly_sweetenergyOnly_dryStalk-otu_hinge_norm.csv",stringsAsFactors = F)
head(thresh_mod)
write.csv(as.character(na.omit(as.character(names(thresh_sub)))),"thresh_height_otus.csv",row.names = F,quote = F)


strong_otus <- intersect(
  as.character(unique(thresh_mod$OTU[thresh_mod$p.value. < 0.05 & thresh_mod$Source != "(Intercept)" & !(thresh_mod$OTU %in% f_geno$OTU[f_geno$p.value< 0.05 & f_geno$Componant != "(Intercept)"])])),
  as.character(unique(nb_mod$OTU[nb_mod$p.value < 0.05 & nb_mod$Componant == "Biomass" & !(nb_mod$OTU %in% f_geno$OTU[f_geno$p.value< 0.05 & f_geno$Componant != "(Intercept)"])]))
)

strong_promoters <- strong_otus[sapply(strong_otus,function(i) thresh_list[[i]]$est[2]>0)]
strong_stunters <- strong_otus[sapply(strong_otus,function(i) thresh_list[[i]]$est[2]<0)]

fa <- seqinr::read.fasta("zotus_rename.fa",forceDNAtolower = F,as.string = T)
strong_promoters_fa <- fa[strong_promoters]
seqinr::write.fasta(strong_promoters_fa,names(strong_promoters_fa),file.out ="drought_strong_promoters_dryStalk.fa")

df <- data.frame(stringsAsFactors = F,
                 "efficacy"=sapply(strong_promoters,function(i)thresh_list[[i]]$thresh[1]),
                 "potency"=sapply(strong_promoters,function(i)thresh_list[[i]]$est[2]),
                 "mls"=sapply(strong_promoters,function(i){percs <- merged[,i]/rowSums(y)*100;merged$my_id[which(percs==max(percs,na.rm = T))]}),
                 "mls_perc"=sapply(strong_promoters,function(i){percs <- merged[,i]/rowSums(y)*100;percs[which(percs==max(percs,na.rm = T))]}),
                 "seq"=sapply(strong_promoters,function(i)as.character(fa[[i]]))
)
head(df)
write.csv(df,"drought_wsOnly_sweetenergyOnly_dryStalk-otu_strong-promoters_summary.csv",quote = F)

pdf("2017_drought_wsOnly_sweetenergyOnly_dryStalk-otu_hinge_strong-promoters.pdf")
for(i in strong_promoters){
  which_otu <- i
  print(
    ggplot(sub[,!duplicated(colnames(sub))],aes_string(which_otu,"dry_weight_stalk_g"))+
      geom_point(size=3,color="black")+
      geom_point(aes(color=geno),size=2.5)+
      geom_vline(data=thresh_list[[which_otu]],aes(xintercept = thresh[1]),linetype="dashed")+
      geom_segment(data=thresh_list[[which_otu]],aes(x=min(sub[,which_otu]),xend=thresh[1],y=est[Source=="(Intercept)"],yend=est[Source=="(Intercept)"]),color="red")+
      geom_segment(data=thresh_list[[which_otu]],aes(x=thresh[1],xend=max(sub[,which_otu]),y=est[Source=="(Intercept)"],yend=est[Source=="(Intercept)"]+est[Source==paste0("(",which_otu,"-chngpt)+")]*(max(sub[,which_otu])-thresh[1])),color="red")+
      theme_light()+
      theme(axis.text = element_text(size = 14),
        axis.title= element_text(size = 18))+
      theme(aspect.ratio=1)
  )
}
dev.off()


#*************************************************************************************************
# Other OTU's in MCL cluster
#*************************************************************************************************
which_otu <- "OTU330"
back <- scan("out.expr72.mci.I60", '', skip = grep(paste0("\\",which_otu,"\\b"),readLines("out.expr72.mci.I60"))-1, nlines = 1, sep = '\t',quiet = T)

tab2 <- table(sapply(strong_stunters,function(i){grep(paste0("\\",i,"\\b"),readLines("out.expr72.mci.I60"))}))
tab1 <- table(sapply(strong_promoters,function(i){grep(paste0("\\",i,"\\b"),readLines("out.expr72.mci.I60"))}))
tab <- table(sapply(thresh_mod$OTU[thresh_mod$p.value. < 0.05 & thresh_mod$Source != "(Intercept)" & !(thresh_mod$OTU %in% f_geno$OTU[f_geno$p.value< 0.05 & f_geno$Componant != "(Intercept)"])],function(i){grep(paste0("\\",i,"\\b"),readLines("out.expr72.mci.I60"))}))

intersect(names(tab2),names(tab1))

max(tab1)
tab1[tab1==2]
tab1

fa <- seqinr::read.fasta("zotus_rename.fa",forceDNAtolower = F,as.string = T)
strong_promoters_fa <- fa[back]
seqinr::write.fasta(strong_promoters_fa,names(strong_promoters_fa),file.out ="drought_MCL-cluster8.fa")


#*************************************************************************************************
# Average Efficacy and Potency within MCL clusters
#*************************************************************************************************
my_list <- list()
con <- file("out.expr72.mci.I60") 
for(i in 1:length(readLines(con))){
  my_list[[i]] <- scan("out.expr72.mci.I60", '', skip = i-1, nlines = 1, sep = '\t',quiet = T)
}

efficacy <- lapply(my_list,function(i){
  effs <- sapply(i,function(j)thresh_list[[j]]$thresh[1])
  effs <- unlist(effs[!is.null(effs)])
  c(mean(effs,na.rm=T),sd(effs,na.rm=T),length(i))
})

potency <- lapply(my_list,function(i){
  c(mean(unlist(sapply(i,function(j)thresh_list[[j]]$est[2])),na.rm=T),
    sd(unlist(sapply(i,function(j)thresh_list[[j]]$est[2])),na.rm=T)
  )
})

df <- cbind(
  setNames(data.frame(do.call(rbind,efficacy)),c("avg_efficacy","sd_efficacy","num_OTUs")),
  setNames(data.frame(do.call(rbind,potency)),c("avg_potency","sd_potency"))
)
df$cluster <- 1:nrow(df)
df$z_eff <- df$avg_efficacy/df$sd_efficacy
df$z_pot <- df$avg_potency/df$sd_potency
df <- na.omit(df[df$z_eff > 1.96 & df$z_pot > 1.96,])
df





#*************************************************************************************************
# WS vs WW OTU association accounting for soil effects
#*************************************************************************************************
sub <- merged[merged$treatment %in% c("WS","WW") & merged$sample_date == "8/1/17",]
y <- sub[,str_detect(colnames(sub),"OTU")]
x <- sub[,colnames(sub) %in% c(colnames(dat),"shannon")]

sub <- x[-c(1:6,7,8,9,10,11,12,14,16,17,21,26,27,28,29)]
pca = PCA(sub, graph = F)
pca_df = data.frame("Samples" = merged$barcode[merged$treatment %in% c("WS","WW") & merged$sample_date == "8/1/17"],
  "PC1" = pca$ind$coord[, 1],
  "PC2" = pca$ind$coord[, 2])
head(pca_df)
varexp = signif(c(pca$eig[1,2], pca$eig[2,2]), 4)
merged1 <- merged[merged$treatment %in% c("WS","WW") & merged$sample_date == "8/1/17",]
merged1$PC1 <- pca_df$PC1
merged1$PC2 <- pca_df$PC2
merged1 <- merged1[,!duplicated(colnames(merged1))]

#plot PCA
p1 = ggplot(merged1, aes(PC1, PC2)) +
  geom_point(aes(color=treatment)) + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_light() +
  xlab(paste("PC1 (", varexp[1], "%)", sep = "")) +
  ylab(paste("PC2 (", varexp[2], "%)", sep = ""))
p1
ggsave("2017_drought_soil-microbiome_effect_08012017_WSonly_PCA.png",p1,width = 5.7,height=4.7,dpi=300)

library(NBZIMM)
sub <- merged1
otu <- sub[,str_detect(colnames(sub),"OTU")]
N = rowSums(otu)
Genotype <- factor(sub$geno,levels=c("G1","E17","E14","S2","S10","E6","E9","E11","S5","E4"))
Barcodes <- sub$my_id
PC1 <- sub$PC1
PC2 <- sub$PC2
unique_id <- sub$unique_id
Drought <- sub$treatment
otu1 <- otu[,3]
f_microbe_abb = glmm.nb(otu1 ~ Drought + PC1 + PC2 + offset(log(N)), random = ~ 1 | unique_id)
out <- anova(f_microbe_abb)

library(foreach)
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl,cores = 8)
test_otus <- foreach(i = colnames(otu),.combine = "rbind",.packages = "NBZIMM") %dopar% {
  otu1 <- otu[,i]
  if(mean(otu1==0,na.rm=T)<0.75){
    tryCatch({
      if(all(na.omit(otu1)>0)){
        f_microbe_abb = glmm.nb(otu1 ~ Drought + PC1 + PC2 + offset(log(N)), random = ~ 1 | unique_id)
        out <- anova(f_microbe_abb)
        out$Componant <- rownames(out)
        rownames(out) <- NULL
        out$OTU <- i
        out$model <- "NB"
        out
      }else{
        f_microbe_abb = glmm.zinb(otu1 ~ Drought + PC1 + PC2 + offset(log(N)), random = ~ 1 | unique_id)
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
write.csv(test_otus,"manual_testing_ws-vs-ww_effect.csv",quote = F,row.names = F)
head(test_otus)


f_geno <- read.csv("manual_testing_ws-vs-ww_effect.csv",stringsAsFactors = F)
head(f_geno)
write.csv(with(f_geno[f_geno$Componant == "Drought",],OTU[p.value < 0.05] ),"drought_otus.csv",row.names = F,quote = F)












