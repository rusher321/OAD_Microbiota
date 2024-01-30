#--------------------- Gee ---------------------------# 
source("../script/fun.R")
library(geeM)
pro <- readRDS("../inputdata/all_pro.rds")
phe <- readRDS("../inputdata/phe_gee.Rds")

phelist <- split(phe, phe$group)

## select species and all the KO 
com_clr <- readRDS("../inputdata/clr_compare.Rds")
load(file = "../inputdata/spe_intersect.Rds")

pro_sub <- pro[, c(downfeature, upfeature, colnames(pro)[136:153])]

reslist <- lapply(phelist,  function(x) {
  phe_tmp <- x[,-2]
  res_tmp <- GEEAna(dataset = pro_sub, metadata = phe_tmp, confounder = c("Age", "Sex", "BMI"),
                    timevar = "Time", scale = F, IDvar = "PID")
  res_tmp
})

colnames(reslist$WuH_2017_Metformin) <- colnames(reslist$GuYY_2017_Acarbose)
qdat <- do.call("rbind", reslist)
study_order <- c("GuYY_2017_Acarbose", "ZhangXY_2021_Acarbose", "ZhangYF_2020_Berberine", 
                 "RenHH_2023_Metformin", "WuH_2017_Metformin", "ZhangXY_2021_Vlidagliptin",
                 "GuYY_2017_Glipizide", "ZhangYF_2020_Placebo")

# compare list 
com_clr_up <- com_clr[com_clr$effect_size < 0 & com_clr$sign_p.value<0.05, ]
com_up <- split(com_clr_up, com_clr_up$group)

com_clr_down <- com_clr[com_clr$effect_size > 0 & com_clr$sign_p.value<0.05, ]
com_down <- split(com_clr_down, com_clr_down$group)

bbr_acar <- intersect(com_up$ZhangYF_2020_Berberine_D84$tax, 
                      com_up$GuYY_2017_Acarbose_D90$tax)
bbr_acar_up <- intersect(upfeature, bbr_acar)
other_up <- setdiff(upfeature, bbr_acar_up)

vild_acar <- intersect(com_down$ZhangXY_2021_Vlidagliptin_D168$tax, 
                      com_down$GuYY_2017_Acarbose_D90$tax)
vild_acar_down <- intersect(downfeature, vild_acar)
other_down <- setdiff(downfeature, vild_acar_down)


qdat_down <- qdat[study_order, 1:24]
qdat_up <- qdat[study_order, 25:90]
qdat_ko <- qdat[study_order, 91:126][, -c(19:24)]
qdat[study_order, 91:126]
pdf(file = "../revision/down_gee.pdf", width = 6, height = 10)
corPlot_tmp(corres = qdat_down, cutoff = 1, adjust = T, tr = T, cluster = T, row_clust = T, col_clust = F)
dev.off()

pdf(file = "../revision/up_gee.pdf", width = 6, height = 10)
corPlot_tmp(corres = qdat_up, cutoff = 1, adjust = T, tr = T, cluster = T, row_clust = T, col_clust = F)
dev.off()

pdf(file = "../revision/ko_gee.pdf", width = 6, height = 10)
corPlot_tmp(corres = qdat_ko, cutoff = 1, adjust = T, tr = T, cluster = F, row_clust = F, col_clust = F)
dev.off()

tax_order <- c(bbr_acar_up, other_up, vild_acar_down, other_down)
tax_order2 <- paste0(rep(tax_order, each = 2), rep(c("wald", "_p.value"), length(tax_order)))
qdat_tax <- qdat[study_order, tax_order2]

pdf(file = "../revision/tax_gee.pdf", width = 6, height = 10)
corPlot_tmp(corres = qdat_tax, cutoff = 1, adjust = T, tr = T, cluster = F, row_clust = F, col_clust = F)
dev.off()

qdat <- rbind(phelist$GuYY_2017_Acarbose, phelist$ZhangXY_2021_Acarbose, phelist$ZhangXY_2021_Vlidagliptin)
qdat_com <- cbind(qdat,  pro[rownames(qdat), ])
qdat_com$group2 <- paste0(qdat_com$group, "_", qdat_com$Time)
qdat_com$Time2 <- ifelse(qdat_com$Time!= "D0", "D90", "D0")

a <- ggplot(qdat_com, aes(x= Bacteroides_xylanisolvens, y= HbA1c))+
  geom_point(aes(color=Time))+facet_grid(.~group, scales = "free")+theme_classic()+
  geom_line(aes(group = PID), color = "lightgrey")

qdat <- rbind(phelist$WuH_2017_Metformin, phelist$RenHH_2023_Metformin , phelist$ZhangYF_2020_Berberine)
qdat_com <- cbind(qdat,  pro[rownames(qdat), ])
qdat_com$group2 <- paste0(qdat_com$group, "_", qdat_com$Time)
qdat_com$Time2 <- ifelse(qdat_com$Time!= "D0", "D90", "D0")

b <- ggplot(qdat_com, aes(x = Clostridium_bartlettii, y= HbA1c))+
  geom_point(aes(color = Time))+facet_grid(.~group, scales = "free")+theme_classic()+
  geom_line(aes(group = PID), color = "lightgrey")

qdat <- rbind(phelist$GuYY_2017_Acarbose, phelist$ZhangXY_2021_Acarbose, phelist$ZhangYF_2020_Berberine)
qdat_com <- cbind(qdat,  pro[rownames(qdat), ])
qdat_com$group2 <- paste0(qdat_com$group, "_", qdat_com$Time)

c <- ggplot(qdat_com, aes(x = `7-beta-hsdh`, y= HbA1c))+
  geom_point(aes(color = Time))+facet_grid(.~group, scales = "free")+theme_classic()+
  geom_line(aes(group = PID), color = "lightgrey")+scale_x_continuous(trans = "log")

d <- ggplot(qdat_com, aes(x = baiI, y= HbA1c))+
  geom_point(aes(color = Time))+facet_grid(.~group, scales = "free")+theme_classic()+
  geom_line(aes(group = PID), color = "lightgrey")+scale_x_continuous(trans = "log")

library(cowplot)
plot_grid(plotlist = list(a,b,c,d), align = "h", ncol = 1)

#-------------------------- Gee v2 -----------------------------------------# 
source("../script/fun.R")
library(geeM)
pro <- readRDS("../inputdata/all_pro.rds")
phe <- readRDS("../inputdata/phe_gee.Rds")

phelist <- split(phe, phe$group)

## select species and all the KO 
com_clr <- readRDS("../inputdata/clr_compare.Rds")
load(file = "../inputdata/spe_intersect.Rds")
pro_sub <- pro[, c(downfeature, upfeature, colnames(pro)[136:153])]
pro_sub$`7-beta-hsdh` -> pro_sub$X7_beta_hsdh

Gee_plot <- function(study, y = "HbA1c", feature = "Bacteroides_xylanisolvens"){
  
  tmpphe <- phelist[[study]]
  qdat <- cbind(tmpphe, pro_sub[rownames(tmpphe), ])
  qdat$Time <- ifelse(qdat$Time == "D0", "Pre", "Post")
  qdat$Time <- factor(qdat$Time, levels = c("Pre", "Post"))
  rmindex <- which(is.na(qdat$HbA1c))
  if(length(rmindex) != 0){
    qdat <- qdat[-rmindex, ]
  }
  p1 <- ggplot(qdat, aes_string(x = "Time", y = y, color = "Time"))+
    geom_boxplot(width = 0.2)+
    geom_point(position=position_dodge(width=0.75), size = 1)+
    geom_line(aes(group=PID), size = 0.3, color = "lightgrey")+
    theme_classic()+theme(legend.position = "none")+
    scale_color_manual(values = c( "#E6A4B4","#8ACDD7"))+xlab("")
  
  p2 <- ggplot(qdat, aes_string(x = "Time", y = feature, color = "Time"))+
    geom_boxplot(width = 0.2)+
    geom_point(position=position_dodge(width=0.75), size =1)+
    geom_line(aes(group=PID), color = "lightgrey")+
    theme_classic()+theme(legend.position = "none")+
    scale_color_manual(values = c( "#E6A4B4","#8ACDD7"))+coord_flip()+ylab("")
  
  p3 <- ggplot(qdat, aes_string(x= feature, y= y))+
    geom_point(aes(color=Time))+theme_classic()+
    scale_color_manual(values = c( "#E6A4B4","#8ACDD7"))+theme(legend.position = "none")
  
  library(aplot)
  qdat <- qdat[, c("HbA1c", feature, "Age", "Sex", "BMI", "PID")]
  qdat <- na.omit(qdat)
  formula_gee <- formula(paste0(feature,"~",y,"+Age+BMI+Sex"))
  geeInd <- geem(formula_gee, id = PID, data = qdat,
                 family = gaussian, corstr = "unstructured")
  res <- summary(geeInd)
  text <- paste0("   GEE, β = ", round(res$beta[2], 2), "\n", 
                 "Rbust SE = ", round(res$se.robust[2],2), ", ",
                 "P = ", res$p[2])
  p3 <- p3+annotate(geom = "text",label = text,  x = -Inf, y = Inf,  hjust = -.2, vjust = 2)+ggtitle(study)
  
  ap <- p3 %>% 
    insert_left(p1, width =.2) %>% 
    insert_bottom(p2, height =.2)
  
  return(ap)
}

plotlist1 <- lapply(c("GuYY_2017_Acarbose", "ZhangXY_2021_Acarbose", "ZhangXY_2021_Vlidagliptin"), 
                    Gee_plot, y = "HbA1c", feature = "Bacteroides_xylanisolvens")

plotlist2 <- lapply(c("WuH_2017_Metformin", "RenHH_2023_Metformin", "ZhangYF_2020_Berberine"),
                    Gee_plot, y = "HbA1c", feature = "Clostridium_bartlettii")

plotlist3 <- lapply(c("GuYY_2017_Acarbose", "ZhangXY_2021_Acarbose", "ZhangYF_2020_Berberine"),
                    Gee_plot, y = "HbA1c", feature = "X7_beta_hsdh")

pdf(file = "gee_scatter.pdf", width = 5, height = 5)
print(plotlist1[[1]])
print(plotlist1[[2]])
print(plotlist1[[3]])
print(plotlist2[[1]])
print(plotlist2[[2]])
print(plotlist2[[3]])
print(plotlist3[[1]])
print(plotlist3[[2]])
print(plotlist3[[3]])
dev.off()


#--------------------------- diversity --------------------------------------# 
pro_spe <- readRDS("../inputdata/species_norm.Rds")[, colnames(pro)[1:135]]
pro_spe_norm <- pro_spe/rowSums(pro_spe)
phe <- readRDS("../inputdata/phenotype_inf.Rds")
phe$group <- paste0(phe$Cohort, "_", phe$Drug)
phe$group2 <- paste0(phe$group, "_", phe$Time)
phelist <- split(phe, phe$group2)
genderate <- function(phe){
  id <- intersect(rownames(phe), rownames(pro_spe_norm))
  res <- diversity_f(pro_spe_norm[id,])
  return(res)
}
diversity <- readRDS("../inputdata/diversity_raw.Rds")

reslist_diversity <- lapply(phelist, genderate)
qdat <- do.call("rbind", reslist_diversity)
rownames(qdat) <- sapply(rownames(qdat), function(x){strsplit(x, split = "[.]")[[1]][2]})
qdat <- cbind(qdat, phe[rownames(qdat), ])
qdat$richness <- diversity[rownames(qdat), "Richness"]
qdat$shannon <- diversity[rownames(qdat), "Shannon"]

std <- function(x) sd(x, na.rm = T)/sqrt(length(x))

levC <- levels(as.factor(qdat$Cohort))
reslist <- list()

for(i in levC){
  met_div <- qdat[qdat$Cohort == i, ]
  reslist[[i]] <- metainf2(met_div)  
}

metainput <- as.data.frame(do.call("rbind", reslist))
namesR <- rownames(metainput) 
metainput <- apply(metainput, 2, function(x){as.numeric(as.character(x))})
rownames(metainput) <- namesR
metainput <- as.data.frame(metainput)
metainput <- metainput[c(1,6,8,4,3,5,7,2),]
rownames(metainput) <- study_name_order[1:8]
metainput$label2 <- rownames(metainput)

library(meta)
library(metafor)
var <- c("richness", "shannon", "Uniqueness_bray", "Uniqueness_kendall",
         "Uniqueness_hellinger", "Uniqueness_JSD", "Uniqueness_spearman")

for(i in var){
  pdf(file = paste0("../revision/", i, "_meta.pdf"), width = 10, height = 6)
  metaplot(metainput, i)
  dev.off()
}

#--------------------------- ko -----------------------------------
source("../script/fun.R")
pro <- readRDS("../inputdata/all_pro.rds")
phe <- readRDS("../inputdata/phenotype_inf.Rds")
featurelist <- c("7-beta-hsdh", "bsh", "baiE", "baiI", 
                 "K00634", "K00929", "K01034", "K01745", 
                 "K17363", "K00864", "K03621", "K00005", 
                 "K05878", "K05879", "K05881")
id <- intersect(rownames(phe), rownames(pro))
tmp <- pro[id, featurelist]
tmp2 <- apply(tmp, 2, function(x){x[x==0] <- min(x[x!=0])/2;x_log <- log(x); scale(x_log)})
rownames(tmp2) <- id

qdat <- cbind(phe[id, ], tmp2)
std <- function(x) sd(x, na.rm = T)/sqrt(length(x))

levC <- levels(as.factor(qdat$Cohort))
reslist <- list()

for(i in levC){
  met_div <- qdat[qdat$Cohort == i, ]
  reslist[[i]] <- metainf_robust(met_div, feature = featurelist)  
}

metainput <- as.data.frame(do.call("rbind", reslist))
namesR <- rownames(metainput) 
metainput <- apply(metainput, 2, function(x){as.numeric(as.character(x))})
rownames(metainput) <- namesR
metainput <- as.data.frame(metainput)
metainput <- metainput[c(1,6,8,4,3,5,7,2),]
rownames(metainput) <- study_name_order[1:8]
metainput$label2 <- rownames(metainput)

library(meta)
library(metafor)

for(i in featurelist){
  pdf(file = paste0("../revision/", i, "_meta.pdf"), width = 10, height = 6)
  metaplot(metainput, i)
  dev.off()
}

#-------------------------- network -------------------------------
pro <- readRDS("../inputdata/all_pro.rds")
pro_spe <- readRDS("../inputdata/species_norm.Rds")[, colnames(pro)[1:135]]
pro_spe_clr <- readRDS("../inputdata/species_clr_pro.Rds")
pro_spe_norm <- pro_spe/rowSums(pro_spe)
phe <- readRDS("../inputdata/phenotype_inf.Rds")
phe$group <- paste0(phe$Cohort, "_", phe$Drug)
phe$group2 <- paste0(phe$group, "_", phe$Time)
phelist <- split(phe, phe$group2)

# raw 
met_d0 <- pro_spe_norm[rownames(phelist$WuH_2017_Metformin_D0), ]
met_d60 <- pro_spe_norm[rownames(phelist$WuH_2017_Metformin_D60), ]
cor_d0 <- cor(met_d0, method = "s")
cor_d60 <- cor(met_d60, method = "s")
diag(cor_d0) <- 0
diag(cor_d60) <- 0
cutoff <- c(0.3, 0.4, 0.5, 0.6, 0.7)
num_d0 <- c()
num_d60 <- c()
for(i in cutoff){
  num_d0 <- c(num_d0, sum(cor_d0 >i, na.rm = T)/2)
  num_d60 <- c(num_d60, sum(cor_d60 >i, na.rm = T)/2)
}

met_d0 <- pro_spe_clr[rownames(phelist$WuH_2017_Metformin_D0), ]
met_d60 <- pro_spe_clr[rownames(phelist$WuH_2017_Metformin_D60), ]
cor_d0 <- cor(met_d0, method = "s")
cor_d60 <- cor(met_d60, method = "s")
diag(cor_d0) <- 0
diag(cor_d60) <- 0
cutoff <- c(0.3, 0.4, 0.5, 0.6, 0.7)
num_d0 <- c()
num_d60 <- c()
for(i in cutoff){
  num_d0 <- c(num_d0, sum(cor_d0 >i, na.rm = T)/2)
  num_d60 <- c(num_d60, sum(cor_d60 >i, na.rm = T)/2)
}

#--------------------------bbr----------------------------------



