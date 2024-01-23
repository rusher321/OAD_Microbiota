# ------------------------------- Fig1 ----------------------------------------# 

## fig1A-B:diversity 

diversity <- readRDS("../inputdata/Fig1_diversity.Rds")

study_name_ord <- c("GuYY_2017_Acarbose_D90", "ZhangXY_2021_Acarbose_D168", 
                    "ZhangYF_2020_Berberine_D84", "WuH_2017_Metformin_D60",
                    "WuH_2017_Metformin_D120",
                    "RenHH_2023_Metformin_D90", 
                    "ZhangXY_2021_Vlidagliptin_D168", 
                    "GuYY_2017_Glipizide_D90",  "ZhangYF_2020_Placebo_D84")

qdat <- diversity
pd <- position_dodge(0.1) ## move them .05 to the left and right
var <- c("Richness", "Shannon")
plist <- list()

for(i in var){
  
  x_var <- paste0(i, "_md")
  sd_var <- paste0(i, "_std")
  p_var <- paste0(i, "_pvalue")
  qdat_tmp <- qdat[, c(x_var, sd_var, p_var, "label")]
  colnames(qdat_tmp)[1:3] <- c("x_var", "sd_var", "p_value") 
  qdat_tmp$shape <- ifelse(qdat_tmp$p_value <0.05, "enrich", "ns")
  qdat_tmp$shape <- factor(qdat_tmp$shape, levels = c("enrich", "ns"))
  qdat_tmp$label <- factor(qdat_tmp$label, levels = rev(study_name_ord))
  
  plist[[i]] <- ggplot(qdat_tmp, aes(x=x_var, y=label)) + 
    geom_errorbar(aes(xmin=x_var-sd_var, xmax=x_var+sd_var), colour="grey", width=.4, position=pd) +
    geom_point(aes(shape = shape, fill = shape),  position=pd, size=6) + # 21 is filled circle
    xlab(paste0("Change of ", i)) +
    ylab("") +
    scale_fill_manual(values = c("#E41A1C",  "#999999")) +
    scale_color_manual(values = c("#E41A1C",  "#999999")) +
    scale_shape_manual(values = c(23, 21))+
    expand_limits(y=0)+                        # Expand y range
    #scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
    theme_classic()+geom_vline(aes(xintercept = 0 ),linetype = "dashed", color= "grey")  
  
}
out <- cowplot::plot_grid(plotlist = plist, nrow = 1, hjust = "h")

ggsave(filename = "../result/Figure1_sfigure1_3/fig1a_b.pdf", plot = out, device = "pdf", width = 9, height = 5)

## fig1C: adonise 

library(hrbrthemes)
library(viridis)
library(ggplot2)

adonise_r2 <- readRDS("../inputdata/Fig1_adonise.Rds")
adonise_r2$label <- factor(adonise_r2$label, rev(study_name_ord))

p1 <- ggplot(adonise_r2, aes(x = label, y= value*100,  color = variable))+
  geom_point(aes(shape=pvalue), size = 5, position = position_dodge(0.3))+
  coord_flip()+ylab("R^2")+xlab("")+ylim(0,10)+
  theme_classic()+
  guides(color=guide_legend(title="Dissimilarty"), size = guide_legend(title = "Pvalue"))+scale_shape_manual(values = c(18, 16))+
  scale_color_ipsum()+
  labs(y="R2 (%)",title="ADONIS")

ggsave(filename = "../result/Figure1_sfigure1_3/fig1c.pdf", plot = p1, device = "pdf", width = 4.5, height = 5.5)

## fig1D: classification matrix 

phe <- readRDS("../inputdata/phenotype_inf.Rds")
pro_tax_com <- readRDS("../inputdata/species_pro.Rds")

library(rmeta) # this package is from my github 

phe_2t <-  phe[phe$Time != "D120", ]
phe_2t$uniq_id <- paste0(phe_2t$PID, "_", phe_2t$Cohort)
phe_2t$uniq_time <- ifelse(phe_2t$Time == "D0", "base", "treat")
phe_2t_match <- matchpairID(phe_2t, "uniq_id", "uniq_time", num = 2)
phe_2t_match$uniq_cohort <- paste0(phe_2t_match$Cohort, "_", phe_2t_match$Drug)
pro_2t_match <- pro_tax_com[rownames(phe_2t_match), ]

phe_3t <- phe[phe$Cohort == "WuH_2017"  & phe$Time != "D60", ]
phe_3t$uniq_id <- paste0(phe_3t$PID, "_", phe_3t$Cohort)
phe_3t$uniq_time <- ifelse(phe_3t$Time == "D0", "base", "treat")
phe_3t$Time <- droplevels(phe_3t$Time)
phe_3t_match <- matchpairID(phe_3t, "PID", "Time", num = 2)
phe_3t_match$uniq_cohort <- paste0(phe_3t_match$Cohort, "_", phe_3t_match$Drug, "_3t")
pro_3t_match <- pro_tax_com[rownames(phe_3t_match), ]

classification_distance = function(pro, phe, cohort){
  ## select the sample  
  if(cohort == "all"){
    pro_sub <- pro[rownames(phe), ]  
    #pro2_sub <- pro2[rownames(phe), ]
  }else{
    pro_sub <- pro[rownames(phe[phe$uniq_cohort == cohort, ]), ]
    phe <- phe[rownames(phe[phe$uniq_cohort == cohort, ]), ]
  }
  base_id <- rownames(phe[phe$Time == "D0", ])
  treat_id <- rownames(phe[phe$Time != "D0", ])
  ## compute the distance 
  distD1 <- as.matrix(vegan::vegdist(vegan::decostand(pro_sub, method = "hellinger"), 
                                     method = "euclidean", binary = 1))
  
  distD2 <- philentropy::JSD(as.matrix(pro_sub))
  rownames(distD2) <- colnames(distD2) <- rownames(pro_sub)
  
  distD3 <- 1 - cor(t(pro_sub), method = "s")
  
  distD4 <- as.matrix(vegan::vegdist(pro_sub, method = "bray"))
  
  distance_res <- list(distD1, distD2, distD3, distD4)
  ## to combine the result  
  ## to init the output 
  result_out <- matrix(NA, nrow = length(base_id), ncol = 4)
  colnames(result_out) <- c("Hellinger", "JSD", "Spearman","Bray")
  rownames(result_out) <- base_id
  
  for(i in 1:4){
    
    ### select the sample 
    dist_tax <- distance_res[[i]]
    dist_tax_select <- dist_tax[base_id, treat_id]
    colnames(dist_tax_select) <- rownames(dist_tax_select)
    sample <- rownames(dist_tax_select)
    result <- NULL
    for(k in 1:length(sample)){
      
      dist = dist_tax_select[grep(sample[k], row.names(dist_tax_select)),]
      intra = dist[grep(sample[k], names(dist))]
      inter = as.vector(dist[grep(sample[k], names(dist),invert = T)])
      if(is.na(intra)==F){
        result = c(result,  1-length(which(inter < intra)))
      }
    }
    result_out[,i] <- result  
  }
  return(result_out)
}


cohort_lev <- names(table(phe_2t_match$uniq_cohort))
cohort_result <- lapply(cohort_lev, function(x){classification_distance(pro = pro_2t_match, 
                                                                        phe = phe_2t_match, cohort = x)})
distance_cohort <- do.call("rbind", cohort_result) #to compute the classification 

cohort_lev <- names(table(phe_3t_match$uniq_cohort))
cohort_result_3t <- lapply(cohort_lev, function(x){classification_distance(pro = pro_3t_match, phe = phe_3t_match, cohort = x)})
distance_cohort_3t <- do.call("rbind", cohort_result_3t) #here to compute the wuhao'3 time point in each cohort 


library(hrbrthemes)
library(viridis)
library(ggplot2)

qdat <- cbind(distance_cohort, phe_2t_match[rownames(distance_cohort), ])
qdat_3t <- cbind(distance_cohort_3t, phe_3t_match[rownames(distance_cohort_3t),])
qdat_3t$uniq_cohort <- paste0(qdat_3t$uniq_cohort, "_D120")
qdat_com <- rbind(qdat, qdat_3t)

qdat_com[,1:4] <- apply(qdat_com[,1:4], 2, function(x){ifelse(x==1, 1, 0)})
class_res <- matrix(NA, nrow = 4, ncol = 9)
rownames(class_res) <- c("Helliger", "JSD", "1-Spearman", "Bray")
colnames(class_res) <- names(table(qdat_com$uniq_cohort))

for(i in 1:4){
  class_res[i, ] <- as.numeric(table(qdat_com[,c(11,i)])[,2]/apply(table(qdat_com[,c(11,i)]), 1, sum))
}

class_res_order <- class_res[, c(1,6,8,7,3,4,5,2,9)]
colnames(class_res_order) <- study_name_ord

library(ComplexHeatmap)
pdf("../result/Figure1_sfigure1_3/fig1d_Classification.distance.pdf", width = 4.5, height = 5)
pheatmap(t(class_res_order), display_numbers  = t(apply(round(class_res_order*100, 2), 2, function(x){paste0(x, "%")})), cluster_rows = F, cluster_cols = F,angle_col = "45", color = rev(colorRampPalette(c("#FBFEF9","#A63446"))(20)), main = "Classification beased on Distance", fontsize_number = 6, fontsize_col = 10)
dev.off()

## Fig1E: compare result 

qdat1 <- readRDS("../inputdata/Fig1_compare.Rds")
taxrank <- readRDS("../inputdata/Fig1_tax_rank.Rds")
taxorder <- readRDS("../inputdata/Fig1_tax_order.Rds")

library(RColorBrewer)
library(circlize)
library(reshape2)

qdat <- dcast(qdat1, tax~group, value.var = "effect_size")
pvalue <- dcast(qdat1, tax~group, value.var = "sign_p.value")
qdat2 <- qdat[,2:10]
rownames(qdat2) <- qdat[, 1]
qdat2[pvalue[,2:10] >0.05] <- 0 

qdat3 <- -qdat2[rownames(taxorder), ]
colnames(qdat3) <- c("Plac_D84", "Glip_D90", "Vild_D168", "Metf_D90",
                     "Metf_D120", "Metf_D60", "BBR_D84", "Acar_D168", "Acar_D90")
split <- taxorder[,1]

circos.clear()

pdf("../result/Figure1_sfigure1_3/Fig1e_species.pdf", width = 7, height = 7)

circos.par(gap.after = c(2, 2, 2, 2, 12))
col_fun1 = colorRamp2(c(-1, 0, 1), c("#003399", "white", "#FF9933"))
circos.heatmap(qdat3, col = col_fun1, split = split, cluster = F, rownames.side = "inside",  
               #bg.border = "grey",
               cell.border = "grey", bg.lwd = 2, bg.lty = 2, show.sector.labels = TRUE, track.height = 0.3)

circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 5) { # the last sector
    cn = rev(colnames(qdat3))
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                1:n - 0.5, cn, 
                cex = 0.5, adj = c(0, 0.5), facing = "inside")
  }
}, bg.border = NA)
library(ComplexHeatmap)
lgd = Legend(title = "Effect size", col_fun = col_fun1)
grid.draw(lgd)
dev.off()


