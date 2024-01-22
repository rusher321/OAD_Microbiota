#----------------------------- sfig3-----------------------------------# 

# sfig3a 
phe <- readRDS("../inputdata/phenotype_inf.Rds")
pro_tax_com <- readRDS("../inputdata/species_pro.Rds")
factor_study <- c("GuYY_2017_Acarbose", "ZhangXY_2021_Acarbose", "ZhangYF_2020_Berberine","WuH_2017_Metformin",
                  "RenHH_2023_Metformin", "ZhangXY_2021_Vlidagliptin","GuYY_2017_Glipizide",  "ZhangYF_2020_Placebo")

distD1 <- as.matrix(vegan::vegdist(vegan::decostand(pro_tax_com, method = "hellinger"), method = "euclidean", binary = 1))
distD2 <- philentropy::JSD(as.matrix(pro_tax_com))
rownames(distD2) <- colnames(distD2) <- rownames(pro_tax_com)

distD3 <- 1 - cor(t(pro_tax_com), method = "s")
distD4 <- as.matrix( vegan::vegdist(pro_tax_com, method = "bray"))

baseid <- rownames(phe[phe$Time == "D0", ])
basephe <- phe[baseid, ]
## all baseline 
distVar <- c("hell", "jds", "spear", "bray")
distList <- list(distD1, distD2, distD3, distD4)

figlist <- list()
for(i in 1:8){
  phe$study <- paste0(phe$Cohort, "_", phe$Drug)
  sub_id <- rownames(phe[phe$study == factor_study[i], ])
  pcoa_sub <- pcoa_base <- dudi.pco(as.dist(distD1[sub_id, sub_id]), scannf = F, nf = 2)
  var1 <- round((pcoa_sub$eig[1]/sum(pcoa_sub$eig)) * 100, 2)
  var2 <- round((pcoa_sub$eig[2]/sum(pcoa_sub$eig)) * 100, 2)
  pcoa_sub_d <- as.data.frame(pcoa_sub$li)
  pcoa_sub_d$Time <- phe[sub_id, "Time"]
  
  figlist[[i]] <- ggplot(pcoa_sub_d, aes(x=A1, y=A2, shape = Time)) +
    geom_point(color = color_study[i]) +
    xlab(paste("PCoA1: ", var1, "%")) +
    ylab(paste("PCoA2: ", var2, "%")) +
    theme_ipsum(base_family = "Helvetica", plot_title_size = 10)+
    ggtitle(factor_study[i])+
    stat_ellipse(level = 0.8)
  
}

all <- cowplot::plot_grid(plotlist = figlist, nrow = 2)

ggsave("../result/Figure1_sfigure1_3/sfig3a_all_prevspost.pcoa.pdf", plot = all, device = "pdf", width = 15, height = 6)


# sfig3b  venn 

com_clr <- readRDS("../inputdata/clr_compare.Rds")
com_clr$enrich3 <- ifelse(com_clr$enrich2 == "None", "none", com_clr$enrich3)
col_clr_enrich <- com_clr[com_clr$enrich3 == "Post", ]
plot_list <- lapply(split(col_clr_enrich, f = col_clr_enrich$group), function(x){x$tax})

metfor_plot_up <- list(plot_list$WuH_2017_Metformin_D60, plot_list$WuH_2017_Metformin_D120, plot_list$RenHH_2023_Metformin_D90, plot_list$ZhangYF_2020_Berberine_D84) 
names(metfor_plot_up) <- c("WuH_2017_Metformin_D60", "WuH_2017_Metformin_D120", "RenHH_2023_Metformin_D90", "ZhangYF_2020_Berberine_D84")

plot_list[["Metformin"]] <- names(which(table(c(plot_list$WuH_2017_Metformin_D60, plot_list$WuH_2017_Metformin_D120, plot_list$RenHH_2023_Metformin_D90)) >=2))
plot_list[["WuH_2017_Metformin_D60"]] <- NULL
plot_list[["WuH_2017_Metformin_D120"]] <- NULL
plot_list[["ZhangYF_2020_Placebo_D84"]] <- NULL
plot_list[["RenHH_2023_Metformin_D90"]] <- NULL
plot_list[["Acarbose"]] <- intersect(plot_list$GuYY_2017_Acarbose_D90, plot_list$ZhangXY_2021_Acarbose_D168)
plot_list[["GuYY_2017_Acarbose_D90"]] <- NULL
plot_list[["ZhangXY_2021_Acarbose_D168"]] <- NULL

library(venn)
library(ggplot2)
library(ggpolypath)
names(plot_list) <- c("Glipizide", "Vildagliptin", "Berberine", "Metformin", "Acarbose")
p1 <- venn(plot_list, ggplot= T, zcolor = brewer.pal(n = 8, name = "Dark2"), ilcs = 2)+theme_classic()+ggtitle("Species increased in abundance\nfollowing OAD treatment")

col_clr_enrich <- com_clr[com_clr$enrich3 == "Pre", ]
plot_list <- lapply(split(col_clr_enrich, f = col_clr_enrich$group), function(x){x$tax})
metfor_plot_down <- list(plot_list$WuH_2017_Metformin_D60, plot_list$WuH_2017_Metformin_D120, plot_list$RenHH_2023_Metformin_D90, plot_list$ZhangYF_2020_Berberine_D84)
names(metfor_plot_down) <- c("WuH_2017_Metformin_D60", "WuH_2017_Metformin_D120", "RenHH_2023_Metformin_D90", "ZhangYF_2020_Berberine_D84")

tmp_acar <- plot_list$GuYY_2017_Acarbose_D90
plot_list[[1]] <- "a"
names(plot_list)[1] <- "GuYY_2017_Glipizide_D90"
plot_list[["Metformin"]] <- names(which(table(c(plot_list$WuH_2017_Metformin_D60, plot_list$WuH_2017_Metformin_D120, plot_list$YangFM_2021_Metformin_D90)) >=2))
plot_list[["WuH_2017_Metformin_D60"]] <- NULL
plot_list[["WuH_2017_Metformin_D120"]] <- NULL
plot_list[["ZhangYF_2020_Placebo_D84"]] <- NULL
plot_list[["RenHH_2023_Metformin_D90"]] <- NULL
plot_list[["Acarbose"]] <- intersect(tmp_acar, plot_list$ZhangXY_2021_Acarbose_D168)
plot_list[["GuYY_2017_Acarbose_D90"]] <- NULL
plot_list[["ZhangXY_2021_Acarbose_D168"]] <- NULL

names(plot_list) <- c("Glipizide", "Vildagliptin", "Berberine", "Metformin", "Acarbose")

p2 <- venn(plot_list, ggplot= T, zcolor = brewer.pal(n = 8, name = "Dark2"), ilcs = 2)+theme_classic()+ggtitle("Species decreased in abundance\nfollowing OAD treatment")

library(patchwork)
out <- p1+p2
ggsave(plot = out, filename = "../result/Figure1_sfigure1_3/sfig3b.pdf", device = "pdf", width = 10, height = 6)

# sfig3c 

p1_metfor <- venn(metfor_plot_up, zcolor = brewer.pal(n = 6, name = "RdPu")[4:6], ggplot= T, ilcs = 2)+theme_classic()
p2_metfor <- venn(metfor_plot_down, zcolor = brewer.pal(n = 6, name = "RdPu")[4:6], ggplot= T, ilcs = 2)+theme_classic()

spe_clr <- readRDS("../inputdata/species_clr_pro.Rds")
comphe <- readRDS("../inputdata/phenotype_inf.Rds")
comphe_sub <- comphe[comphe$Drug == "Metformin" | comphe$Drug == "Berberine", ]
spe_clr_sub <- spe_clr[rownames(comphe_sub), ]
qdat <- as.data.frame(cbind(spe_clr_sub, comphe_sub))
qdat$Time <- factor(qdat$Time , levels = c("D0", "D60", "D90", "D120"))

box_plot <- function(qdat, x){
  
  p1 <- ggplot(qdat, aes_string(x = "Time", y = x, color = "Cohort")) +
    geom_rect(xmin = 0.4, xmax = 2.5,
              ymin = -Inf, ymax = Inf,
              fill ='white',
              inherit.aes = F)+
    #geom_jitter()+
    geom_line(aes(group=PID), color = "grey", alpha = 0.6)+
    geom_boxplot(alpha = 0.5) +
    scale_color_manual(values = c("#7570B3", "#1F78B4", "#D95F02"))+
    stat_compare_means(comparisons = list(1,2))+
    theme_bw() + 
    xlab("") +
    ylab(paste0(x, " \n(CLR-transformed RA)"))+scale_y_continuous(expand = c(0.1,0.1))+
    theme(panel.grid=element_blank(),
          legend.position = "none",
          axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+facet_wrap(.~Cohort, scales = "free", nrow = 1)
  return(p1)
}

up_list <- Reduce(intersect, metfor_plot_up)

down_list <- Reduce(intersect, metfor_plot_down)

test <- lapply(c(up_list, down_list), box_plot, qdat = qdat)
library(cowplot)
test2 <- plot_grid(plotlist = test, nrow = 2)
#### output 
library(patchwork)
out <- p1_metfor+p2_metfor
ggsave(filename = "../result/Figure1_sfigure1_3/sFig3_left.pdf", plot = out, device = "pdf", width = 8, height = 8)
ggsave(filename = "../result/Figure1_sfigure1_3/sfig3_right.pdf", plot = test2, device = "pdf", width = 4, height = 6)
