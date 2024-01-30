#-----------------------------Fig3--------------------------------------------#

# fig3a
source("../script/fun.R")
library(ggplot2)
fun_com <- read.csv("../inputdata/bafun.compare.csv",row.names = 1)
fun_com$group <- gsub("YangFM_2021_Metformin_D90", "RenHH_2023_Metformin_D90", fun_com$group)
qdat <- fun_com
library(tidyr)
all_fun <- c("bsh", "baiE", "baiI", "7-beta-hsdh")
effect_size <- spread(qdat[, c(2,5,6)], group, effect_size)
rownames(effect_size) <- effect_size$tax
qdat_com <- -effect_size[all_fun, -1]
qvalue <- spread(qdat[, c(3,5,6)], group, p.adjust)
rownames(qvalue) <- qvalue$tax
qvalue_com <- qvalue[all_fun, -1]
pvalue <- spread(qdat[, c(1,5,6)], group, sign_p.value)
rownames(pvalue) <- pvalue$tax
pvalue_com <- pvalue[all_fun, -1]

com_plot <- qdat_com
com_plot$KO <- rownames(com_plot)

library(reshape2)
com_plot2 <- melt(com_plot)
com_pvalue <- melt(pvalue_com)
com_qvalue <- melt(qvalue_com)

com_plot2$dir <- ifelse(com_plot2$value < 0 , "base", "post") 
com_plot2$variable <- factor(com_plot2$variable, levels = rev(study_name_order))
com_plot2$pvalue <- com_pvalue$value
com_plot2$qvalue <- com_qvalue$value
com_plot2$enrich <- ifelse(com_plot2$qvalue < 0.05 , "*", ifelse(com_plot2$qvalue >0.05 & com_plot2$pvalue <0.05, "#", " "))
com_plot2$KO <- factor(com_plot2$KO, levels = rev(all_fun))
com_plot2 <- na.omit(com_plot2)

#com_plot2$variable <- factor(com_plot2$variable, levels = study_name2)
p_ba <- ggplot(com_plot2, aes(y = KO, x = value, color= dir)) + 
  scale_x_continuous(breaks=seq(-1,1,0.2)) +
  geom_vline(xintercept = c(-0.3, 0.3),linetype="dashed",color="grey50")+
  geom_segment(xend=0,aes(yend=KO),size=5, alpha = 0.6)  +
  geom_vline(xintercept = 0) +
  theme_minimal() +
  geom_text(aes(x = value, y = KO, label = enrich, color = "black"))+
  scale_color_manual(values = c(post="#FF9933", base = "#003399"))+
  theme(panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        strip.text.y = element_blank(),
        axis.title=element_text(size=10,face="bold"),
        strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(color = "black", size= 6),
        axis.text.y = element_text(color = "black", size= 10))+
  xlab("Effect size")+
  ylab("")+facet_grid(.~variable, scales = "free_y", space = "free")

ggsave(plot = p_ba, filename = "../result/Figure3/fig3a_v1.pdf", device = "pdf", width = 11, height = 2)

# fig3d_v2 
com_plot2$value <- -com_plot2$value
com_plot2$shape <- ifelse(com_plot2$value > 0 , "+", ifelse(com_plot2$value < 0 , "-", "o"))
com_plot2$size <- abs(com_plot2$value)
com_plot2$sig <- ifelse(com_plot2$qvalue < 0.05, 1, 0.9)

colt <- c("#4C38CB", "#9191C8", "#DADAEC", 
          "#F0C1C1", "#E28383", "#D44545", "#CD2626")

color <- colorRampPalette(colt)(20)
com_plot2$sig2 <- ifelse(com_plot2$pvalue<0.05, "P<0.05", "P>=0.05")
com_plot2$shape2 <- ifelse(com_plot2$value>0, 2, 6)
qdat <- com_plot2

qdat$KO <- factor(qdat$KO, levels = rev(all_fun))

qdat$group <- factor(qdat$variable, levels = rev(study_name_order))

p_ba2 <- ggplot(qdat, aes( x = KO, y = group))+
  geom_point(aes(size= size, color=sig2, shape = shape))+
  geom_point(aes(size= size+0.2,fill = -value, alpha= sig, color=sig2,  shape = shape),stroke=0)+
  scale_fill_gradientn(
    colours = color,
    guide=guide_colourbar(ticks=T,nbin=50,barheight=.5, label=T,barwidth=10)
  )+
  scale_color_manual(values = c("black","grey"))+
  scale_shape_manual(values = c(24, 25, 1))+
  theme_minimal()+
  nature_theme2+
  theme(legend.position = "top")+scale_size_continuous(range = c(1,3))

ggsave(plot = p_ba2, filename = "../result/Figure3/fig3a_v2.pdf", device = "pdf", width = 4.5, height = 2.8)

# fig3b
source("../script/fun.R")
Meta_fun <- maaslin2_heatmap("../../result/04.bileacid/fig3/gene_meta_maaslin.res.tsv")
colnames(Meta_fun) <- c("BaiE", "BaiI", "Bsh", "7_beta_hsdh")
out <- Meta_fun
max_value <- ceiling(max(out))
min_value <- ceiling(min(out))
range_value <- max(c(abs(max_value),abs(min_value)))
breaks <- seq(-1*range_value, range_value, by = 1)
rownames(out) <- gsub("X", "", rownames(out))
color <- colorRampPalette(c("darkblue", "grey90", "darkred"))
out <- t(out[c("UDCA", "GUDCA", "TUDCA", "DCA", "GDCA", "TDCA", "LCA", "GLCA", "TLCA"), c(4,1,2,3)])


p <-   pheatmap::pheatmap(
  out,
  cellwidth = 5,
  cellheight = 5,
  # changed to 3
  fontsize = 6,
  kmeans_k = NA,
  border = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  legend = TRUE,
  border_color = "grey93",
  color = color(range_value*2),
  breaks = breaks,
  treeheight_row = 0,
  treeheight_col = 0,
  #annotation_row = ann_row,
  display_numbers = matrix(ifelse(
    out > 0.0, "+", ifelse(out < 0.0, "-", "")), nrow(out)),
  silent = TRUE
)
pdf(file = "../result/Figure3/fig3b.pdf", width = 6, height = 4)
p
dev.off()


# fig3c 
com_ba_phe <- readRDS("../inputdata/ba_inf.Rds")
com_ba_phe$Sex <- ifelse(com_ba_phe$Sex == "female", 1, 0)
library(MASS)
library(ggplot2)
library(patchwork)
p1 <- scatterplotWithSCC(com_ba_phe, x = "baiI", y = "DCA", group = "Cohort", PCC = T, adj = c("Age", "Sex"), color = c("#1B9E77",  "#D95F02"))

p2 <- scatterplotWithSCC(com_ba_phe, x = "X7.beta.hsdh", y = "GUDCA", group = "Cohort", PCC = T, adj = c("Age", "Sex"), color = c("#1B9E77",  "#D95F02"))

out1 <- p1+p2

ggsave(filename = "../result/Figure3/fig3c.pdf", plot = out1, device = "pdf", height = 4, width = 10)

# fig3d 
all_fun <- c("K01034", "K00634", "K00929", "K01745", 
             "K17363", "K05878", "K05879", "K05881",
             "K00005", "K00864", "K03621" )

fun_com <- read.csv("../inputdata/all.fun.compare_0807.csv", row.names = 1)
fun_com$group <- gsub("YangFM_2021_Metformin_D90", "RenHH_2023_Metformin_D90", fun_com$group)
qdat <- fun_com
anno <- data.frame( SCFA = c(rep("Butyrate", 3), rep("IMP",2), 
                             rep("Glycerolipid", 6)))
library(tidyr)

effect_size <- spread(qdat[, c(2,5,6)], group, effect_size)
rownames(effect_size) <- effect_size$tax
qdat_com <- -effect_size[all_fun, -1]
qvalue <- spread(qdat[, c(3,5,6)], group, p.adjust)
rownames(qvalue) <- qvalue$tax
qvalue_com <- qvalue[all_fun, -1]
pvalue <- spread(qdat[, c(1,5,6)], group, sign_p.value)
rownames(pvalue) <- pvalue$tax
pvalue_com <- pvalue[all_fun, -1]
anno <- data.frame( SCFA = c(rep("Butyrate", 3), rep("IMP",2), 
                             rep("Glycerolipid", 6)))
rownames(anno) <- all_fun

com_plot <- cbind(qdat_com, anno)
com_plot$KO <- rownames(com_plot)

library(reshape2)
com_plot2 <- melt(com_plot)
com_pvalue <- melt(pvalue_com)
com_qvalue <- melt(qvalue_com)

com_plot2$dir <- ifelse(com_plot2$value < 0 , "base", "post") 
com_plot2$variable <- factor(com_plot2$variable, levels = study_name_order)
com_plot2$pvalue <- com_pvalue$value
com_plot2$qvalue <- com_qvalue$value
com_plot2$enrich <- ifelse(com_plot2$qvalue < 0.05 , "*", ifelse(com_plot2$qvalue >0.05 & com_plot2$pvalue <0.05, "#", " "))
com_plot2$KO <- factor(com_plot2$KO, levels = all_fun)
com_plot2$SCFA <- factor(com_plot2$SCFA, levels = c("Butyrate", "IMP", 
                                                    "Glycerolipid"))
com_plot2 <- na.omit(com_plot2)
#com_plot2$variable <- factor(com_plot2$variable, levels = study_name2)
p_scfa <- ggplot(com_plot2, aes(y = KO, x = value, color= dir)) + 
  scale_x_continuous(breaks=seq(-1,1,0.2)) +
  geom_vline(xintercept = c(-0.3, 0.3),linetype="dashed",color="grey50")+
  geom_segment(xend=0,aes(yend=KO),size=5, alpha = 0.6)  +
  geom_vline(xintercept = 0) +
  theme_minimal() +
  geom_text(aes(x = value, y = KO, label = enrich, color = "black"))+
  scale_color_manual(values = c(post="#FF9933", base = "#003399"))+
  theme(panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        strip.text.y = element_blank(),
        axis.title=element_text(size=10,face="bold"),
        strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(color = "black", size= 6),
        axis.text.y = element_text(color = "black", size= 10))+
  xlab("Effect size")+
  ylab("")+facet_grid(SCFA~variable, scales = "free_y", space = "free")
ggsave(plot = p_scfa, filename = "../result/Figure3/fig3d_v1.pdf", device = "pdf", width = 11, height = 5)

# fig3d_v2 
com_plot2$value <- -com_plot2$value
com_plot2$shape <- ifelse(com_plot2$value > 0 , "+", ifelse(com_plot2$value < 0 , "-", "o"))
com_plot2$size <- abs(com_plot2$value)
com_plot2$sig <- ifelse(com_plot2$qvalue < 0.05, 1, 0.9)

colt <- c("#4C38CB", "#9191C8", "#DADAEC", 
          "#F0C1C1", "#E28383", "#D44545", "#CD2626")

color <- colorRampPalette(colt)(20)
com_plot2$sig2 <- ifelse(com_plot2$pvalue<0.05, "P<0.05", "P>=0.05")
com_plot2$shape2 <- ifelse(com_plot2$value>0, 2, 6)
qdat <- com_plot2

qdat$KO <- factor(qdat$KO, levels = c("K01034", "K00634", "K00929", "K01745", 
                                      "K17363", "K05878", "K05879", "K05881",
                                      "K00005", "K00864", "K03621" ))

qdat$group <- factor(qdat$variable, levels = rev(study_name_order))

p_scfa2 <- ggplot(qdat, aes( x = KO, y = group))+
  geom_point(aes(size= size+0.2, color=sig2, shape = shape))+
  geom_point(aes(size= size+0.2,fill = -value, alpha= sig, color=sig2,  shape = shape),stroke=0)+
  scale_fill_gradientn(
    colours = color,
    guide=guide_colourbar(ticks=T,nbin=50,barheight=.5, label=T,barwidth=10)
  )+
  scale_color_manual(values = c("black","grey"))+
  scale_shape_manual(values = c(24, 25, 1))+
  theme_minimal()+
  nature_theme2+
  theme(legend.position = "top")+scale_size_continuous(range = c(1,3))+theme(axis.text.y = element_blank())+ylab("")

ggsave(plot = p_scfa2, filename = "../result/Figure3/fig3d_v2.pdf", device = "pdf", width = 4.5, height = 4.5)

library(aplot)
out <- p_scfa2 %>% insert_left(p_ba2, width = 0.3)
ggsave(plot = out, filename = "../result/Figure3/fig3d_v3.pdf", device = "pdf", width = 8, height = 4.5)

# fig3e
ba_com <- readRDS("../inputdata/acar_ba.rds")
ba_com_sub <- ba_com[grep("Acar", ba_com$Group) ,]
ba_com_sub2 <- ba_com_sub[ba_com_sub$Gene %in% c("7-beta-hsdh", "baiI", "baiE", "bsh"), ]
ba_com_sub2$shape <- ifelse(ba_com_sub2$Effect_size > 0 , "-", ifelse(ba_com_sub2$Effect_size < 0 , "+", "o"))
ba_com_sub2$size <- abs(ba_com_sub2$Effect_size)
ba_com_sub2$sig <- ifelse(ba_com_sub2$BH.adjusted.Pvalue < 0.05, 1, 0.3)
colt <- c("#4C38CB", "#9191C8", "#DADAEC", 
          "#F0C1C1", "#E28383", "#D44545", "#CD2626")
color <- colorRampPalette(colt)(20)
ba_com_sub2$sig2 <- ifelse(ba_com_sub2$P.value..Wilcox.sign.test.<0.05,"P<0.05","P>=0.05")
ba_com_sub2$shape2 <- ifelse(ba_com_sub2$Effect_size > 0, 2, 6)
qdat <- ba_com_sub2
qdat$group <- factor(qdat$Group, levels = c("ZhangXY_2022_Acarbose_D168","GuYY_2017_Acarbose_D90", "ZhaoLP_2018_Acarbose_D84", "ZhaoLP_2018_Acarbose_D56", "ZhaoLP_2018_Acarbose_D28",  "ZhaoLP_2018_Acar+WTP_D84", "ZhaoLP_2018_Acar+WTP_D56", "ZhaoLP_2018_Acar+WTP_D28"))
qdat$Gene <- factor(qdat$Gene, levels = c("7-beta-hsdh", "baiI", "baiE", "bsh"))

p1 <- ggplot(qdat, aes( x = Gene, y = group))+
  geom_point(aes(size= size, color=sig2, shape = shape))+
  geom_point(aes(size= size, fill = Effect_size, alpha= sig, color=sig2,  shape = shape),stroke=0)+
  scale_fill_gradientn(
    colours = color,
    guide=guide_colourbar(ticks=T,nbin=50,barheight=.5, label=T,barwidth=10)
  )+scale_color_manual(values = c("black","grey"))+
  scale_shape_manual(values = c(24, 25, 1))+
  theme_minimal()+
  nature_theme2+
  theme(legend.position  = "none")+scale_size_continuous(range = c(1,5))

qdat_sub <- readRDS("../inputdata/acar_sp.rds")
qdat_sub$group <- factor(qdat_sub$group, levels = c("ZhangXY_2022_Acarbose_D168","GuYY_2017_Acarbose_D90", "ZhaoLP_2018_Acarbose_D84", "ZhaoLP_2018_Acarbose_D56", "ZhaoLP_2018_Acarbose_D28",  "ZhaoLP_2018_Acar+WTP_D84", "ZhaoLP_2018_Acar+WTP_D56", "ZhaoLP_2018_Acar+WTP_D28"))

p2 <- ggplot(qdat_sub, aes( x = tax, y = group))+
  geom_point(aes(size= size, color=sig2, shape = shape))+
  geom_point(aes(size= size, fill = effect_size, alpha= sig, color=sig2,  shape = shape),stroke=0)+
  scale_fill_gradientn(
    colours = color,
    guide=guide_colourbar(ticks=T,nbin=50,barheight=.5, label=T,barwidth=10)
  )+scale_color_manual(values = c("black","grey"))+
  scale_shape_manual(values = c(25, 24, 1))+
  theme_minimal()+
  nature_theme2+
  theme(legend.position  = "none", axis.text.y = element_text())+scale_size_continuous(range = c(1,5))

library(aplot)
out <- p2 %>% insert_left(p1, width = 0.3)

ggsave("../result/Figure3/fig3e.pdf",plot = out, device  = "pdf", width = 8, height = 4)

# fig3f




