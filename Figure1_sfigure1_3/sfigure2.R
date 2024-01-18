#---------------------------- sfig2 ------------------------------------------# 
library(ggpubr)
library(ggplot2)

## sfigure_a 
phe <- readRDS("../inputdata/phenotype_inf.Rds")
pro_tax_com <- readRDS("../inputdata/species_pro.Rds")

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
plist <- list()
for(i in 1:4){
  library(ape)
  library(ade4)
  pcoa_base <- dudi.pco(as.dist(distList[[i]][baseid, baseid]), scannf = F, nf = 2)
  var1 <- round((pcoa_base$eig[1]/sum(pcoa_base$eig)) * 100, 2)
  var2 <- round((pcoa_base$eig[2]/sum(pcoa_base$eig)) * 100, 2)
  pcoa_base_d <- as.data.frame(pcoa_base$li)
  pcoa_base_d$study <- paste0(basephe$Cohort, "_", basephe$Drug)
  factor_study <- c("GuYY_2017_Acarbose", "ZhangXY_2021_Acarbose", "ZhangYF_2020_Berberine","WuH_2017_Metformin",
                    "RenHH_2023_Metformin", "ZhangXY_2021_Vlidagliptin","GuYY_2017_Glipizide",  "ZhangYF_2020_Placebo")
  
  #color_study <-  c("#963F2D", "#FCCE8E", "#EFC94A", "4B778D", "#0CAFA9",  "#D2E69C","#F6A3BF", "#D891F2","#E1BC91", "#E3D0B9")
  color_study <- c("#1B9E77", "#A6761D", "#D95F02", "#7570B3",  "#1F78B4", "#E6AB02", "#66A61E","#E7298A")
  
  pcoa_base_d$study <- factor(pcoa_base_d$study, levels = factor_study)
  
  plist[[i]] <- ggplot(pcoa_base_d, aes(x=A1, y=A2,  color=study)) +
    geom_point() +
    xlab(paste("PCoA1: ", var1, "%")) +
    ylab(paste("PCoA2: ", var2, "%")) +
    theme_ipsum(base_family = "Helvetica")+
    ggtitle(paste0("baseline based ", distVar[i]))+scale_color_manual(values = color_study)+
    stat_ellipse(level = 0.8)
  
}

out <- cowplot::plot_grid(plotlist = plist, nrow = 2)

ggsave("../result/Figure1_sfigure1_3/sfig2a_all_base.pcoa.pdf", plot = out, device = "pdf", width = 15, height = 10)


## sfigure_b_c
base_diver <- readRDS("../inputdata/base_diver.Rds")
color_value <-  c("#1B9E77", "#66A61E", "#7570B3", "#1F78B4",  "#A6761D", 
                  "#E6AB02", "#D95F02", "#E7298A")

p_shanno <- ggplot(base_diver, aes(x = group, y = Shannon, fill = group))+
  geom_boxplot()+stat_compare_means(label = "p.signif", hide.ns = TRUE, ref.group = ".all.")+
  scale_fill_manual(values = color_value)+
  theme_classic()+
  theme(axis.title = element_text(size = 13,color="black"),
        axis.text = element_text(size = 12,color = "black"),
        #axis.text.y = element_text(size = 11,color = rev(dat3$Color)),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.title = element_text(size = 13,color = "black"),
        legend.text = element_text(size = 12,color = "black"))+guides(fill="none")

ggsave(filename = "../result/Figure1_sfigure1_3/sfig2b_base.shannon.pdf", plot = p_shanno, 
       device = "pdf", width = 7, height = 4)

p_rich <- ggplot(base_diver, aes(x = group, y = Richness, fill = group))+
  geom_boxplot()+stat_compare_means(label = "p.signif", hide.ns = TRUE, ref.group = ".all.")+
  scale_fill_manual(values = color_value)+
  theme_classic()+
  theme(axis.title = element_text(size = 13,color="black"),
        axis.text = element_text(size = 12,color = "black"),
        #axis.text.y = element_text(size = 11,color = rev(dat3$Color)),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.title = element_text(size = 13,color = "black"),
        legend.text = element_text(size = 12,color = "black"))+guides(fill="none")

ggsave(filename = "../result/Figure1_sfigure1_3/sfig2c_base.rich.pdf", plot = p_rich, 
       device = "pdf", width = 7, height = 4)



## sfigure_d 
load("../inputdata/sfigure2c.Rdata")

netlist <- list()
var <- c("helli", "jsd", "spear" , "bray")

library(igraph)
library(ggraph)

for(i in 1:length(var)){
  
  output_plot <- output_list[[i]][-c(3,5), -c(3,5)]
  colnames(output_plot)[4] <- rownames(output_plot)[4] <- "RenHH_2023_Metformin"
  qdat1 <- melt(output_plot)
  colnames(qdat1) <- c("from", "to", "r2")
  
  output_plot_p <- outputp_list[[i]][-c(3,5), -c(3,5)]
  colnames(output_plot_p)[4] <- rownames(output_plot_p)[4] <- "RenHH_2023_Metformin"
  qdat1_p <- melt(output_plot_p)
  qdat1$p <- ifelse(qdat1_p[,3] < 0.05, "pink", "grey")
  
  nodes <- as.character(unique(qdat1$from))
  node_Data <- data.frame(node = nodes, drug = sapply(strsplit(nodes, split = "_"), function(x){x[3]}))
  igraph_data <- graph_from_data_frame(d = qdat1, directed = F, vertices = node_Data)
  lay = create_layout(igraph_data, layout = "circle")
  study_name_ord <- c("GuYY_2017_Acarbose_D90", "ZhangXY_2021_Acarbose_D168", 
                      "ZhangYF_2020_Berberine_D84", "WuH_2017_Metformin_D60",
                      "RenHH_2023_Metformin_D90", 
                      "ZhangXY_2021_Vlidagliptin_D168", 
                      "GuYY_2017_Glipizide_D90",  "ZhangYF_2020_Placebo_D84")
  
  color_name_ord <- c("#1B9E77", "#A6761D", "#D95F02", "#7570B3",  "#1F78B4", "#E6AB02", "#66A61E","#E7298A")
  
  color_value <- color_name_ord[pmatch(lay$name, study_name_ord)]
  
  netlist[[var[i]]] <- ggraph(lay) + 
    geom_edge_link(aes(alpha = r2, width = r2, color = p)) + 
    scale_edge_color_manual(values = c("#ffb6b9", "#bbded6"))+
    geom_node_point(aes(color = name), size = 15) +
    theme_graph(base_family =  "Helvetica")+geom_node_text(aes(label = name), repel=TRUE, nodge_x = 0.5)+scale_color_manual(values = color_value)
}

for(j in var){
  ggsave(filename = paste0("../result/Figure1_sfigure1_3/sfig2_base_adonis.",j,".pdf"),  netlist[[j]], device = "pdf", width = 10, height = 6 )
}

