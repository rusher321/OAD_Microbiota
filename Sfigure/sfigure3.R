#---------------------------- sfig3 ------------------------------------------#
library(ggpubr)
library(ggplot2)


## sfigure_a
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

ggsave(filename = "../result/Sfig/sfig3b_base.shannon.pdf", plot = p_shanno,
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

ggsave(filename = "../result/Sfig/sfig3a_base.rich.pdf", plot = p_rich,
       device = "pdf", width = 7, height = 4)



## sfigure_b
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
  ggsave(filename = paste0("../result/Sfig/sfig3b_base_adonis.",j,".pdf"),  netlist[[j]], device = "pdf", width = 10, height = 6 )
}

