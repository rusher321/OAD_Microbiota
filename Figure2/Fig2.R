#---------------------------- Fig2 network ----------------------------------# 
phe <- readRDS("../inputdata/phenotype_inf.Rds")
phe$group <-  paste0(phe$Cohort, "_", phe$Drug)
pro <- readRDS("../inputdata/species_pro.Rds")
out_rank <- readRDS("../inputdata/Fig1_tax_rank.Rds")

library(NetCoMi)
outlist <- readRDS(file = "../inputdata/diff.stat.RDS")

# Fig2b 

qdat <- do.call("cbind", list(outlist$GuYY_2017_Acarbose[[2]], 
                              outlist$ZhangXY_2021_Acarbose[[2]], 
                              outlist$ZhangYF_2020_Berberine[[2]]))
qdat <- log2(qdat+1)
qdat <- qdat[rowSums(qdat) >2, ]
qdat[,c(2,4,6)] <- -qdat[,c(2,4,6)]
qdat <- qdat[-29, ] # remove the phlyum Candidatus_Saccharibacteria

species_family <- out_rank[rownames(qdat), 2, drop = F] 
species_family[,1] <- gsub("p__", "", species_family[,1])
bac_list <- rownames(species_family[species_family$family == "Firmicutes", ,drop=F])

library(RColorBrewer)

colt <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(20)
colt2 <- c(colt[1:10], "#f5f5f5", "#f5f5f5", colt[11:20])
tmp_colt <- colorRampPalette(c("#187498", "#f5f5f5", "#EB5353"))(30)
tmp_colt2 <- c(tmp_colt[5:14],"#f5f5f5", "#f5f5f5", "#f5f5f5", "#f5f5f5", tmp_colt[18:30])

colnames(species_family) <- "phylum"
ann_col <- list(phylum = c(Actinobacteria = "#E7298A80", Bacteroidetes = "#D95F0280", Firmicutes = "#1B9E7780",Proteobacteria = "#7570B380"))

qdat <- as.data.frame(qdat)
qdat$phylum <- species_family$phylum
qdat$phylum <- factor(qdat$phylum, levels = c("Firmicutes", "Bacteroidetes", "Actinobacteria", "Proteobacteria"))
qdat <- qdat[order(qdat[,5], decreasing = T), ]
qdat <- qdat[order(qdat$phylum), ]

pdf("../result/Figure2/Fig2b.diff_degree.pdf", width = 20, height = 5)
pheatmap(t(qdat[,c(1,3,5,2,4,6)]), display_numbers = F, cluster_rows = F,  color = tmp_colt2, cluster_cols = F, annotation_col  = species_family, angle_col = "45", annotation_colors = ann_col, gaps_row = c(3), gaps_col = c(63, 91, 99))
dev.off()

# Fig2c  generate the input of cytoscape 

namelist <- names(outlist)
for(i in namelist){
  
  netPlot <- outlist[[i]][[1]]
  netPlot[lower.tri(netPlot)] <- 0
  netPlot2 <- melt(netPlot)
  netPlot2 <- netPlot2[netPlot2$value !=0 , ]
  colnames(netPlot2) <- c("from", "to", "size")
  netPlot2$from <- as.character(netPlot2$from)
  netPlot2$to <- as.character(netPlot2$to)
  
  vertices <- outlist[[i]][[2]]
  netlist <- list(vertices = data.frame(label = rownames(vertices), value = apply(vertices, 1, sum)), edges= netPlot2)
  
  rownames(netlist$vertices) <- netlist$vertices$label
  ann_col <- out_rank[rownames(vertices), 2, drop= F]
  colnames(ann_col) <- "Phlyum"
  netlist$vertices$phylum <- gsub("p__", "",out_rank[rownames(netlist$vertices), 2])
  
  out_vertice <- paste0("../result/Figure2/", i, ".vertices.csv")
  write.csv(netlist$vertices, out_vertice, quote = F)
  
  netlist$edges$size2 <- ifelse(netlist$edges$size >0, "pos", "neg")
  netlist$edges$size3 <- abs(netlist$edges$size)
  out_edge <- paste0("../result/Figure2/", i, ".edges.csv")
  
  write.csv(netlist$edges, out_edge, quote = F)
  
}

