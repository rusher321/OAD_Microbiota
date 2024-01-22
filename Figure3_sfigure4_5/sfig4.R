#-------------------------------- sfig4------------------------------------#
# sfig4a
qdat <- readRDS("../inputdata/BA_meta.Rds")
box_plot <- function(qdat, x){
  
  qdat$Time <- ifelse(qdat$Time == "D84", "D90", qdat$Time)
  qdat$Cohort <- factor(qdat$Cohort, levels = c("GuYY_2017_Acarbose", "ZhangYF_2020_Berberine", "GuYY_2017_Glipizide"))
  qdat[,x] <- qdat[,x]*100
  p1 <- ggplot(qdat, aes_string(x = "Time", y = x, color = "Cohort")) +
    geom_rect(xmin = 0.4, xmax = 2.5,
              ymin = -Inf, ymax = Inf,
              fill ='white',
              inherit.aes = F)+
    geom_jitter(alpha = 0.6)+
    geom_boxplot(alpha = 0.5) +
    scale_color_manual(values = c("#1B9E77","#D95F02", "#66A61E"))+
    #stat_compare_means(comparisons = my_comparisons,
    #                   method = "wilcox.test", paired = T)+
    theme_bw() + 
    xlab("") +
    ylab(paste0(x, " (%)"))+scale_y_continuous(expand = c(0.1,0.1))+
    theme(panel.grid=element_blank(),
          legend.position = "none",
          axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+facet_wrap(.~Cohort, scales = "free", nrow = 1)
  return(p1)
  
}

bavar <- c("DCA", "UDCA")

plist <- lapply(bavar, box_plot, qdat= qdat)

out <- cowplot::plot_grid(plotlist = plist, align = "h", nrow = 2)

ggsave(plot = out, filename = "../result/Figure3/sfig4a.pdf", device = "pdf",
       width = 4, height = 6)

# sfig4b 
spe_com <- readRDS("../inputdata/clr_compare.Rds")
spe_com <- spe_com[spe_com$drug %in% c("Acarbose", "Berberine"), ]
spe_com$shape <- ifelse(spe_com$effect_size > 0 , "+", ifelse(spe_com$effect_size < 0 , "-", "o"))
spe_com$size <- abs(spe_com$effect_size)
spe_com$sig <- ifelse(spe_com$p.adjust < 0.05, 1, 0.9)

colt <- c("#4C38CB", "#9191C8", "#DADAEC", 
          "#F0C1C1", "#E28383", "#D44545", "#CD2626")

color <- colorRampPalette(colt)(20)

nature_theme2 <- ggplot2::theme(
  axis.text.x = ggplot2::element_text(size = 8, vjust = 1, hjust = 1, angle = 45),
  axis.text.y = ggplot2::element_text(size = 8, hjust = 1),
  legend.title = ggplot2::element_text(size = 6, face = 'bold'),
  legend.text = ggplot2::element_text(size = 6),
  axis.line = ggplot2::element_line(colour = 'black', size = .25),
  axis.line.x = ggplot2::element_line(colour = 'black', size = .25),
  axis.line.y = ggplot2::element_line(colour = 'black', size = .25),
  panel.border = ggplot2::element_blank(),
  panel.grid.major = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank())

spe_com$sig2 <- ifelse(spe_com$sign_p.value<0.05,"P<0.05","P>=0.05")
spe_com$shape2 <- ifelse(spe_com$effect_size>0, 2, 6)
spe_fun_sort <- read.csv("../inputdata/fig3.species.sort.tmp.csv",row.names = 1)
spe_sort <- rownames(spe_fun_sort)
spe_sort2 <- spe_sort[c(41:49, 1:40,50)]

qdat <- spe_com[spe_com$tax %in% spe_sort2,]
qdat$tax <- as.factor(qdat$tax)
qdat$tax <- factor(qdat$tax, levels = rev(spe_sort2))
qdat$group <- factor(qdat$group, levels = study_name_v2)

p1 <- ggplot(qdat, aes( y = tax, x = group))+
  geom_point(aes(size= size, color=sig2, shape = shape))+
  geom_point(aes(size= size,fill = -effect_size, alpha= sig, color=sig2,  shape = shape),stroke=0)+
  scale_fill_gradientn(
    colours = color,
    guide=guide_colourbar(ticks=T,nbin=50,barheight=.5, label=T,barwidth=10)
  )+
  scale_color_manual(values = c("black","grey"))+
  scale_shape_manual(values = c(24, 25, 1))+
  theme_minimal()+
  nature_theme2+
  theme(legend.position = "top")+scale_size_continuous(range = c(1,3))

############# sandy plot 
library(ggalluvial)
library(reshape2)
library(hash)

spe_fun_sort$species <- rownames(spe_fun_sort)
qdat <- melt(spe_fun_sort)
qdat <- qdat[qdat$value !=0, ]
qdat$direction <- ifelse(qdat$value >0 , "pos", "neg")
tmp_stat <- 1/table(qdat$species)
h <- hash(names(tmp_stat), tmp_stat)
qdat$fre <- unlist(lapply(qdat$species, function(x){h[[x]]}))

qdat$variable <- factor(qdat$variable, levels = c("BaiE", "BaiI"))
qdat$species <- factor(qdat$species, rownames(spe_fun_sort)[c(41:49,1:40,50)])
tmp <- data.frame(species = rownames(spe_fun_sort), variable = "tmp", value = 1, direction = "tmp", fre = 0.01)
tmp_qdat <- rbind(tmp, qdat)
tmp_qdat$species <- factor(tmp_qdat$species, rownames(spe_fun_sort)[c(41:49,1:40,50)])
qdat$variable <- factor(qdat$variable, levels = c("BaiI", "BaiE"))
p2 <- ggplot(qdat,
             aes(y = fre,
                 axis2 = species, axis1 = variable)) +
  geom_flow(aes(fill = direction)) +
  geom_stratum() +
  geom_text(stat = "stratum", size = 2,aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("species", "variable"))+scale_fill_manual(values = c("#ceefe4", "#d9d9f3", "white"))+theme_classic()+
  theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank())

library(patchwork)
out <- p2+p1
ggsave(plot = out,filename =  "../result/Figure3/sfig4b.pdf", device = "pdf", width = 12, height = 8)


