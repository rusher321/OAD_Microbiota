#------------------------------- Fig2 --------------------------------------#
# figa
qdat <- readRDS("../inputdata/fig4_num.rds")
qdat$Var2 <- factor(qdat$Var2, levels = rev(c("GuYY_2017_Acarbose", "ZhangXY_2021_Acarbose", "ZhangYF_2020_Berberine",
                                              "YangFM_2021_Metformin", "ZhangXY_2021_Vlidagliptin" ,"GuYY_2017_Glipizide", "ZhangYF_2020_Placebo")))
qdat$Var1 <- factor(qdat$Var1, levels = c("HbA1c", "FPG", "PPG120", "Ins0", "Ins120", "HOMA_IR"))

qdat_fig <- ggplot(qdat, aes(x= Var1, y = Var2, fill= Freq)) +
  geom_tile() +
  scale_fill_distiller(palette = "Purples", direction = 1) +
  theme_ipsum(base_family = "Helvetica")+
  geom_text(label = round(qdat$Freq,0), size = round(qdat$Freq/3,0))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "../result/Figure2/Fig2a_left.pdf", qdat_fig,
       device = "pdf", width = 6, height = 4)

out_plot <- readRDS("../inputdata/Fig4_r2.rds")

out_fig <- ggplot(out_plot, aes(x= V2, y = V1, fill= V4)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdPu", direction = 1) +
  theme_ipsum(base_family = "Helvetica")+
  geom_text(label = round(out_plot$V4,2), size = round(out_plot$V4*8,0))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(filename = "../result/Figure2/Fig2a_right.pdf", out_fig,
       device = "pdf", width = 6, height = 4)

# figb

qdat <- readRDS("../inputdata/Fig4_bacteroides.rds")

p <-
  ggplot(qdat, aes(x=bacteria, y= OR, group = drug, fill = drug))+
  geom_errorbar(aes(ymin=lower, ymax=upper), linetype = "dashed", width = .2, position=position_dodge(width=c(0.6)))+
  geom_point(size=4, shape = 23, position=position_dodge(width=c(0.6)))+
  geom_hline(yintercept = 1,linetype="dashed")+
  scale_fill_manual(values = c("lightblue","#F4ABC4"))+
  xlab("")+theme_classic()+
  theme(
    axis.text = element_text(size = 12,color=1),
    axis.text.y = element_text(face = "italic"),
    axis.title = element_text(size = 14,color=1),
    strip.text = element_text(size = 13,color=1),
    strip.background = element_rect(fill="#4393C3")
  )+scale_y_continuous(breaks = c(0.5,1,1.6), labels = c(0.6, 1, 1.6))+coord_flip()

ggsave(filename = "../result/Figure2/Fig2b_Bacteria.or.pdf", plot = p, device = "pdf", width = 6, height = 4)


# figc/d
changePlist <- readRDS("../inputdata/Fig4_changeP.rds")
pro_com <- readRDS("../inputdata/all_pro.rds")
dat_list <- lapply(changePlist, function(x){y = cbind(x, pro_com[rownames(x),]);y})
names(dat_list)

library(cowplot)
scatterplotWithSCC <- function (dat, x, y, group = NULL,PCC=F,adj=NULL,title)
{
  dat <- dat[!is.na(dat[, x]), , drop = F]
  s0 <- cor.test(dat[,x],dat[,y],method = "s")
  if(x=="7-beta-hsdh"){
    colnames(dat)[which(colnames(dat)==x)] <- "beta_hsdh"
    x <- "beta_hsdh"
  }
  if(PCC){
    if(is.null(adj)){
      stop("If perform PCC, adj must be support")
    }
    dat <- dat[,c(x,y,adj,group),drop=F]
    dat <- na.omit(dat)
    adj <- dat[,adj,drop=F]
    for(i in 1:ncol(adj)){
      adj[,i] <- as.numeric(adj[,i])
    }

    s0 <- ppcor::pcor.test(dat[, x,drop=F], dat[, y,drop=F],adj, method = "s")
  }

  if (is.null(group)) {
    lab <- paste0("scc rho = ", round(s0$estimate, 3), "; p = ",
                  formatC(s0$p.value, digits = 2))
    p <- ggplot(dat, aes_string(x, y)) +
      geom_point(size = 4, alpha = 0.55, colour = "White",shape = 21,fill = "Red") +
      geom_smooth(method = "lm", colour = "white", alpha = 0.2) +
      annotate("text", x = -Inf, y = Inf, vjust = 1.2,
               hjust = 0, label = lab, size = 3) +
      theme_bw()+ggtitle(label = title)+scale_x_log10()
  }
  else {
    lst_levels <- levels(as.factor(dat[, group]))
    lab <- c()
    for(l in lst_levels){
      id <- dat[,group]==l
      a <- cor.test(dat[id,x],dat[id,y],method = "s")
      if(PCC){
        a <- ppcor::pcor.test(dat[id,x,drop=F],dat[id,y,drop=F],adj[id,,drop=F],method = "s")
      }
      lab <- c(lab,paste0(l," scc rho=",round(a$estimate,3),"; p=",formatC(a$p.value, digits = 2)))
    }
    lab <- paste(lab,collapse = "\n")
    lab <- paste0(paste0("Totol rho = ", round(s0$estimate, 3), "; p = ",formatC(s0$p.value, digits = 2)),"\n",lab)
    p <- ggplot(dat, aes_string(x, y, color = group)) +
      geom_point(size = 4, alpha = 0.5) +
      #geom_smooth(method = "loess",se = F) +
      geom_smooth(data = dat, method = "lm",se = F, aes_string(x, y), color = "black") +
      annotate("text", x = -Inf, y = Inf, vjust = 1.2, hjust = 0, label = lab,size = 3) +
      theme_bw()+ coord_trans(x = "log2")
  }
  p
}

# bsh
p1 <- scatterplotWithSCC(dat = dat_list$GuYY_2017_Acarbose, x = "bsh", y = "HOMA_IR",title= "GuYY_2017_Acarbose")
p2 <- scatterplotWithSCC(dat = dat_list$YangFM_2021_Metformin, x = "bsh", y = "HbA1c",title= "RenHH_2023_Metformin")

# K00005 K00864 K05878 K05879 K05881
kolist <- c("K00005",  "K05878", "K05879")
acar_ko <- lapply(kolist, scatterplotWithSCC, dat= dat_list$YangFM_2021_Metformin, y = "HbA1c", title = "RenHH_2023_Metformin")
acar_ko_p <- plot_grid(plotlist = list(p2, acar_ko[[1]], acar_ko[[2]], acar_ko[[3]]), align = "h", nrow = 1)

vild_ko <- lapply(kolist, scatterplotWithSCC, dat= dat_list$ZhangXY_2021_Vlidagliptin, y = "HbA1c", title = "ZhangXY_2022_Vild")
vild_ko_p <- plot_grid(plotlist = list(p1, vild_ko[[1]], vild_ko[[2]], vild_ko[[3]]), align = "h", nrow = 1)

out <- acar_ko_p/vild_ko_p
ggsave(filename = "../result/Figure2/Fig2c_d.pdf", plot = out,
       device = "pdf", width = 10, height = 6)







