#------------------------------ sfig6-----------------------------------------#
feature_ass <- readRDS("../../../OAD_Microbiota/inputdata/sfig6.rds")
feature_ass2 <- feature_ass[feature_ass$var %in% c("HbA1c", "FPG", "HOMA_IR"), ]
feature_post <- feature_ass2[feature_ass2$estimate > 0 & feature_ass2$group != "ZhangYF_2020_Placebo", ]
feature_post$group <- gsub("YangFM_2021_Metformin", "RenHH_2023_Metformin", feature_post$group)
plot_list <- lapply(split(feature_post, f = feature_post$group), function(x){unique(x$variable)})

# detail
com_pos <- combn(names(plot_list), 2)
interall <- function(x){
  for(i in 1:ncol(x)){
    study1 <- x[1, i]
    study2 <- x[2, i]
    inter_bac <- intersect(plot_list[[study1]], plot_list[[study2]])
    if(length(inter_bac)!=0){
      print(paste0(study1, " ", study2, ": ",paste0(inter_bac, "; ")))
    }
  }
}

interall(com_pos)

library(venn)
library(ggplot2)
library(ggpolypath)

p1 <- venn(plot_list, ggplot= T, zcolor = brewer.pal(n = 8, name = "Dark2"), ilcs = 2)+theme_classic()
ggsave("../result/Sfig/sfig9_pcor.positive.venn.pdf", p1, device = "pdf", width = 6, height = 6)

feature_neg <- feature_ass2[feature_ass2$estimate < 0 & feature_ass2$group != "ZhangYF_2020_Placebo", ]
plot_list <- lapply(split(feature_neg, f = feature_neg$group), function(x){unique(x$variable)})
com_neg <- combn(names(plot_list),2)
interall(com_neg)

p2 <- venn(plot_list, ggplot= T, zcolor = brewer.pal(n = 8, name = "Dark2"), ilcs = 2)+theme_classic()
ggsave("../result/Sfig/sfig9_pcor.negative.venn.pdf", p2, device = "pdf", width = 6, height = 6)
