#----------------------------- sfig5 ------------------------------------#
# sfig5a 
qvalue <- readRDS("../inputdata/sfig5a.rds")
lt <- apply(qvalue, 2, function(x){rownames(qvalue_plot)[x<0.05]})

library(ComplexHeatmap)
m <- m_spe <- make_comb_mat(lt)
ss = set_size(m)
cs = comb_size(m)

order <- c("ZhaoLP2_2018_Acarbose_D28", "ZhaoLP2_2018_Acarbose_D56", "ZhaoLP2_2018_Acarbose_D84", "GuYY_2017_Acarbose_D90", "ZhangXY_2021_Acarbose_D168", 
           "ZhaoLP2_2018_Acar+WTP_D28", "ZhaoLP2_2018_Acar+WTP_D56", "ZhaoLP2_2018_Acar+WTP_D84")

pdf("../result/Figure3/sfig5a.acarbose.compare.upset.pdf", width = 15, height = 4)
ht = UpSet(m, 
           #set_order = order(ss),
           set_order = pmatch(order, names(ss)),
           comb_order = order(-cs),
           top_annotation = HeatmapAnnotation(
             "Genre Intersections" = anno_barplot(cs, 
                                                  ylim = c(0, max(cs)*1.1),
                                                  border = FALSE, 
                                                  gp = gpar(fill = "black"), 
                                                  height = unit(4, "cm")
             ), 
             annotation_name_side = "left", 
             annotation_name_rot = 90),
           left_annotation = rowAnnotation(
             "Significant Species" = anno_barplot(-ss, 
                                                  baseline = 0,
                                                  axis_param = list(
                                                    at = c(0, -20, -40, -60, -80),
                                                    labels = c(0, 20, 40, 60, 80),
                                                    labels_rot = 0),
                                                  border = FALSE, 
                                                  gp = gpar(fill = "black"), 
                                                  width = unit(4, "cm")
             ),
             set_name = anno_text(set_name(m), 
                                  location = 0.5, 
                                  just = "center",
                                  width = max_text_width(set_name(m)) + unit(4, "mm"))
           ), 
           right_annotation = NULL,
           show_row_names = FALSE)
ht = draw(ht)
od = column_order(ht)
decorate_annotation("Genre Intersections", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
})
dev.off()

# sfig5b
library(RColorBrewer)
library(ggplot2)
load(file = "../inputdata/sfig5b.Rdata")
colt <- colorRampPalette(rev(brewer.pal(n = 7, name ="PiYG")))(20)

p <- ComplexHeatmap::pheatmap(as.matrix(qdat2)[, c(6:8,1:5)], display_numbers = qdatNum[,c(6:8,1:5)], cluster_rows = F,  color =colt, cluster_cols = F, fontsize_number = 20, number_color = "White", angle_col = "90", annotation_row  = con_acar[all_dietspe, ,drop=F], annotation_colors = list(V2 = c(negative = "#1B9E77", positive = "#D95F02")), gaps_col = c(5))

pdf("../result/Figure3/sfig5b.acar.spe.pdf", width = 5, height = 8)
print(p)
dev.off()


