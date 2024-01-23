#--------------------- Gee ---------------------------# 
source("../script/fun.R")
library(geeM)
pro <- readRDS("../inputdata/all_pro.rds")
phe <- readRDS("../inputdata/phe_gee.Rds")

phelist <- split(phe, phe$group)

## select species and all the KO 
com_clr <- readRDS("../inputdata/clr_compare.Rds")
load(file = "../inputdata/spe_intersect.Rds")

pro_sub <- pro[, c(downfeature, upfeature, colnames(pro)[136:153])]

reslist <- lapply(phelist,  function(x) {
  phe_tmp <- x[,-2]
  res_tmp <- GEEAna(dataset = pro_sub, metadata = phe_tmp, confounder = c("Age", "Sex", "BMI"),
                    timevar = "Time", scale = F, IDvar = "PID")
  res_tmp
})

colnames(reslist$WuH_2017_Metformin) <- colnames(reslist$GuYY_2017_Acarbose)
qdat <- do.call("rbind", reslist)
study_order <- c("GuYY_2017_Acarbose", "ZhangXY_2021_Acarbose", "ZhangYF_2020_Berberine", 
                 "RenHH_2023_Metformin", "WuH_2017_Metformin", "ZhangXY_2021_Vlidagliptin",
                 "GuYY_2017_Glipizide", "ZhangYF_2020_Placebo")

qdat_down <- qdat[study_order, 1:24]
qdat_up <- qdat[study_order, 25:90]
qdat_ko <- qdat[study_order, 91:126][, -c(19:24)]
qdat[study_order, 91:126]
pdf(file = "../revision/down_gee.pdf", width = 6, height = 10)
corPlot_tmp(corres = qdat_down, cutoff = 1, adjust = T, tr = T, cluster = T, row_clust = T, col_clust = F)
dev.off()

pdf(file = "../revision/up_gee.pdf", width = 6, height = 10)
corPlot_tmp(corres = qdat_up, cutoff = 1, adjust = T, tr = T, cluster = T, row_clust = T, col_clust = F)
dev.off()

pdf(file = "../revision/ko_gee.pdf", width = 6, height = 10)
corPlot_tmp(corres = qdat_ko, cutoff = 1, adjust = T, tr = T, cluster = T, row_clust = T, col_clust = F)
dev.off()



