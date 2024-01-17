qdat <- readRDS("../inputdata/sfigure1.RDS")

library(meta)
library(metafor)

metaplot <- function(qdat, var){
  
  meta_res <- metagen(qdat[, paste0(var, "_md")],
                      qdat[, paste0(var, "_std")],
                      data = qdat,
                      studlab = paste(label2),
                      comb.fixed = FALSE,
                      comb.random = TRUE,
                      method.tau = "DL",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "SMD") 
  forest.meta(meta_res,  xlab = var)
  
}

### res 
var <- c("HbA1c", "FPG", "PPG120", "Ins0", "Ins120", "HOMA_IR")
dir.create("../result/Figure1_sfigure1_3", recursive = T)

### plot 

for(i in var){
  pdf(file = paste0("../result/Figure1_sfigure1_3/", i, "meta.pdf"), width = 10, height = 6)
  metaplot(qdat, i)
  dev.off()
}
