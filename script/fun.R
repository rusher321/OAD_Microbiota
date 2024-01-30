maaslin2_heatmap <-
  function(
    output_results,
    title = NA,
    cell_value = 'qval',
    data_label = 'data',
    metadata_label = 'metadata',
    border_color = 'grey93',
    color = colorRampPalette(c("darkblue", "grey90", "darkred")),
    col_rotate = 90,
    first_n = 50) {
    
    # read MaAsLin output
    df <- read.table(
      output_results,
      header = TRUE,
      sep = "\t",
      fill = TRUE,
      comment.char = "" ,
      check.names = FALSE
    )
    
    title_additional <- ""
    
    title_additional <- ""
    if (!is.na(first_n) & first_n > 0 & first_n < dim(df)[1]) {
      if (cell_value == 'coef') {
        df <- df[order(-abs(df[[cell_value]])) , ]
      } else{
        df <- df[order(df[[cell_value]]), ]
      }
      # get the top n features with significant associations
      df_sub <- df[1:first_n,]
      for (first_n_index in seq(first_n, dim(df)[1]))
      {
        if (length(unique(df_sub$feature)) == first_n)
        {
          break
        }
        df_sub <- df[1:first_n_index,]
      }
      # get all rows that have the top N features
      df <- df[which(df$feature %in% df_sub$feature),]
      title_additional <- paste("Top", first_n, sep=" ")
    }
    
    if (dim(df)[1] < 2) {
      print('There are no associations to plot!')
      return(NULL)
    }
    
    metadata <- df$metadata
    data <- df$feature
    dfvalue <- df$value
    value <- NA
    
    # values to use for coloring the heatmap
    # and set the colorbar boundaries
    if (cell_value == "pval") {
      value <- -log(df$pval) * sign(df$coef)
      value <- pmax(-20, pmin(20, value))
      if (is.null(title))
        title <- "(-log(pval)*sign(coeff))"
    } else if (cell_value == "qval") {
      value <- -log(df$qval) * sign(df$coef)
      value <- pmax(-20, pmin(20, value))
      if (is.null(title))
        title <- "(-log(qval)*sign(coeff))"
    } else if (cell_value == "coef") {
      value <- df$coef
      if (is.null(title))
        title <- "(coeff)"
    }
    
    if (title_additional!="") {
      title <- paste(title_additional, "features with significant associations", title, sep=" ")
    } else {
      title <- paste("Significant associations", title, sep=" ")
    }
    
    # identify variables with more than one level present
    verbose_metadata <- c()
    metadata_multi_level <- c()
    for (i in unique(metadata)) {
      levels <- unique(df$value[df$metadata == i])
      if (length(levels) > 1) {
        metadata_multi_level <- c(metadata_multi_level, i)
        for (j in levels) {
          verbose_metadata <- c(verbose_metadata, paste(i, j))
        }
      } else {
        verbose_metadata <- c(verbose_metadata, i)
      }
    }
    
    n <- length(unique(data))
    m <- length(unique(verbose_metadata))
    
    if (n < 2) {
      print(
        paste(
          "There is not enough features in the associations",
          "to create a heatmap plot.",
          "Please review the associations in text output file.")
      )
      return(NULL)
    }
    
    if (m < 2) {
      print(
        paste(
          "There is not enough metadata in the associations",
          "to create a heatmap plot.",
          "Please review the associations in text output file.")
      )
      return(NULL)
    }
    
    a = matrix(0, nrow = n, ncol = m)
    a <- as.data.frame(a)
    
    rownames(a) <- unique(data)
    colnames(a) <- unique(verbose_metadata)
    
    for (i in seq_len(dim(df)[1])) {
      current_metadata <- metadata[i]
      if (current_metadata %in% metadata_multi_level) {
        current_metadata <- paste(metadata[i], dfvalue[i])
      }
      if (abs(a[as.character(data[i]), 
                as.character(current_metadata)]) > abs(value[i]))
        next
      a[as.character(data[i]), as.character(current_metadata)] <- value[i]
    }
    
    # get the range for the colorbar
    return(a)
  }



scatterplotWithSCC <- function (dat, x, y, group = NULL,PCC=F,adj=NULL, color = NULL)
{
  dat <- dat[!is.na(dat[, x]), , drop = F]
  s0 <- cor.test(dat[,x],dat[,y],method = "s")
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
      geom_point(size = 0.8, alpha = 0.5) +
      geom_smooth(method = "rlm", se = F) +
      annotate("text", x = -Inf, y = Inf, vjust = 1.2,
               hjust = 0, label = lab, size = 3) + theme_bw()
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
      geom_point(size = 0.8, alpha = 0.5) +
      geom_smooth(method = "lm",se = F) +
      geom_smooth(data = dat, method = "rlm",se = F, aes_string(x, y), color = "black") +
      annotate("text", x = -Inf, y = Inf, vjust = 1.2, hjust = 0, label = lab,size = 3) +
      theme_bw()
    if(!is.null(color)){
      p <- p+scale_color_manual(values = color)
    }
    
  }
  p
}

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


study_name_v2 <- c("GuYY_2017_Acarbose_D90", "ZhangXY_2021_Acarbose_D168", 
                   "ZhangYF_2020_Berberine_D84", "WuH_2017_Metformin_D60",
                   "WuH_2017_Metformin_D120", "YangFM_2021_Metformin_D90", 
                   "ZhangXY_2021_Vlidagliptin_D168", 
                   "GuYY_2017_Glipizide_D90",  "ZhangYF_2020_Placebo_D84")

study_name_order <- c("GuYY_2017_Acarbose_D90", "ZhangXY_2021_Acarbose_D168", 
                   "ZhangYF_2020_Berberine_D84", "WuH_2017_Metformin_D60",
                   "WuH_2017_Metformin_D120", "RenHH_2023_Metformin_D90", 
                   "ZhangXY_2021_Vlidagliptin_D168", 
                   "GuYY_2017_Glipizide_D90",  "ZhangYF_2020_Placebo_D84")

GEEAna <- function (dataset, metadata, confounder, timevar, scale, IDvar,
                    ...)
{
  id <- intersect(rownames(dataset), rownames(metadata))
  dataset <- dataset[id, ]
  metadata <- metadata[id, ]
  print("confirm the sample ID is order by time")
  metadata <- metadata[order(metadata[, timevar]), ]
  dataset <- dataset[rownames(metadata), ]
  confounderindex <- which(colnames(metadata) %in% c(confounder,
                                                     timevar, IDvar))
  
  datacon <- metadata[, confounderindex]
  metadatafilter <- metadata[, -confounderindex, drop=F]
  result <- matrix(NA, nrow = ncol(metadatafilter), ncol = ncol(dataset) *
                     2)
  result <- as.data.frame(result)
  rownames(result) <- colnames(metadatafilter)
  for (i in 1:c(ncol(metadatafilter))) {
    for (j in 1:ncol(dataset)) {
      dat_com <- data.frame(x = metadatafilter[, i], y = dataset[,
                                                                 j], datacon, PatientID = metadata[, IDvar])
      formula <- formula(paste0("y~x+", paste(confounder,
                                              collapse = "+")))
      dat_com <- dat_com[!apply(dat_com, 1, function(x) {
        any(is.na(x))
      }), ]
      if (scale) {
        dat_com$x <- scale(invt(dat_com$x))
        dat_com$y <- scale(invt(dat_com$y))
      }
       tryCatch({geeInd <- geem(formula, id = PatientID, data = dat_com,
                     family = gaussian, corstr = "unstructured")
       tmp <- summary(geeInd)
       result[i, c((2 * j - 1):(2 * j))] <- c(tmp$wald.test[2],
                                              tmp$p[2])
       colnames(result)[c((2 * j - 1):(2 * j))] <- paste0(colnames(dataset)[j],
                                                          c("wald", "_p.value"))
       },# notice
                     error = function(e){
                       result[i, c((2 * j - 1):(2 * j))] <- c(0, 1)
                       colnames(result)[c((2 * j - 1):(2 * j))] <- paste0(colnames(dataset)[j],
                                                                          c("wald", "_p.value"))
                     })
      
    }
  }
  return(result)
}

corPlot_tmp <- function (corres, cutoff, adjust, tr,  cluster=F, row_clust = F, col_clust=F)
{
  trans <- function(x) {
    if (x <= 0.05 & x > 0.01) {
      out <- "*"
    }
    else if (x <= 0.01 & x > 0.001) {
      out <- "**"
    }
    else if (x <= 0.001) {
      out <- "***"
    }
    else {
      out <- " "
    }
    return(out)
  }
  trans2 <- function(x){
    if(x<=0.05){
      out <- "#"
    }else{
      out <- " "
    }
    return(out)
  }
  trans3 <- function(x, y){
    if(x != " "){
      out <- x
    }else{
      out <- y
    }
    return(out)
  }
  
  
  xname <- rownames(corres)
  corres <- apply(corres, 2, as.numeric)
  sp.corr.t <- corres
  rownames(sp.corr.t) <- xname
  index <- 2 * c(1:(ncol(sp.corr.t)/2))
  dat.pvalue <- sp.corr.t[, index]
  dat.pvalue[is.na(dat.pvalue)] <- 1
  dat.pvalue.tmp <- dat.pvalue
  dat.cor <- sp.corr.t[, -index]
  dat.cor[is.na(dat.cor)] <- 0
  colnames(dat.cor) <- gsub("_p.value", "", colnames(dat.pvalue))
  if (adjust) {
    dat.raw <- dat.pvalue
    dat.pvalue <- apply(dat.pvalue, 2, p.adjust, method = "BH")
  }
  rmindex <-     pvalue.index2 <- apply(dat.pvalue, 1, function(x) any(x < 
                                                                         cutoff))
  rmindex2 <-     pvalue.index2 <- apply(dat.pvalue, 2, function(x) any(x < 
                                                                          cutoff)) 
  dat.cor.cle <- dat.cor[, rmindex2]
  dat.pva.cle <- dat.pvalue[, rmindex2]
  dat.raw.cle <- dat.raw[, rmindex2]
  #dat.pvalue.tmp <- dat.pvalue.tmp[rmindex, paste0(selectF, "_p.value")]
  
  num <- matrix(NA, nrow = nrow(dat.pva.cle), ncol = ncol(dat.pva.cle))
  for (i in 1:ncol(dat.pva.cle)) {
    num[, i] <- mapply(trans, dat.pva.cle[, i])
  }
  
  num2 <- matrix(NA, nrow = nrow(dat.pva.cle), ncol = ncol(dat.pva.cle))
  for (i in 1:ncol(dat.pva.cle)) {
    num2[, i] <- mapply(trans2, dat.raw.cle[, i])
  }
  
  num3 <- matrix(NA, nrow = nrow(dat.pva.cle), ncol = ncol(dat.pva.cle))
  for (i in 1:ncol(dat.pva.cle)) {
    num3[, i] <- mapply(trans3, num[, i], num2[,i])
  }
  num <- num3
  
  colt <- c("#4C38CB", "#9191C8", "#DADAEC", 
            "#F0C1C1", "#E28383", "#D44545", "#CD2626")
  colt <- colorRampPalette(colt)(9)
  #colt <- c("#87CEEB", "#FFFFFF", "#FF69B4")
  if (tr) {
    dat.cor.cle <- t(dat.cor.cle)
    num <- t(num3)
  }
  
  wald <- max(max(dat.cor.cle), abs(min(dat.cor.cle)))
  wald <- 8
  gapwald1 <- seq(0, wald,  wald/4)
  gapwald2 <- seq(-wald, 0, wald/4)
  
  if(!cluster){
      ComplexHeatmap::pheatmap(dat.cor.cle, treeheight_row = 43, treeheight_col = 23,
                       cellwidth = 20, cellheight = 8, cluster_cols = col_clust, cluster_rows = row_clust,
                       fontsize_row = 8, fontsize_col = 13, show_colnames = T,
                       display_numbers = num, color =  colt,breaks = c(gapwald2,
                                                                       gapwald1[-1]) ,number_color = "black")
  }else{
    ComplexHeatmap::pheatmap(dat.cor.cle, treeheight_row = 43, treeheight_col = 23,
                       cellwidth = 20, cellheight = 8, cluster_cols = col_clust, cluster_rows = row_clust,
                       fontsize_row = 8, fontsize_col = 13, show_colnames = T,
                       display_numbers = num, color =  colt,breaks = c(gapwald2,
                                                                       gapwald1[-1]) ,number_color = "black")
  }
}

Uniqueness <- function(d){
  d <- as.matrix(d)
  u <- matrix(NA,nrow(d),1)
  u <- as.data.frame(u)
  rownames(u) <- rownames(d)
  colnames(u) <- "Uniqueness"
  for(x in 1:nrow(d)){
    a <- as.numeric(d[x,-x])
    u[x,1] <- min(a)
  }
  u
}

diversity_f <- function(dat){
  
  richness <- apply(dat, 1, function(x){sum(x!=0)})
  shannon <- diversity(dat, index = "shannon")
  invsimpson <- diversity(dat, index = "invsimpson")
  # transform clr
  
  distD1 <- vegan::vegdist(vegan::decostand(dat, method = "hellinger"),
                           method = "euclidean", binary = 1)
  distD2 <- philentropy::JSD(as.matrix(dat))
  rownames(distD2) <- colnames(distD2) <- rownames(dat)
  distD2 <- as.dist(distD2)
  
  distD3 <- as.dist(1 - cor(t(dat), method = "s"))
  
  d.kendall <- 1-cor.fk(t(dat))
  d.bray <- vegdist(dat, method = "bray")
  #dat_t <- zCompositions::cmultRepl(dat)
  #dat_clr <- robCompositions::cenLR(dat_t)$x.clr
  #d.aitchison <- vegdist(dat_clr, method="euclidean")
  
  uniq_bray <- Uniqueness(d.bray)
  uniq_kendall <- Uniqueness(d.kendall)
  #uniq_aitchison <- Uniqueness(d.aitchison)
  uniq_hell <- Uniqueness(distD1)
  uniq_jsd <- Uniqueness(distD2)
  uniq_spe <- Uniqueness(distD3)
  
  out <- data.frame(richness, shannon, invsimpson, uniq_bray, uniq_kendall,
                    uniq_hell, uniq_jsd, uniq_spe)
  rownames(out) <- rownames(dat)
  colnames(out)[4:8] <- c("Uniqueness_bray", "Uniqueness_kendall",
                          "Uniqueness_hellinger", "Uniqueness_JSD", "Uniqueness_spearman")
  return(out)
  
}

std_robust <- function(x){sd = jointseg::estimateSd(x, method = "von Neumann"); return(sd/sqrt(length(x[!is.na(x)])))}

sd_robust <- function(x){sd = jointseg::estimateSd(x, method = "von Neumann");sd}

metainf2 <- function(phe, feature = c("richness", "shannon", "Uniqueness_bray", "Uniqueness_kendall",
                                      "Uniqueness_hellinger", "Uniqueness_JSD", "Uniqueness_spearman")
){
  
  durgL <- unique(phe$Drug)
  timeL <- sort(unique(phe$Time))
  timeN <- length(timeL)
  #feature <-  c("richness", "shannon", "Uniqueness_bray", "Uniqueness_kendall",
  #              "Uniqueness_hellinger", "Uniqueness_JSD", "Uniqueness_spearman")
  
  
  out <- matrix(" ", nrow = length(durgL)*(length(timeL)-1), ncol = length(feature)*3)
  rownames(out) <- paste0(rep(durgL,each = length(timeL)-1), "_", timeL[-1])
  
  colnames(out) <- paste0(rep(feature, each = 3), c("_md", "_std", "_pvalue") )
  
  for(i in 1:length(durgL)){
    
    tmpphe <- phe[phe$Drug == durgL[i], ]
    tmpphe_match <- matchpairID(configdat = tmpphe, ID = "PID", Time = "Time", num = length(timeL))
    
    for(j in 1:length(feature)){
      for(z in 2:length(timeL)){      
        base <- tmpphe_match[tmpphe_match$Time == timeL[1], feature[j]]
        treat <- tmpphe_match[tmpphe_match$Time == timeL[z], feature[j]]
        value1 <- treat - base
        value2 <- round(median(value1, na.rm = T), 2)
        if(all(is.na(value1))){
          pvalue <- 1
        }else{
          pvalue <- wilcox.test(base, treat, paired= T)$p.value
        }
        #value3 <- round(mad(value1, na.rm = T), 2)
        value3 <- round(std_robust(value1), 2)
        
        if(timeN == 2){
          out[i+z-2, (3*j-2):(3*j)] <- c(value2, value3, pvalue)
        }else{
          out[2*i+z-3, (3*j-2):(3*j)] <- c(value2, value3, pvalue) 
        }      
        
        #out[i, (3*j-2):(3*j)] <- c(value2, value3, pvalue)
        
      }
    }
  }
  # colnames(out) <- paste0(rep(feature, each = 2), c("_md", "_std") )
  # for(i in 1:length(durgL)){
  #   
  #   tmpphe <- phe[phe$Drug == durgL[i], ]
  #   tmpphe_match <- matchpairID(configdat = tmpphe, ID = "PID", Time = "Time", num = length(timeL))
  #   
  #   for(j in 1:length(feature)){
  #     for(z in 2:length(timeL)){
  #       base <- tmpphe_match[tmpphe_match$Time == timeL[1], feature[j]]
  #       treat <- tmpphe_match[tmpphe_match$Time == timeL[z], feature[j]]
  #       value1 <- treat - base
  #       value2 <- round(median(value1, na.rm = T), 2)
  #       value3 <- round(std(value1), 2)
  #       if(timeN == 2){
  #         out[i+z-2, (2*j-1):(2*j)] <- c(value2, value3)
  #       }else{
  #         out[2*i+z-3, (2*j-1):(2*j)] <- c(value2, value3) 
  #       }
  #     }
  #   }
  # }
  return(out)
}





metainf_robust <- function(phe, feature = NULL){
  
  durgL <- unique(phe$Drug)
  timeL <- sort(unique(phe$Time))
  timeN <- length(timeL)
  #feature <-  c("richness", "shannon", "Uniqueness_bray", "Uniqueness_kendall",
  #              "Uniqueness_hellinger", "Uniqueness_JSD", "Uniqueness_spearman")
  
  if(is.null(feature)){
    stop("give some feature!")
  }
  out <- matrix(" ", nrow = length(durgL)*(length(timeL)-1), ncol = length(feature)*3)
  rownames(out) <- paste0(rep(durgL,each = length(timeL)-1), "_", timeL[-1])
  
  colnames(out) <- paste0(rep(feature, each = 3), c("_md", "_std", "_pvalue") )
  
  for(i in 1:length(durgL)){
    
    tmpphe <- phe[phe$Drug == durgL[i], ]
    tmpphe_match <- matchpairID(configdat = tmpphe, ID = "PID", Time = "Time", num = length(timeL))
    
    for(j in 1:length(feature)){
      for(z in 2:length(timeL)){      
        base <- tmpphe_match[tmpphe_match$Time == timeL[1], feature[j]]
        treat <- tmpphe_match[tmpphe_match$Time == timeL[z], feature[j]]
        value1 <- treat - base
        value2 <- round(median(value1, na.rm = T)/sd_robust(value1), 2)
        if(all(is.na(value1))){
          pvalue <- 1
        }else{
          pvalue <- wilcox.test(base, treat, paired= T)$p.value
        }
        #value3 <- round(mad(value1, na.rm = T), 2)
        value3 <- round(std_robust(value1), 2)
        
        if(timeN == 2){
          out[i+z-2, (3*j-2):(3*j)] <- c(value2, value3, pvalue)
        }else{
          out[2*i+z-3, (3*j-2):(3*j)] <- c(value2, value3, pvalue) 
        }      
        
        #out[i, (3*j-2):(3*j)] <- c(value2, value3, pvalue)
        
      }
    }
  }
  return(out)
}
