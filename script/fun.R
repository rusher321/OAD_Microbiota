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
