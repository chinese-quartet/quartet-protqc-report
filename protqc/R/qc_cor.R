#' Calculating RC value; Plotting a scatterplot
#'
#' @param expr_dt A expression table file (at peptide level)
#' @param meta_dt A metadata file
#' @param output_dir A directory of the output file(s)
#' @import ggplot2
#' @import data.table
#' @importFrom dplyr %>%
#' @importFrom ggthemes theme_few
#' @importFrom edgeR DGEList
#' @importFrom edgeR filterByExpr
#' @importFrom edgeR calcNormFactors
#' @importFrom limma voom
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @export

DEPanalysis <- function(expr, group){

  dge <- DGEList(counts = expr)
  design <- model.matrix(~group)

  keep <- filterByExpr(dge, design)
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)

  v <- voom(dge, design, plot=F)
  fit <- lmFit(v, design)

  fit <- eBayes(fit)
  result <- topTable(fit, coef=ncol(design), sort.by = 'logFC', number = Inf)
  result$Sequence = rownames(result)
  result$Sequence.Number = nrow(result)
  result$Sample1 =  levels(group)[1]
  result$Sample2 =  levels(group)[2]
  result$Sample.Pair = paste(levels(group)[2], levels(group)[1], sep = '/')
  return(result)
}

qc_cor <- function(expr_dt, meta_dt, output_dir = NULL, show_sample_pairs = F){

  ref_dt_dir <- file.path(
    system.file(package = "protqc"), "data/reference_dataset.rds")
  ref_dt <- readRDS(ref_dt_dir)

  group <- factor(meta_dt$sample)
  expr_matrix <- as.matrix(data.frame(
    expr_dt[,2:ncol(expr_dt)],row.names = expr_dt[, 1]))
  expr_matrix[is.na(expr_matrix)] <- 0

  result_final <- c()
  samples <- c('D5','D6','F7','M8')
  for(j in c(1,3,4)){
    sample_pair <- paste(samples[j], samples[2], sep = '/')
    ref_tmp <- ref_dt[ref_dt$Sample.Pair %in% sample_pair, ]

    col1<-grep(samples[j],colnames(expr_matrix))
    col2<-grep(samples[2],colnames(expr_matrix))

    grouped_expr <- as.matrix(expr_matrix[,c(col1,col2)])
    grouped_expr <- grouped_expr[
      apply(grouped_expr, 1, function(x) length(which(x == 0)) < 3), ]
    grouped_expr <- grouped_expr[rownames(grouped_expr) %in% ref_tmp$Sequence,]

    sample_pairs <- factor(rep(samples[c(j, 2)], each = 3),
                           levels = c(samples[2], samples[j]), ordered = T)
    result_tmp <- DEPanalysis(grouped_expr,group = sample_pairs)
    result_tmp <- result_tmp[result_tmp$adj.P.Val < 0.05, ]

    result_final <- rbind(result_final, result_tmp)
  }
  result_final <- as.data.table(result_final)
  result_trim <- result_final[, c(7, 11, 1)]
  result_trim$name <- apply(result_trim, 1, function(x) paste(x[1],x[2]))
  ref_dt$name <- apply(ref_dt, 1, function(x) paste(x[1],x[2]))
  result_withref <- merge(result_trim,ref_dt,by = 'name')

  df_test <- data.frame(result_withref[, c(1,2,3,4,7)])
  colnames(df_test) <- c(
    'Name', 'Sequence', 'Sample.Pair', 'logFC.Test', 'logFC.Reference')
  cor_value <- round(cor(df_test$logFC.Test, df_test$logFC.Reference), 3)

  p <- ggplot(df_test,aes(x=logFC.Reference, y=logFC.Test))+
    theme_few()+
    theme(
      axis.text = element_text(family="Arial",size=16,
                               face="plain",color = "black"),
      axis.title = element_text(family="Arial",size=16,
                                face='plain',color = "black"),
      legend.text = element_text(family="Arial",size=16,
                                 color = "black"),
      plot.subtitle = element_text(family="Arial",size=16,
                                   face='plain',color = "black",hjust=0.5),
      plot.title = element_text(family="Arial",size=16,
                                face='plain',color = "black",hjust=0.5)
    )+
    labs(
      y = 'Test Dataset', x = 'Reference Datasets',
      title=paste("Correlation = ",cor_value,sep=""),
      subtitle = paste("(Number of peptides = ",nrow(df_test),')',sep="")
    )+
    coord_fixed(
      xlim = c(-max(df_test$logFC.Reference),max(df_test$logFC.Reference)),
      ylim = c(-max(df_test$logFC.Reference),max(df_test$logFC.Reference))
    )

  if(show_sample_pairs == T) {
    colors.sample.Quartet<- c(
      'D5/D6' = '#4CC3D9', 'F7/D6' = '#FFC65D', 'M8/D6' = '#F16745')
    p <- p + geom_point(aes(color = Sample.Pair),size=2.5, alpha=.5) +
      scale_color_manual(values = colors.sample.Quartet) +
      scale_fill_manual(values = colors.sample.Quartet)
  }else {
    p <- p + geom_point(color = 'steelblue4',size=2.5, alpha=.1)
  }

  if(!is.null(output_dir)){
    output_dir_final1 <- file.path(output_dir,'corr_plot.png')
    output_dir_final2 <- file.path(output_dir,'deps_table.tsv')
    output_dir_final3 <- file.path(output_dir,'corr_table.tsv')
    ggsave(output_dir_final1,p,height = 5.5,width = 5.5)
    write.table(result_final,output_dir_final2,sep = '\t',row.names = F)
    write.table(df_test,output_dir_final3,sep = '\t',row.names = F)
  }

  return(list(plot = p, DEPs = result_final,
              logfc = df_test, COR = cor_value))
}
