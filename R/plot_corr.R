#' Plotting a scatterplot
#'
#' @param expr_dt_path A file path of the expression table file
#' @param meta_dt_path A file path of the metadata file
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

plot_corr <- function(expr_dt_path,meta_dt_path,output_dir){

  expr_dt <- read.csv(expr_dt_path)
  meta_dt <- read.csv(meta_dt_path)

  ref_snrcorr_dir <- paste(system.file(package = "ProtQC"), "/data/ref_snrcorr.rds", sep = "")
  snrcorr <- readRDS(ref_snrcorr_dir)

  ref_dt_dir <- paste(system.file(package = "ProtQC"), "/data/example_ref_dt.rds", sep = "")
  ref_dt <- readRDS(ref_dt_dir)

  group <- factor(meta_dt$sample)
  expr_matrix <- as.matrix(data.frame(expr_dt[,2:ncol(expr_dt)],row.names = expr_dt$rowname))
  expr_matrix[is.na(expr_matrix)] <- 0

  result_final <- c()
  sample_pairs <- combn(levels(group),2)
  for(i in 1:ncol(sample_pairs)){

    v1 <- expr_matrix[,grep(sample_pairs[1,i],group)]
    v2 <- expr_matrix[,grep(sample_pairs[2,i],group)]

    dge <- DGEList(counts = cbind(v1,v2))
    sample_pair <- factor(rep(1:2,c(ncol(v1),ncol(v2))))
    design <- model.matrix(~sample_pair)

    keep <- filterByExpr(dge, design)
    dge <- dge[keep,,keep.lib.sizes=FALSE]
    dge <- calcNormFactors(dge)

    v <- voom(dge, design, plot=F)
    fit <- lmFit(v, design)

    fit <- eBayes(fit)
    result <- topTable(fit, coef=ncol(design), sort.by = 'logFC', number = Inf)
    result$gene = rownames(result)
    result$sample_pair1 =  sample_pairs[1,i]
    result$sample_pair2 =  sample_pairs[2,i]
    result$sample_pair = paste(sample_pairs[2,i],sample_pairs[1,i],sep = '.')

    result_final <-rbind(result_final,result)
  }
  test_dt <- as.data.table(result_final)

  sample_pairs <- unique(test_dt$sample_pair)

  df_cor <- c()
  for(i in 1:length(sample_pairs)){
    df <- test_dt[test_dt$sample_pair %in% sample_pairs[i],]
    df_perpair <- df[df$adj.P.Val<0.05,c(7,1)]

    df_ref <- ref_dt[ref_dt$sample_pair %in% sample_pairs[i],]
    df_ref_perpair <- df_ref[df_ref$adj.P.Val<0.05,c(7,1)]

    df_test_perpair <- merge(df_perpair,df_ref_perpair,by.x = 'gene',by.y = 'gene')

    REC_value <- cor(df_test_perpair$logFC.x,df_test_perpair$logFC.y)

    df_cor <- rbind(df_cor,data.frame(sample_pair = sample_pairs[i],REC = REC_value))
  }

  sample_pair <- df_cor[order(df_cor$REC)[3],1]
  cor_value <- df_cor[order(df_cor$REC)[3],2]
  cor_value <- round(cor_value,3)

  cor_value_rank_length <- length(rank(c(cor_value,snrcorr$Corr)))
  cor_value_rank <- c(cor_value_rank_length-rank(c(cor_value,snrcorr$Corr))[1]+1)

  output_cor_value <- data.table(
    "Quality Metrics" = c("Correlation with Reference Datasets"),
    Value = c(cor_value),
    "Historical value(mean ± SD)" = c('0.698 ± 0.302'),
    Rank = c(paste(as.character(cor_value_rank),'/',cor_value_rank_length,sep = ''))
  )

  test_dt <- data.frame(test_dt)
  df <- test_dt[test_dt$sample_pair %in% sample_pair,]
  df_perpair <- df[df$adj.P.Val<0.05,c(7,1)]

  ref_dt <- data.frame(ref_dt)
  df_ref <- ref_dt[ref_dt$sample_pair %in% sample_pair,]
  df_ref_perpair <- df_ref[df_ref$adj.P.Val<0.05,c(7,1)]

  df_test_perpair <- merge(df_perpair,df_ref_perpair,by.x = 'gene',by.y = 'gene')

  p <- ggplot(df_test_perpair,aes(x=logFC.x, y=logFC.y))+
    geom_point(color='steelblue4', size=2.5, alpha=.1)+
    scale_fill_brewer(palette = 'Blues')+
    theme_few()+
    theme(axis.text = element_text(family="Arial",size=16,face="plain",color = "black"),
          axis.title = element_text(family="Arial",size=16,face='plain',color = "black"),
          legend.text = element_text(family="Arial",size=16,color = "black"),
          plot.subtitle = element_text(family="Arial",face='plain',color = "black",hjust=0.5,size=16),
          plot.title = element_text(family="Arial",face='plain',color = "black",hjust=0.5,size=16))+
    labs(x='Reference Datasets',
         y='Test Dataset',
         title=paste("Correlation = ",cor_value,sep=""),
         subtitle = paste("(N = ",nrow(df_test_perpair),')',sep=""))+
    coord_fixed(
      xlim = c(-max(df_test_perpair$logFC.y),max(df_test_perpair$logFC.y)),
      ylim = c(-max(df_test_perpair$logFC.y),max(df_test_perpair$logFC.y)))

  output_dir_final1 <- paste(output_dir,'corr_plot.png',sep = '')
  output_dir_final2 <- paste(output_dir,'deps_table.tsv',sep = '')
  output_dir_final3 <- paste(output_dir,'corr_table.tsv',sep = '')

  ggsave(output_dir_final1,p,height = 5.5,width = 5.5)
  write.csv(test_dt,output_dir_final2,row.names = F)
  write.csv(df_test_perpair,output_dir_final3,row.names = F)

  return(output_cor_value)
}
