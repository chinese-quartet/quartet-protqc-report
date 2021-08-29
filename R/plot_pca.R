#' Plotting a PCA panel
#'
#' @param expr_dt_path A file path of the expression table file
#' @param meta_dt_path A file path of the metadata file
#' @param output_dir A directory of the output file(s)
#' @import data.table
#' @importFrom dplyr %>%
#' @import ggplot2
#' @export

plot_pca <- function(expr_dt_path,meta_dt_path,output_dir){

  ref_snrcorr_dir <- paste(system.file(package = "ProtQC"), "/data/ref_snrcorr.rds", sep = "")
  snrcorr <- readRDS(ref_snrcorr_dir)

  expr_dt_path

  expr_dt <- read.csv(expr_dt_path)
  meta_dt <- read.csv(meta_dt_path)

  expr_df <- data.frame(expr_dt[,2:ncol(expr_dt)],row.names = expr_dt$rowname)
  expr_df_t<-t(na.omit(expr_df))

  IDs <- rownames(expr_df_t)
  group <- meta_dt$sample
  IDs.group.mat<-data.table(IDs=IDs,group=group)

  pca_prcomp <- prcomp(expr_df_t,retx=T,scale. = T)
  pcs <- as.data.frame(predict(pca_prcomp))

  pcs$Sample_id <- rownames(pcs)
  pcs$group <- meta_dt$sample

  dt.perc.pcs <- data.table(
    PCX=1:nrow(pcs),
    Percent=summary(pca_prcomp)$importance[2,],
    AccumPercent=summary(pca_prcomp)$importance[3,])

  dt.dist <- data.table(
    ID.A = rep(IDs,each=length(IDs)),
    ID.B = rep(IDs,time=length(IDs)))

  dt.dist$group.A <- IDs.group.mat[match(dt.dist$ID.A,IDs.group.mat$IDs)]$group
  dt.dist$group.B <- IDs.group.mat[match(dt.dist$ID.B,IDs.group.mat$IDs)]$group

  dt.dist[,Type:=ifelse(ID.A==ID.B,'Same',
                        ifelse(group.A==group.B,'Intra','Inter'))]

  dt.dist[,Dist:=(dt.perc.pcs[1]$Percent*(pcs[ID.A,1]-pcs[ID.B,1])^2+
                    dt.perc.pcs[2]$Percent*(pcs[ID.A,2]-pcs[ID.B,2])^2)]

  dt.dist.stats <- dt.dist[,.(Avg.Dist=mean(Dist)),by=.(Type)]
  setkey(dt.dist.stats,Type)
  signoise <- dt.dist.stats['Inter']$Avg.Dist/dt.dist.stats['Intra']$Avg.Dist

  signoise_db <- round(10*log10(signoise),2)
  signoise_db_rank_length <- length(rank(c(signoise_db,snrcorr$SNR)))
  signoise_db_rank <- c(signoise_db_rank_length-rank(c(signoise_db,snrcorr$SNR))[1]+1)

  output_signoise_db <- data.table(
    "Quality Metrics" = c("Signal-to-Noise Ratio (SNR)"),
    Value = c(signoise_db),
    "Historical value(mean ± SD)" = c('21.64 ± 6.21'),
    Rank = c(paste(as.character(signoise_db_rank),'/',signoise_db_rank_length,sep = ''))
  )

  colors.sample.Quartet<- c('D5' = '#4CC3D9','D6' = '#7BC8A4','F7' = '#FFC65D','M8' = '#F16745')

  scale.axis.x <- c( min(pcs$PC1), max(pcs$PC1))
  scale.axis.y <- c( min(pcs$PC2), max(pcs$PC2))

  p<-ggplot(pcs, aes(x=PC1, y=PC2))+
    geom_point(aes(color=group),size=8)+
    # geom_encircle(aes(color=group), expand=0.02, alpha=.8)+
    theme_few()+
    labs(x = sprintf("PC1(%.2f%%)", summary(pca_prcomp)$importance[2,1]*100),
         y = sprintf("PC2(%.2f%%)", summary(pca_prcomp)$importance[2,2]*100),
         title=paste("SNR = ",signoise_db,sep=""),
         subtitle = paste("(N = ",nrow(expr_dt),')',sep=""))+
    theme(axis.title = element_text(family="Arial",size=16,face="plain",color = "black"),
          axis.text = element_text(family="Arial",size=16,face="plain",color = "black"),
          title = element_text(size = 12,hjust=0.5),
          legend.title = element_text(size=16),
          legend.text = element_text(size=16,color='gray40'),
          plot.subtitle = element_text(hjust=0.5,size=16),
          plot.title = element_text(hjust=0.5,size=16))+
    scale_color_manual(values = colors.sample.Quartet)+
    scale_fill_manual(values = colors.sample.Quartet)+
    scale_x_continuous(limits=c(1.1*scale.axis.x[1], 1.1*scale.axis.x[2]))+
    scale_y_continuous(limits=c(1.1*scale.axis.y[1], 1.1*scale.axis.y[2]))+
    guides(colour = guide_legend(override.aes = list(size=2)))+
    guides(shape=guide_legend(override.aes = list(size=3)))

  output <- data.table(
    sample_id = meta_dt$name,
    group = meta_dt$sample,
    pcs[,1:(ncol(pcs)-2)]
  )

  output_dir_final1 <- paste(output_dir,'pca_plot.png',sep = '')
  ggsave(output_dir_final1,p,width = 6,height = 5.5)

  output_dir_final2 <- paste(output_dir,'pca_table.csv',sep = '')
  write.csv(output,output_dir_final2,row.names = F)

  return(output_signoise_db)
}
