#' Calculating SNR value; Plotting a PCA panel
#'
#' @param expr_dt A expression profile (at protein level)
#' @param meta_dt A metadata file
#' @param output_dir A directory of the output file(s)
#' @import data.table
#' @importFrom dplyr %>%
#' @import ggplot2
#' @export

qc_snr <- function(expr_dt, meta_dt, output_dir = NULL) {

  expr_df <- data.frame(expr_dt[,2:ncol(expr_dt)],row.names = expr_dt[, 1])
  expr_df[is.na(expr_df)] <- 0

  expr_df_t <- t(expr_df)
  IDs <- rownames(expr_df_t)
  group <- meta_dt$sample
  IDs.group.mat<-data.table(IDs=IDs,group=group)

  pca_prcomp <- prcomp(expr_df_t,retx=T,scale. = T)
  pcs <- as.data.frame(predict(pca_prcomp))

  pcs$Sample.ID <- rownames(pcs)
  pcs$Sample <- meta_dt$sample

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

  signoise_db <- round(10*log10(signoise),3)

  colors.sample.Quartet<- c(
    'D5' = '#4CC3D9','D6' = '#7BC8A4','F7' = '#FFC65D','M8' = '#F16745')

  scale.axis.x <- c( min(pcs$PC1), max(pcs$PC1))
  scale.axis.y <- c( min(pcs$PC2), max(pcs$PC2))

  p<-ggplot(pcs, aes(x=PC1, y=PC2))+
    geom_point(aes(color=Sample),size=8)+
    # geom_encircle(aes(color=Sample), expand=0.02, alpha=.8)+
    theme_few()+
    labs(x = sprintf("PC1(%.2f%%)", summary(pca_prcomp)$importance[2,1]*100),
         y = sprintf("PC2(%.2f%%)", summary(pca_prcomp)$importance[2,2]*100),
         title=paste("SNR = ",signoise_db,sep=""),
         subtitle = paste("(Number of proteins = ",nrow(expr_dt),')',sep=""))+
    theme(axis.title = element_text(family="Arial",size=16,
                                    face="plain",color = "black"),
          axis.text = element_text(family="Arial",size=16,
                                   face="plain",color = "black"),
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

  output <- data.table(pcs[,c((ncol(pcs)-1):(ncol(pcs)), 1:(ncol(pcs)-2))])

  if(!is.null(output_dir)){
    output_dir_final1 <- file.path(output_dir,'pca_plot.png')
    ggsave(output_dir_final1,p,width = 6,height = 5.5)
    output_dir_final2 <- file.path(output_dir,'pca_table.tsv')
    write.table(output,output_dir_final2,sep = '\t',row.names = F)
  }

  return(list(plot = p, table = output, SNR = signoise_db))
}
