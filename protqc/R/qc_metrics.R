#' Statistics for basic information
#' @param expr_dt A expression profile
#' @param meta_dt A metadata file
#' @import stats
#' @importFrom psych corr.test
#' @importFrom reshape2 melt
#' @export

qc_info <- function(expr_dt, meta_dt) {

  # Load data --------------------------------
  m <- meta_dt[, colnames(meta_dt) %in% c("library", "sample")]
  d <- expr_dt
  s <- meta_dt$sample

  # Replace zero by NA -----------------------
  d[d == 0] <- NA

  # Statistics: number of features -----------
  uniq_pro <- unique(d[, 1])
  stat_num <- length(uniq_pro)

  # Statistics: missing percentage -----------
  d_all_num <- nrow(d) * (ncol(d) - 1)
  d_missing_num <- length(which(is.na(d)))
  prop_missing <- d_missing_num * 100 / d_all_num
  stat_missing <- round(prop_missing, 3)

  # Check if replicates available ------------
  samples <- table(s)
  rep_samples <- samples[samples > 1]
  rep_num <- length(rep_samples)
  if (rep_num == 0) {
    stop("No replicates are available.")
  } else {
    # Calculating: absolute correlation ------
    d_mtx <- d[, 2:ncol(d)]
    d_cortest <- corr.test(d_mtx, method = "pearson", adjust = "fdr")
    d_pmtx <- d_cortest$p
    d_cormtx <- d_cortest$r
    d_cormtx[d_pmtx > 0.05] <- 0
    d_cordf <- melt(d_cormtx)
    d_cordf <- d_cordf[d_cordf$Var2 != d_cordf$Var1, ]
    d_cordf <- merge(d_cordf, m, by.x = "Var1", by.y = "library")
    d_cordf <- merge(d_cordf, m, by.x = "Var2", by.y = "library")
    d_cordf <- d_cordf[d_cordf$sample.x == d_cordf$sample.y, ]
    cor_value <- median(d_cordf$value)
    stat_acor <- round(cor_value, 3)
  }

  # Calculating: CV --------------------------
  d_long <- melt(d)
  d_long <- na.omit(d_long)
  if (length(d_long$value[d_long$value < 0])) {
    d_long$value <- 2 ^ (d_long$value)
    message("All values were squared to avoid negative values.")
  }
  d_long <- merge(d_long, m, by.x = "variable", by.y = "library")
  colnames(d_long) <- c("library", "feature", "value", "sample")
  d_cv <- aggregate(
    value ~ feature + sample, data = d_long,
    FUN = function(x) sd(x) / mean(x))
  stat_cv <- round(median(d_cv$value, na.rm = T) * 100, 3)

  # Output -----------------------------------
  stat_all <- c(stat_num, stat_missing, stat_acor, stat_cv)

  return(stat_all)
}

#' Calculating SNR value; Plotting a PCA panel
#' @param expr_dt A expression profile (at protein level)
#' @param meta_dt A metadata file
#' @param output_dir A directory of the output file(s)
#' @param plot if True, a plot will be output.
#' @import stats
#' @import utils
#' @importFrom rlang :=
#' @importFrom data.table data.table
#' @importFrom data.table setkey
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 ggsave
#' @importFrom ggthemes theme_few
#' @export

qc_snr <- function(expr_dt, meta_dt, output_dir=NULL, plot=TRUE) {

  # Load data --------------------------------------
  expr_ncol <- ncol(expr_dt)
  expr_df <- data.frame(expr_dt[, 2:expr_ncol], row.names = expr_dt[, 1])

  # Replace NA by zero -----------------------------
  expr_df[is.na(expr_df)] <- 0

  # Label the grouping info ------------------------
  ids <- colnames(expr_df)
  group <- meta_dt$sample
  ids_group_mat <- data.table(id = ids, group = group)

  # PCA --------------------------------------------
  expr_df_t <- t(expr_df)
  pca_prcomp <- prcomp(expr_df_t, retx = T, scale. = T)
  pcs <- as.data.frame(predict(pca_prcomp))
  pcs$sample_id <- rownames(pcs)
  pcs$sample <- meta_dt$sample

  # Calculating: SNR -------------------------------
  dt_perc_pcs <- data.table(PCX = 1:nrow(pcs),
                            Percent = summary(pca_prcomp)$importance[2, ],
                            AccumPercent = summary(pca_prcomp)$importance[3, ])

  dt_dist <- data.table(id_a = rep(ids, each = length(ids)),
                        id_b = rep(ids, time = length(ids)))

  dt_dist$group_a <- ids_group_mat[match(dt_dist$id_a, ids_group_mat$id)]$group
  dt_dist$group_b <- ids_group_mat[match(dt_dist$id_b, ids_group_mat$id)]$group

  dt_dist[, type := ifelse(id_a == id_b, "Same",
                           ifelse(group_a == group_b, "Intra", "Inter"))]

  dt_dist[, dist := (dt_perc_pcs[1]$Percent * (pcs[id_a, 1] - pcs[id_b, 1])^2 +
                      dt_perc_pcs[2]$Percent * (pcs[id_a, 2] - pcs[id_b, 2])^2)]

  dt_dist_stats <- dt_dist[, list(avg_dist = mean(dist)), by = list(type)]
  setkey(dt_dist_stats, type)
  signoise <- dt_dist_stats["Inter"]$avg_dist / dt_dist_stats["Intra"]$avg_dist
  signoise_db <- round(10 * log10(signoise), 3)

  # Plot -------------------------------------------
  if (plot) {
    colors_custom <- c("D5" = "#4CC3D9",
                     "D6" = "#7BC8A4",
                     "F7" = "#FFC65D",
                     "M8" = "#F16745")
    text_custom_theme <- element_text(size = 16,
                                      face = "plain",
                                      color = "black",
                                      hjust = 0.5)
    scale_axis_x <- c(min(pcs$PC1), max(pcs$PC1))
    scale_axis_y <- c(min(pcs$PC2), max(pcs$PC2))

    pc1_prop <- summary(pca_prcomp)$importance[2, 1]
    pc2_prop <- summary(pca_prcomp)$importance[2, 2]
    text_axis_x <- sprintf("PC1(%.2f%%)", pc1_prop * 100)
    text_axis_y <- sprintf("PC2(%.2f%%)", pc1_prop * 100)
    limit_x <- c(1.1 * scale_axis_x[1], 1.1 * scale_axis_x[2])
    limit_y <- c(1.1 * scale_axis_y[1], 1.1 * scale_axis_y[2])

    p_title <- paste("SNR = ", signoise_db, sep = "")
    p_subtitle <- paste("(Number of proteins = ", nrow(expr_dt), ")", sep = "")
    p <- ggplot(pcs, aes(x = .data$PC1, y = .data$PC2)) +
      geom_point(aes(color = sample), size = 8) +
      theme_few() +
      theme(plot.title = text_custom_theme,
            plot.subtitle = text_custom_theme,
            axis.title = text_custom_theme,
            axis.text = text_custom_theme,
            legend.title = text_custom_theme,
            legend.text = element_text(size = 16, color = "gray40")) +
      labs(x = text_axis_x,
          y = text_axis_y,
          title = p_title,
          subtitle = p_subtitle) +
      scale_color_manual(values = colors_custom) +
      scale_x_continuous(limits = limit_x) +
      scale_y_continuous(limits = limit_y) +
      guides(colour = guide_legend(override.aes = list(size = 2))) +
      guides(shape = guide_legend(override.aes = list(size = 3)))
  }

  pc_num <- ncol(pcs)
  output <- data.table(pcs[, c((pc_num - 1):pc_num, 1:(pc_num - 2))])

  # Save & Output ----------------------------------
  if (!is.null(output_dir)) {
    if (plot) {
      output_dir_final1 <- file.path(output_dir, "pca_plot.png")
      ggsave(output_dir_final1, p, width = 6, height = 5.5)
    }
    output_dir_final2 <- file.path(output_dir, "pca_table.tsv")
    write.table(output, output_dir_final2, sep = "\t", row.names = F)
  }

  return(list(table = output, SNR = signoise_db))
}

#' Analysis: differential expression
#' @param expr A expression table file (at peptide level)
#' @param group The grouping info
#' @import stats
#' @importFrom edgeR DGEList
#' @importFrom edgeR filterByExpr
#' @importFrom edgeR calcNormFactors
#' @importFrom limma voom
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @export

dep_analysis <- function(expr, group) {

  dge <- DGEList(counts = expr)
  design <- model.matrix(~ group)

  keep <- filterByExpr(dge, design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)

  v <- voom(dge, design, plot = F)
  fit <- lmFit(v, design)

  fit <- eBayes(fit)
  result <- topTable(fit, coef = ncol(design), sort.by = "logFC", number = Inf)
  result$Sequence <- rownames(result)
  result$Sequence.Number <- nrow(result)
  result$Sample1 <- levels(group)[1]
  result$Sample2 <- levels(group)[2]
  result$Sample.Pair <- paste(levels(group)[2], levels(group)[1], sep = "/")

  return(result)

}

#' Calculating RC value; Plotting a scatterplot
#' @param expr_dt A expression table file (at peptide level)
#' @param meta_dt A metadata file
#' @param output_dir A directory of the output file(s)
#' @param plot if True, a plot will be output.
#' @param show_sample_pairs if True, samples in plot will be labeled.
#' @import stats
#' @import utils
#' @importFrom rlang .data
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggthemes theme_few
#' @importFrom data.table as.data.table
#' @importFrom ggplot2 ggsave
#' @export

qc_cor <- function(expr_dt, meta_dt,
                   output_dir=NULL, plot=FALSE, show_sample_pairs=FALSE) {

  # Load data ------------------------------------------------------
  load(system.file("data/reference_dataset.rda", package = "protqc"))
  ref_dt <- reference_dataset
  expr_ncol <- ncol(expr_dt)
  expr_df <- data.frame(expr_dt[, 2:expr_ncol], row.names = expr_dt[, 1])
  expr_matrix <- as.matrix(expr_df)

  # Replace NA by zero ---------------------------------------------
  expr_matrix[is.na(expr_matrix)] <- 0

  # Check the grouping info ----------------------------------------
  samples <- unique(meta_dt$sample)
  if (length(samples) > 1) {
    check_d6 <- "D6" %in% samples
    if (!check_d6) {
      stop("No D6 samples are available.")
    } else {
      samples <- c("D6", samples[!samples %in% "D6"])
      pair_num <- length(samples[!samples %in% "D6"])
    }
  } else {
    stop("No grouping info.")
  }

  # Analysis: Differential expression ------------------------------
  result_final <- c()
  for (j in 2:(pair_num + 1)) {
    sample_pair <- paste(samples[j], "D6", sep = "/")
    ref_tmp <- ref_dt[ref_dt$Sample.Pair %in% sample_pair, ]

    col1 <- grep(samples[j], colnames(expr_matrix))
    col2 <- grep("D6", colnames(expr_matrix))

    e_tmp <- expr_matrix[, c(col1, col2)]
    e_tmp <- e_tmp[apply(e_tmp, 1, function(x) length(which(x == 0)) < 3), ]
    expr_grouped <- e_tmp[rownames(e_tmp) %in% ref_tmp$Sequence, ]

    sample_pairs <- factor(x = rep(c(samples[j], "D6"), each = 3),
                           levels = c("D6", samples[j]),
                           ordered = T)
    result_tmp <- dep_analysis(expr = expr_grouped, group = sample_pairs)
    result_tmp <- result_tmp[result_tmp$adj.P.Val < 0.05, ]

    result_final <- rbind(result_final, result_tmp)

  }

  # Calculating: RC -----------------------------------------------
  result_final <- as.data.table(result_final)
  result_trim <- result_final[, c(7, 11, 1)]
  result_trim$name <- apply(result_trim, 1, function(x) paste(x[1], x[2]))
  ref_dt$name <- apply(ref_dt, 1, function(x) paste(x[1], x[2]))
  result_withref <- merge(result_trim, ref_dt, by = "name")

  df_test <- data.frame(result_withref[, c(1, 2, 3, 4, 7)])
  colnames(df_test) <- c("Name", "Sequence", "Sample.Pair",
                         "logFC.Test", "logFC.Reference")

  cor_value <- cor(x = df_test$logFC.Test, y = df_test$logFC.Reference)
  cor_value <- round(cor_value, 3)

  # Plot ----------------------------------------------------------
  if (plot) {
    text_custom_theme <- element_text(size = 16,
                                      face = "plain",
                                      color = "black",
                                      hjust = 0.5)

    scale_axis_r <- c(min(df_test$logFC.Reference),
                      max(df_test$logFC.Reference))
    scale_axis_t <- c(min(df_test$logFC.Test),
                      max(df_test$logFC.Test))
    limit <- max(abs(c(scale_axis_r, scale_axis_t)))
    limit_axis <- c(- limit, limit)

    plot_title <- paste("RC = ", cor_value, sep = "")
    plot_subtitle <- paste("(Number of peptides = ", nrow(df_test), ")", sep = "")

    p <- ggplot(df_test, aes(x = .data$logFC.Reference, y = .data$logFC.Test)) +
      theme_few() +
      theme(plot.title = text_custom_theme,
            plot.subtitle = text_custom_theme,
            axis.title = text_custom_theme,
            axis.text = text_custom_theme,
            legend.title = text_custom_theme,
            legend.text = element_text(size = 16, color = "gray40")) +
      labs(y = "log2FC (Test Dataset)",
          x = "log2FC (Reference Datasets)",
          title = plot_title,
          subtitle = plot_subtitle) +
      coord_fixed(xlim = limit_axis, ylim = limit_axis)

    if (show_sample_pairs == T) {
      colors_custom <- c("D5/D6" = "#4CC3D9",
                        "F7/D6" = "#FFC65D",
                        "M8/D6" = "#F16745")
      p <- p +
        geom_point(aes(color = .data$Sample.Pair), size = 2.5, alpha = .5) +
        scale_color_manual(values = colors_custom)
    }else {
      p <- p + geom_point(color = "steelblue4", size = 2.5, alpha = .1)
    }

  }

  # Save & Output -------------------------------------------------
  if (!is.null(output_dir)) {
    if (plot) {
      output_dir_final1 <- file.path(output_dir, "corr_plot.png")
      ggsave(output_dir_final1, p, height = 5.5, width = 5.5)
    }
    output_dir_final2 <- file.path(output_dir, "deps_table.tsv")
    output_dir_final3 <- file.path(output_dir, "corr_table.tsv")
    write.table(result_final, output_dir_final2, sep = "\t", row.names = F)
    write.table(df_test, output_dir_final3, sep = "\t", row.names = F)
  }

  output_list <- list(DEPs = result_final,
                      logfc = df_test,
                      COR = cor_value)

  return(output_list)
}
