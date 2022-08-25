pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "SingleCellExperiment", 
          "vroom", "jhtools", "glue", "openxlsx", "ggsci", "patchwork", "cowplot", "tidyverse", "dplyr",
          "ComplexHeatmap", "viridis", "factoextra", "pheatmap", "RColorBrewer", "circlize",
          "jhuanglabHyperion")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
conflicted::conflict_prefer("mutate", "dplyr")
project <- "hyperion"
dataset <- "qzhang"
species <- "human"
workdir <- glue::glue("~/projects/{project}/analysis/{dataset}/{species}/workdir")
setwd(workdir)

#dir define
sce_dir <- "/cluster/home/jhuang/projects/hyperion/analysis/qzhang/human/steinbock/rds/combined"
fibscore_dir <- "/cluster/home/yjliu_jh/share"
sparse_dir <- "/cluster/home/yjliu_jh/share"

#data prepare
sce <- readr::read_rds(glue::glue("{sce_dir}/all_anno.rds"))
dot_colors <- c(Surgery_without_chemo = "#B61932", Paracancerous = "#01A0A6",
                Punc_liver = "#5589C2", Punc_pancreas = "#BCBDDC",
                Surgery_after_chemo = "#FFE475", Normal = "#5D9D58")


colors <- alpha(dot_colors, 0.6)
cld <- colData(sce) %>% as.data.frame()
cld$log2area <- log2(cld$area)
fibscore <- read_rds(glue::glue("{fibscore_dir}/fibscore_roi_new.rds")) %>% dplyr::select(sample_tiff_id, score_col)
fibscore_fil <- fibscore[fibscore$sample_tiff_id %in% cld$sample_tiff_id,] %>% dplyr::distinct()
score <- read_rds(glue::glue("{sparse_dir}/scores.rds")) %>% dplyr::select(TIFF, SparseScore)


#sfig2 updated code from @rujia zheng
df <- cld %>% dplyr::mutate(ROI_area = width_px*height_px) %>% group_by(sample_tiff_id) %>%
  dplyr::mutate(cell_num = n(), ROI_area = ROI_area) %>% mutate(cell_density = 100*cell_num/ROI_area) %>%
  dplyr::select(sample_tiff_id, cell_density, stype3) %>% dplyr::distinct() %>%
  left_join(fibscore_fil, by = c("sample_tiff_id" = "sample_tiff_id")) %>%
  left_join(score, by = c("sample_tiff_id" = "TIFF"))
df$SparseScore <- as.double(df$SparseScore)

axistheme1 <- theme(axis.title.y = element_text(size = 7),
                    axis.text.y = element_text(size = 7),
                    axis.title.x = element_blank(),
                    axis.ticks.x =  element_blank(),
                    axis.text.x = element_blank())
plot.theme <- theme_bw() + theme(panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank())
my_comparisons <- list(c("Normal","Paracancerous"), c("Normal", "Surgery_without_chemo"),
                       c("Paracancerous", "Surgery_without_chemo"),
                       c("Surgery_without_chemo", "Surgery_after_chemo"),
                       c("Punc_pancreas", "Punc_liver"))
cld$stype3<- factor(cld$stype3, levels = c("Surgery_without_chemo", "Surgery_after_chemo", "Punc_pancreas","Punc_liver", "Paracancerous", "Normal"))
df$stype3<- factor(df$stype3, levels = c("Surgery_without_chemo", "Surgery_after_chemo", "Punc_pancreas","Punc_liver", "Paracancerous", "Normal"))

p1 <- (ggplot(cld, aes_string(x = "stype3", y = "log2area")) +
         geom_boxplot(outlier.shape= NA, lwd = 0.2, color = "#000000", fill = dot_colors) +
         ylim(2.5, 14.5) +
         plot.theme +
         theme(legend.position = "none") +
         axistheme1+
         stat_compare_means(method = "wilcox.test",
                            comparisons = my_comparisons,
                            hide.ns = TRUE,
                            bracket.size = 0.2,
                            vjust = 0.8,
                            aes(label = ..p.signif..)))
p2 <- (ggplot(cld, aes_string(x = "stype3", y = "eccentricity")) +
         geom_boxplot(outlier.shape= NA, lwd = 0.2, color = "#000000", fill = dot_colors) +
         ylim(0.3, 1.6) +
         plot.theme +
         theme(legend.position = "none") +
         axistheme1+
         stat_compare_means(method = "wilcox.test",
                            comparisons = my_comparisons,
                            hide.ns = TRUE,
                            bracket.size = 0.2,
                            vjust = 0.8,
                            aes(label = ..p.signif..)))
p3 <- (ggplot(df, aes_string(x = "stype3", y = "cell_density")) +
         geom_boxplot(outlier.shape= NA, lwd = 0.2, color = "#000000", fill = dot_colors) +
         ylim(0, 1.9) +
         plot.theme +
         theme(legend.position = "none") +
         axistheme1+
         stat_compare_means(method = "wilcox.test",
                            comparisons = my_comparisons,
                            hide.ns = TRUE,
                            bracket.size = 0.2,
                            vjust = 0.8,
                            aes(label = ..p.signif..)))
p4 <- (ggplot(df, aes_string(x = "stype3", y = "score_col")) +
         geom_boxplot(outlier.shape= NA, lwd = 0.2, color = "#000000", fill = dot_colors) +
         ylim(-2, 8) +
         plot.theme +
         theme(legend.position = "none") +
         axistheme1+
         stat_compare_means(method = "wilcox.test",
                            comparisons = my_comparisons,
                            hide.ns = TRUE,
                            bracket.size = 0.2,
                            vjust = 0.8,
                            aes(label = ..p.signif..)))
p5 <- (ggplot(df, aes_string(x = "stype3", y = "SparseScore")) +
         geom_boxplot(outlier.shape= NA, lwd = 0.2, color = "#000000", fill = dot_colors) +
         ylim(0, 1) +
         plot.theme +
         theme(legend.position = "none") +
         axistheme1+
         stat_compare_means(method = "wilcox.test",
                            comparisons = my_comparisons,
                            hide.ns = TRUE,
                            bracket.size = 0.2,
                            vjust = 0.8,
                            aes(label = ..p.signif..)))
p <- list(p1, p2, p3, p4, p5)

multi_plot(fig_fn = "sfig2_fig1_stype3_area_ecce_boxplot_swilcox.pdf", width = 7, height = 1.5, p_list = p,
           ncol = 5, nrow = 1, tag_levels = NULL)
