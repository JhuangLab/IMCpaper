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

#data prepare
sce <- readr::read_rds(glue::glue("{sce_dir}/all_anno.rds"))

# fig1a updated code from @rujia zheng ------------------------------------
sp_info <- colData(sce) |> as.data.frame() |> dplyr::select(sample_id, stype3) |> 
  distinct() |> remove_rownames() |> mutate(group = ifelse(stype3 != "Paracancerous" & stype3 != "Normal", "Tumor", stype3))

color <- c(Surgery_without_chemo = "#B61932", Paracancerous = "#01A0A6", 
           Punc_liver = "#5589C2", Punc_pancreas = "#BCBDDC", 
           Surgery_after_chemo = "#FFE475", Normal = "#5D9D58")
sp_info$group <- factor(sp_info$group, levels = c("Tumor", "Paracancerous", "Normal"))
sp_info$stype3 <- factor(sp_info$stype3, levels = c("Punc_liver", "Punc_pancreas", 
                                                    "Surgery_after_chemo", "Surgery_without_chemo",
                                                    "Paracancerous", "Normal"))
# plot
p <- ggplot(sp_info, aes(x = group, fill = stype3)) +
  geom_bar(position = "stack", width = 0.8) +
  scale_fill_manual(values = color) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank()) +
  geom_text(aes(label = ..count..), stat = "count", position = position_stack(vjust = 0.5))
ggsave("sample_info_fig1a.pdf", width = 4, height = 6)

set.seed(2022)
stype3_color <- c(Surgery_without_chemo = "#B61932", Paracancerous = "#01A0A6",
                  Punc_liver = "#5589C2", Punc_pancreas = "#BCBDDC",
                  Surgery_after_chemo = "#FFE475", Normal = "#5D9D58")

sce_sub <- subset_hyp(sce, cell_number = 300000)

png("tsne_stype3_300000_fig1a.png", res = 300, width = 7, height = 5, units = 'in')
print(plot2d(sce_sub, color.by = "stype3", show.cluser.id = F, size = 0.01,
             item.use = c("tSNE1", "tSNE2")) + scale_color_manual(values = stype3_color)) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
dev.off()

# fig1d updated code from @rujia zheng ------------------------------------
meta_all <- metadata(sce)$sampleinfo
coldata <- as.data.frame(colData(sce))
c_mat_lst <- jhuanglabHyperion::fetch_heatmap_dat(sce, cluster.id = "cell_type10")
c_mat <- c_mat_lst[[1]]
cluster_prop <- c_mat_lst[[2]]

# cell grouping and ordering data for heatmap
# the ordering should NOT be changed to avoid bugs of ComplexHeatmap package!!
epithelial_cells <- c("Epithelial_tumor", "Epithelial_normal")
lymph_cells <- c("CD8T", "CD4T", "Bcell")
myeloid_cells <- c("Macrophage_HLADRn_CD163n", "Macrophage_CD163p", "Macrophage_HLADRp",
                   "Monocyte", "M_MDSC", "Macrophage_HLADRp_CD163p")
other_cells <- c("B7_H4_cell", "CCR6_cell",  "Ki67_cell", "Unknown",
                 "MMT_HLADRn_CD163n", "MMT_HLADRp", "MMT_HLADRp_CD163p",
                 "DC", "Endothelial")
stromal_cells <- c("CAF_col1p", "Vimentin_cell", "mCAF", "PSC")
grouping_celltype <- data.frame(cell_type10 = c(epithelial_cells, lymph_cells,
                                                myeloid_cells, other_cells, stromal_cells),
                                meta_celltype = c(rep("Epithelial", 2), rep("Lymphocyte", 3),
                                                  rep("Myeloid", 6), rep("Other", 9), rep("Stroma", 4)))
heatmap_dat <- c_mat[grouping_celltype$cell_type10, ]
split_row <- grouping_celltype$meta_celltype

# set colors and block annotation
colors <- jhtools::show_me_the_colors("hyperion_cells")
meta_colors <- c("#B61932", "#7E9A51", "#BEA278", "#BAA4CB", "#568AC2")
labels_meta <- c("Epi", "Lymph", "Myeloid", "Others", "Stroma")
labels_meta2 <- c("Epithelial Cells", "Lymphoid Cells", "Myeloid Cells", "Other Cells", "Stromal Cells")
block_anno = rowAnnotation(Type = anno_block(gp = gpar(fill = meta_colors),
                                             labels = labels_meta,
                                             labels_gp = gpar(col = "white", fontsize = 10)))

col_fun = circlize::colorRamp2(c(-1, 0, 1, 2, 3), c("#5658A5", "#8BCDA3", "#FBF4AA", "#F16943", "#9D1A44"))
prop_data <- cluster_prop[ ,3]
names(prop_data) <- cluster_prop[ ,1]

# bar annotation
meta_all1 <- meta_all %>% mutate(base_excision_eval = case_when(
  nchar(base_excision_eval) == 4 ~ "BRPC_LAPC",
  TRUE ~ base_excision_eval))
coldat <- left_join(coldata, meta_all1[, c("sample_id", "base_excision_eval")])
group_dat <- coldat %>% group_by(cell_type10) %>%
  mutate(allcount = n()) %>%
  group_by(cell_type10, base_excision_eval, allcount) %>%
  summarise(count = n()) %>%
  na.omit() %>% mutate(frac = count / allcount)
group_dat$cell_type_order <- factor(group_dat$cell_type10,
                                    levels = rev(grouping_celltype$cell_type10))
prop_data <- prop_data[match(rownames(heatmap_dat), names(prop_data))]
cell_amount_anno = ComplexHeatmap::rowAnnotation(
  frac = anno_barplot(prop_data,
                      gp = gpar(fill = colors[match(rownames(heatmap_dat), names(colors))])))
text_color <- scales::muted(colors[match(rev(grouping_celltype$cell_type10), names(colors))], l = 50)

# bubble plot
# group_dat <- group_dat %>% mutate(base_excision_eval = case_when(nchar(base_excision_eval) == 4 ~ "BRPC_LAPC", TRUE ~ base_excision_eval))
group_dat$base_order <- factor(group_dat$base_excision_eval,
                               levels = c("RPC", "BRPC_LAPC", "MPC"))
pcirc <- ggplot(group_dat, aes(y = cell_type_order, x = base_order)) +
  geom_point(aes(color = cell_type_order, alpha = frac, size = count)) +
  scale_color_manual(values = colors) +
  scale_size(range = c(1, 15), name = "Fraction of clinical in celltype group") +
  theme(
    axis.text.y = element_text(color = text_color),
    # not officially supported but works anyway, one may also seek for ggtext solution
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  expand_limits(y = c(1, length(levels(group_dat$cell_type_order)) + 0.8))

# main heatmap
ht_main <- Heatmap(heatmap_dat, name = "matrix", heatmap_legend_param = list(title = "exp"),
                   cluster_row_slices = F, col = col_fun, row_km = 0,
                   cluster_rows = F, show_heatmap_legend = F, left_annotation = block_anno,
                   split = split_row, show_row_names = F, row_title = NULL,
                   show_row_dend = F, show_column_dend = F)

ht <- attach_annotation(ht_main, cell_amount_anno, side = "right")
ht <- draw(ht)
grob <- grid.grabExpr(draw(ht))

# legends
heat_legend <- Legend(title = "Intensity", at = -1:3, col_fun = col_fun)
meta_legend <- Legend(labels = labels_meta2, title = "Cell Categories",
                      legend_gp = gpar(fill = meta_colors))
merged_legends <- packLegend(heat_legend, meta_legend)

# merged plot
p <- ggdraw() +
  draw_plot(pcirc, 0, 0.033, 0.27, 0.95) +   ## x, y, width, height
  draw_plot(grob, 0.275, 0.001, 0.52, 0.98) +
  draw_plot(grid.grabExpr(draw(merged_legends)), 0.84, 0.6, 0.02, 0.1)
ggsave(glue::glue("heatmap_main_fig1d.pdf"), p, height = 7, width = 13)


# fig1e updated code from @rujia zheng ------------------------------------
celltype10 <- names(show_me_the_colors("hyperion_cells"))
cld <- colData(sce) %>% as.data.frame()
cld <- cld %>% group_by(stype3) %>% mutate(n_cell = n()) %>%
  group_by(cell_type10, .add = T) %>% mutate(cell_ratio = n()/n_cell)
df <- cld %>% ungroup() %>% dplyr::select(stype3, cell_type10, cell_ratio) %>%
  dplyr::distinct() %>% pivot_wider(names_from = cell_type10, values_from = cell_ratio) %>%
  mutate(across(where(is.numeric), replace_na, 0)) %>%
  pivot_longer(celltype10, names_to = "cell_type", values_to = "ratio")
df$stype3<- factor(df$stype3, levels = c("Surgery_without_chemo", "Surgery_after_chemo", "Punc_pancreas", "Punc_liver", "Paracancerous", "Normal"))
stack <- function (df = df) {
  p <- ggplot(df, mapping = aes(stype3, ratio*100, fill = cell_type)) +
    geom_bar(stat = 'identity', position = position_stack(), width = 0.8) +
    labs(x = '', y = 'Cell percentage') +
    theme_bw() + scale_fill_manual(values = show_me_the_colors("hyperion_cells")) +
    theme(axis.title = element_text(size = 13),
          axis.text.x = element_text(angle = 60, hjust = 1),
          axis.title.x = element_blank(),
          axis.ticks.x =  element_blank(),
          axis.text.y = element_text(color = "black"),
          axis.ticks = element_line(size = 0.25),
          strip.text = element_text(size = 8),
          strip.background = element_rect(fill = NA, color = NA),
          legend.position = "none") +
    scale_y_continuous(expand = c(0,0))
  return(p)
}
p <- stack(df = df)
plot_fn <- "stype2_stackplot_fig1e.pdf"
ggsave(plot_fn, p, width = 2, height = 3, dpi = 300)


# fig1f ------------------------------------





