pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "SingleCellExperiment", "BiocNeighbors",
          "vroom", "jhtools", "glue", "openxlsx", "ggsci", "patchwork", "cowplot", "tidyverse", "dplyr", "survminer", "survival")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "hyperion"
dataset <- "qzhang"
species <- "human"
workdir <- glue::glue("~/projects/{project}/output/{dataset}/human/figure4")
setwd(workdir)

outdir <- glue::glue("~/projects/{project}/output/{dataset}/human/figure4_v1")

#read in data
sce <- readr::read_rds("/cluster/home/jhuang/projects/hyperion/analysis/qzhang/human/steinbock/rds/combined/all_anno.rds")
cl <- readr::read_rds("/cluster/home/yjliu_jh/share/meta_clu3.rds") %>% as_tibble()
sample_chemo_type_list <- readr::read_rds("/cluster/home/yjliu_jh/share/sample_chemo_type_list.rds")
cld <- colData(sce) %>% as_tibble()
sinfo <- metadata(sce)$sampleinfo %>% as_tibble()

#generate cell_type11 & meta_merge
cld <- cld %>% left_join(sinfo %>% dplyr::select(sample_id, base_chemotherapy), by = "sample_id") %>%
  dplyr::mutate(chemo = case_when(stype2 %in% c("puncture_liver", "puncture_pdac","before_chemo") ~ "before_ca",
                                  (stype2 %in% c("tumor", "after_chemo")) & (base_chemotherapy == "no") ~ "before_ca",
                                  TRUE ~ "after_ca")) %>% dplyr::select(-base_chemotherapy) %>%
  dplyr::mutate(cell_type11 = case_when(cell_type10 %in%
                                          c("Macrophage_HLADRp", "Macrophage_HLADRp_CD163p", "Macrophage_CD163p", "Macrophage_HLADRn_CD163n") ~ "Macrophage",
                                        cell_type10 %in% c("MMT_HLADRp_CD163p", "MMT_HLADRp", "MMT_HLADRn_CD163n") ~ "MMT",
                                        TRUE ~ cell_type10)) %>%
  dplyr::mutate(meta_merge = case_when(meta_cluster %notin% c("Tumor_boundary","Stroma_Macro") ~ "Other_metas",
                                       TRUE ~ meta_cluster))

# distance compare --------------------------------------------------------
meta_colors <- c(Tumor_boundary = "#BC3C29FF", Stroma_Macro = "#0072B5FF", Other_metas = "#E18727FF")
datadir <- "/cluster/home/ylxie_jh/projects/hyperion/analysis/qzhang/human/steinbock/downAnalysis/chemo_or_not_DNmac_compare"

# Epithelial_tumor_to_Macrophage_HLADRn_CD163n ----------------------------
df_closecell_distance_SM <- readr::read_rds(glue::glue("{datadir}/SM_community_closecell_dist_k1_Epithelial_tumor_to_Macrophage_HLADRn_CD163n.rds")) %>%
  do.call("rbind", .) %>% as_tibble() %>% dplyr::mutate(meta_merge = "Stroma_Macro")
df_closecell_distance_TB <- readr::read_rds(glue::glue("{datadir}/TB_community_closecell_dist_k1_Epithelial_tumor_to_Macrophage_HLADRn_CD163n.rds")) %>%
  do.call("rbind", .) %>% as_tibble() %>% dplyr::mutate(meta_merge = "Tumor_boundary")
df_closecell_distance_othersnoBT <- readr::read_rds(glue::glue("{datadir}/OthersnoBT_community_closecell_dist_k1_Epithelial_tumor_to_Macrophage_HLADRn_CD163n.rds")) %>%
  do.call("rbind", .) %>% as_tibble() %>% dplyr::mutate(meta_merge = "Other_metas")

df_closecell_distance <- rbind(df_closecell_distance_SM, df_closecell_distance_TB, df_closecell_distance_othersnoBT) %>% as_tibble()
df_closecell_distance <- df_closecell_distance %>% 
  left_join(cld %>% dplyr::rename(from_cell = cell_id), by = c("from_cell", "meta_merge")) %>% 
  dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all)
my_comparisons <- list(c("Stroma_Macro", "Tumor_boundary"), c("Stroma_Macro", "Other_metas"), c("Tumor_boundary", "Other_metas"))
df_closecell_distance$meta_merge <- factor(df_closecell_distance$meta_merge, levels = c("Tumor_boundary", "Stroma_Macro", "Other_metas"))
p1 <- ggboxplot(df_closecell_distance, x = "meta_merge", y = "distance",
               color = "meta_merge",
               outlier.shape = NA,
               size = 0.2, palette = meta_colors) + ylim(0,200) +
  stat_compare_means(aes(label = ..p.signif..),comparisons = my_comparisons, method = "wilcox.test", label.y = c(140, 160, 180), size = 1) +
  theme(axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm')) + rotate_x_text(30)
ggsave("Epithelial_tumor_to_Macrophage_Dn_minDis_meta3_compare_fig4.pdf", p1, width = 1.8, height = 3)
readr::write_rds(p1, "Epithelial_tumor_to_Macrophage_Dn_minDis_meta3_compare_fig4.rds")


# Macrophage_HLADRn_CD163n_to_Epithelial_tumor ----------------------------
df_closecell_distance_SM <- readr::read_rds(glue::glue("{datadir}/SM_community_closecell_dist_k1_Macrophage_HLADRn_CD163n_to_Epithelial_tumor.rds")) %>%
  do.call("rbind", .) %>% as_tibble() %>% dplyr::mutate(meta_merge = "Stroma_Macro")
df_closecell_distance_TB <- readr::read_rds(glue::glue("{datadir}/TB_community_closecell_dist_k1_Macrophage_HLADRn_CD163n_to_Epithelial_tumor.rds")) %>%
  do.call("rbind", .) %>% as_tibble() %>% dplyr::mutate(meta_merge = "Tumor_boundary")
df_closecell_distance_others <- readr::read_rds(glue::glue("{datadir}/OthersnoBT_community_closecell_dist_k1_Macrophage_HLADRn_CD163n_to_Epithelial_tumor.rds")) %>%
  do.call("rbind", .) %>% as_tibble() %>% dplyr::mutate(meta_merge = "Other_metas")

df_closecell_distance <- rbind(df_closecell_distance_SM, df_closecell_distance_TB, df_closecell_distance_others) %>% as_tibble()
df_closecell_distance <- df_closecell_distance %>% 
  left_join(cld %>% dplyr::rename(from_cell = cell_id), by = c("from_cell", "meta_merge")) %>% 
  dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all)
my_comparisons <- list(c("Stroma_Macro", "Tumor_boundary"), c("Stroma_Macro", "Other_metas"), c("Tumor_boundary", "Other_metas"))
df_closecell_distance$meta_merge <- factor(df_closecell_distance$meta_merge, levels = c("Tumor_boundary", "Stroma_Macro", "Other_metas"))

p2 <- ggboxplot(df_closecell_distance, x = "meta_merge", y = "distance",
               color = "meta_merge",
               outlier.shape = NA,
               size = 0.2, palette = meta_colors) + ylim(0,200) +
  stat_compare_means(aes(label = ..p.signif..),comparisons = my_comparisons, method = "wilcox.test", label.y = c(120, 140, 160), size = 1) +
  theme(axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm')) + rotate_x_text(30)
ggsave("Macrophage_Dn_to_Epithelial_tumor_minDis_meta3_compare_fig4.pdf", p2, width = 1.5, height = 3)
readr::write_rds(p2, "Macrophage_Dn_to_Epithelial_tumor_minDis_meta3_compare_fig4.rds")



# os of DN-mac ------------------------------------------------------------
datadir <- glue::glue("~/projects/{project}/analysis/{dataset}/human/steinbock/downAnalysis/cell_type_prop_os")
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}


#summary cell type prop
cell_type_meta3_total <- cld %>% dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>%
  dplyr::filter(meta_cluster != "Bulk_tumor") %>%
  group_by(sample_id, meta_merge) %>% summarise(nt = n()) %>% ungroup()

cell_type_meta3 <- cld %>% dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>%
  dplyr::filter(meta_cluster != "Bulk_tumor") %>%
  group_by(sample_id, cell_type10, meta_merge) %>% summarise(nc = n()) %>% ungroup()

cell_type_meta3$cell_type10 <- factor(cell_type_meta3$cell_type10)
cell_type_meta3 <- cell_type_meta3 %>% tidyr::complete(sample_id, cell_type10, meta_merge, fill = list(nc = 0))
cell_type_meta3 <- left_join(cell_type_meta3_total, cell_type_meta3, by = c("sample_id", "meta_merge")) %>%
  dplyr::mutate(prop = nc/nt) %>% dplyr::select(-c(nt, nc))


#add os info
cell_type_meta3 <- left_join(cell_type_meta3, sinfo, by = "sample_id") %>%
  dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>%
  dplyr::select(sample_id, cell_type10, meta_merge, prop, pfs_state, pfs_month, os_state, os_month)

cell_type_meta3_tmp <- cell_type_meta3 %>% dplyr::filter(cell_type10 == "Macrophage_HLADRn_CD163n")
cell_type_meta3_tmp1 <- cell_type_meta3_tmp %>% dplyr::filter(meta_merge == "Tumor_boundary") %>%
  dplyr::mutate(group_mean = case_when(prop >= mean(prop) ~ "high", prop < mean(prop) ~ "low"))

cell_type_meta3_tmp1_os <- cell_type_meta3_tmp1 %>% dplyr::select(sample_id, os_state, os_month, group_mean) %>% na.omit()

#mean group
psurvx <- ggsurvplot(surv_fit(Surv(os_month, os_state) ~ group_mean, data = cell_type_meta3_tmp1_os), palette = ggsci::pal_nejm("default")(2),
                     legend.labs = levels(droplevels(as.factor(unlist(cell_type_meta3_tmp1_os[, "group_mean"])))),
                     pval=T, risk.table = F, xlim = c(0,75))
ggsave(glue::glue("cell_type_meta3_DN_mac_Tumor_boundary_os_mean_fig4.pdf"),
       plot = psurvx, width = 7, height = 7)
readr::write_rds(psurvx, "cell_type_meta3_DN_mac_Tumor_boundary_os_mean_fig4.rds")


# DN_marker Other_marker top10 pos cell compare ---------------------------------
#total region
sample_chemo_type_list <- readr::read_rds("/cluster/home/yjliu_jh/share/sample_chemo_type_list.rds")
pos_cells10 <- readr::read_rds("/cluster/home/yjliu_jh/share/pos_cell_list_10percent.rds")

Macrophage_Dn_total <- cld %>% dplyr::filter(cell_type10 == "Macrophage_HLADRn_CD163n") %>%
  dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>%
  group_by(sample_id) %>% summarise(nt = n()) %>% ungroup()

Macrophage_Other_total <- cld %>% dplyr::filter(cell_type10 %in% c("Macrophage_HLADRp",
                                                                   "Macrophage_HLADRp_CD163p",
                                                                   "Macrophage_CD163p")) %>%
  dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>%
  group_by(sample_id) %>% summarise(nt = n()) %>% ungroup()

macros <- c("Ki67_pos_Macrophage", "EGFR_pos_Macrophage", "CCR6_pos_Macrophage")
plist <- list()
#top10 prop
for (i in 1:length(macros)) {
  Macrophage_Dn_pos <- cld %>%
    dplyr::filter(cell_type10 == "Macrophage_HLADRn_CD163n") %>%
    dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>%
    dplyr::filter(cell_id %in% pos_cells10[[macros[i]]]) %>%
    group_by(sample_id) %>% summarise(np = n()) %>% ungroup()
  Macrophage_Dn_pos <- left_join(Macrophage_Dn_total, Macrophage_Dn_pos,
                                 by = c("sample_id")) %>%
    replace_na(list(np = 0)) %>% dplyr::mutate(prop = np/nt) %>%
    ungroup() %>% dplyr::select(-c(nt, np)) %>%
    dplyr::mutate(group = "DN-Mp")
  
  Macrophage_Other_pos <- cld %>%
    dplyr::filter(cell_type10 %in% c("Macrophage_HLADRp",
                                     "Macrophage_HLADRp_CD163p",
                                     "Macrophage_CD163p")) %>%
    dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>%
    dplyr::filter(cell_id %in% pos_cells10[[macros[i]]]) %>%
    group_by(sample_id) %>% summarise(np = n()) %>% ungroup()
  Macrophage_Other_pos <- left_join(Macrophage_Other_total, Macrophage_Other_pos,
                                    by = c("sample_id")) %>%
    replace_na(list(np = 0)) %>% dplyr::mutate(prop = np/nt) %>%
    ungroup() %>% dplyr::select(-c(nt, np)) %>%
    dplyr::mutate(group = "Other-Mp")
  
  Macrophage_top10_pos_com <- rbind(Macrophage_Dn_pos, Macrophage_Other_pos)
  
  plist[[i]] <- ggboxplot(Macrophage_top10_pos_com, x = "group", y = "prop",
                          color = "group", palette = pal_nejm("default")(2),
                          add = "jitter", size = 0.2, add.params = list(size = 0.2)) +
    ylab(macros[i]) +
    stat_compare_means(aes(group = group),label = "p.signif", label.y = 0.6) +
    theme(axis.title.y = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.title = element_blank(),
          legend.key.size = unit(0.1, 'cm'))
}
# fil_name <- "Macrophage_DN_compare_with_Others_prop_perSample_totalRegion_fig4.pdf"
# multi_plot(fig_fn = fil_name, width = 8, height = 7,
#           p_list = plist, ncol = 5, nrow = 3, tag_levels = NULL)
readr::write_rds(plist, "Ki67_EGFR_CCR6_DN_mac_compare_with_Others_prop_perSample_totalRegion_compare_fig4.rds")



# PD_L1/EGFR pos DN-mac in 3 metas ----------------------------------------
pos_cells10 <- readr::read_rds("/cluster/home/yjliu_jh/share/pos_cell_list_10percent.rds")
meta_colors <- c(Tumor_boundary = "#BC3C29FF", Stroma_Macro = "#0072B5FF", Other_metas = "#E18727FF")

Macrophage_Dn_total <- cld %>% dplyr::filter(cell_type10 == "Macrophage_HLADRn_CD163n") %>%
  dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>%
  group_by(sample_id, meta_merge) %>% summarise(nt = n()) %>% ungroup()

#PD-L1
Macrophage_Dn_pos <- cld %>% dplyr::filter(cell_type10 == "Macrophage_HLADRn_CD163n") %>%
  dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>%
  dplyr::filter(cell_id %in% pos_cells10[["PD_L1_pos_Macrophage"]]) %>%
  group_by(sample_id, meta_merge) %>% summarise(np = n()) %>% ungroup()
Macrophage_Dn_pos <- left_join(Macrophage_Dn_total, Macrophage_Dn_pos, by = c("sample_id", "meta_merge")) %>%
  replace_na(list(np = 0)) %>% dplyr::mutate(prop = np/nt) %>%
  ungroup() %>% dplyr::select(-c(nt, np))
Macrophage_Dn_pos$meta_merge <- factor(Macrophage_Dn_pos$meta_merge, levels = c("Tumor_boundary", "Stroma_Macro", "Other_metas"))
my_comparisons <- list(c("Stroma_Macro", "Tumor_boundary"), c("Stroma_Macro", "Other_metas"), c("Tumor_boundary", "Other_metas"))
p <- ggboxplot(Macrophage_Dn_pos, x = "meta_merge", y = "prop",
               color = "meta_merge",
               outlier.shape = NA,
               size = .2,
               palette = meta_colors) + ylab(glue::glue("Proportion of PD-L1+ DN-mp")) + ylim(0, .5) +
  stat_compare_means(aes(label = ..p.signif..),comparisons = my_comparisons, method = "wilcox.test", size = 1, label.y = c(0.5,0.55,0.6,0.65)) +
  stat_summary(fun.data = function(x) data.frame(y = 0.2, label = paste("mean=", round(mean(x),2))), geom="text", size = 1) +
  theme(axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm')) + rotate_x_text(30)
ggsave(glue::glue("fig4_PD_L1_pos_Macrophage_Dn_sample_meta3_compare.pdf"), p, width = 2, height = 2)
readr::write_rds(p, "fig4_PD_L1_pos_Macrophage_Dn_sample_meta3_compare.rds")


#EGFR
Macrophage_Dn_pos <- cld %>% dplyr::filter(cell_type10 == "Macrophage_HLADRn_CD163n") %>%
  dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>%
  dplyr::filter(cell_id %in% pos_cells10[["EGFR_pos_Macrophage"]]) %>%
  group_by(sample_id, meta_merge) %>% summarise(np = n()) %>% ungroup()
Macrophage_Dn_pos <- left_join(Macrophage_Dn_total, Macrophage_Dn_pos, by = c("sample_id", "meta_merge")) %>%
  replace_na(list(np = 0)) %>% dplyr::mutate(prop = np/nt) %>%
  ungroup() %>% dplyr::select(-c(nt, np))
Macrophage_Dn_pos$meta_merge <- factor(Macrophage_Dn_pos$meta_merge, levels = c("Tumor_boundary", "Stroma_Macro", "Other_metas"))
my_comparisons <- list(c("Stroma_Macro", "Tumor_boundary"), c("Stroma_Macro", "Other_metas"), c("Tumor_boundary", "Other_metas"))
p <- ggboxplot(Macrophage_Dn_pos, x = "meta_merge", y = "prop",
               color = "meta_merge",
               outlier.shape = NA,
               size = .2,
               palette = meta_colors) + ylab(glue::glue("Proportion of EGFR+ DN-mp")) + ylim(0, .5) +
  stat_compare_means(aes(label = ..p.signif..),comparisons = my_comparisons, method = "wilcox.test", size = 1, label.y = c(0.5,0.55,0.6,0.65)) +
  stat_summary(fun.data = function(x) data.frame(y = 0.2, label = paste("mean=", round(mean(x),2))), geom="text", size = 1) +
  theme(axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm')) + rotate_x_text(30)
ggsave(glue::glue("fig4_EGFR_pos_Macrophage_Dn_sample_meta3_compare.pdf"), p, width = 2, height = 2)
readr::write_rds(p, "fig4_EGFR_pos_Macrophage_Dn_sample_meta3_compare.rds")



p1 <- readr::read_rds("cell_type_meta3_DN_mac_Tumor_boundary_os_mean_fig4.rds")
p2 <- readr::read_rds("Epithelial_tumor_to_Macrophage_Dn_minDis_meta3_compare_fig4.rds")
p3 <- readr::read_rds("Macrophage_Dn_to_Epithelial_tumor_minDis_meta3_compare_fig4.rds")
p1$plot <- p1$plot + theme(axis.title.y = element_text(size = 6),
                           axis.text.y = element_text(size = 6),
                           axis.title.x = element_blank(),
                           axis.text.x = element_text(size = 6),
                           legend.text = element_text(size = 6),
                           legend.title = element_text(size = 6),
                           legend.key.size = unit(0.1, 'cm'))
p <- p1$plot + p2 + p3 + plot_layout(guides = "collect") 
ggsave(glue::glue("{outdir}/fig4_panel1.pdf"), p, height = 2, width = 9)


p4 <- readr::read_rds("fig4_PD_L1_pos_Macrophage_Dn_sample_meta3_compare.rds")
p5 <- readr::read_rds("fig4_EGFR_pos_Macrophage_Dn_sample_meta3_compare.rds")
p6 <- readr::read_rds("Ki67_EGFR_CCR6_DN_mac_compare_with_Others_prop_perSample_totalRegion_compare_fig4.rds")

p <- p4|p5|p6[[1]]|p6[[2]]|plot_layout(guides = "collect") 
ggsave(glue::glue("{outdir}/fig4_panel2.pdf"), p, height = 2.5, width = 11.5)



