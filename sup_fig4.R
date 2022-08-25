pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "SingleCellExperiment", "BiocNeighbors",
          "vroom", "jhtools", "glue", "openxlsx", "ggsci", "patchwork", "cowplot", "tidyverse", "dplyr", "survminer", "survival")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "hyperion"
dataset <- "qzhang"
species <- "human"
workdir <- glue::glue("~/projects/{project}/output/{dataset}/{species}/sfigure4")
setwd(workdir)
outdir <- glue::glue("~/projects/{project}/output/{dataset}/{species}/sfigure4_v1")

# os of DN-mac in SM and Others------------------------------------------------------------
datadir <- glue::glue("~/projects/{project}/analysis/{dataset}/{species}/steinbock/downAnalysis/cell_type_prop_os")
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

#read in data
sce <- readr::read_rds("/cluster/home/jhuang/projects/hyperion/analysis/qzhang/human/steinbock/rds/combined/all_anno.rds")
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

#Stroma_Macro
cell_type_meta3_tmp1 <- cell_type_meta3_tmp %>% dplyr::filter(meta_merge == "Stroma_Macro") %>%
  dplyr::mutate(group_mean = case_when(prop >= mean(prop) ~ "high", prop < mean(prop) ~ "low"))

cell_type_meta3_tmp1_os <- cell_type_meta3_tmp1 %>% dplyr::select(sample_id, os_state, os_month, group_mean) %>% na.omit()

#mean group
psurvx <- ggsurvplot(surv_fit(Surv(os_month, os_state) ~ group_mean, data = cell_type_meta3_tmp1_os), palette = ggsci::pal_nejm("default")(2),
                     legend.labs = levels(droplevels(as.factor(unlist(cell_type_meta3_tmp1_os[, "group_mean"])))),
                     pval=T, risk.table = F, xlim = c(0,75), censor.size=0.3, size = 0.3, font.tickslab = c(6))
ggsave(glue::glue("cell_type_meta3_DN_mac_Stroma_Macro_os_mean_sfig4.pdf"),
       plot = psurvx, width = 7, height = 7)
readr::write_rds(psurvx, "cell_type_meta3_DN_mac_Stroma_Macro_os_mean_sfig4.rds")


#Other_metas
cell_type_meta3_tmp1 <- cell_type_meta3_tmp %>% dplyr::filter(meta_merge == "Other_metas") %>%
  dplyr::mutate(group_mean = case_when(prop >= mean(prop) ~ "high", prop < mean(prop) ~ "low"))

cell_type_meta3_tmp1_os <- cell_type_meta3_tmp1 %>% dplyr::select(sample_id, os_state, os_month, group_mean) %>% na.omit()

#mean group
psurvx <- ggsurvplot(surv_fit(Surv(os_month, os_state) ~ group_mean, data = cell_type_meta3_tmp1_os), palette = ggsci::pal_nejm("default")(2),
                     legend.labs = levels(droplevels(as.factor(unlist(cell_type_meta3_tmp1_os[, "group_mean"])))),
                     pval=T, risk.table = F, xlim = c(0,75), censor.size=0.3, size = 0.3, font.tickslab = c(6))
ggsave(glue::glue("cell_type_meta3_DN_mac_Other_metas_os_mean_sfig4.pdf"),
       plot = psurvx, width = 7, height = 7)
readr::write_rds(psurvx, "cell_type_meta3_DN_mac_Other_metas_os_mean_sfig4.rds")


###########################################################################
# PD-L1/EGFR intensity of DN-mac in different distance from DN-mac to tumor  ----------------------------------
#Macrophage_HLADRn_CD163n to Epithelial_tumor
target_cell <- "Epithelial_tumor"
query_cell <- "Macrophage_HLADRn_CD163n"
datadir <- glue::glue("~/projects/{project}/analysis/{dataset}/{species}/steinbock/downAnalysis/macro_Dn_epi_dis")

list_closecell_distance <- readr::read_rds(glue::glue("{datadir}/TB_list_closecell_dist_k{k}_{query_cell}_to_{target_cell}.rds"))

df_closecell_distance <- do.call("rbind", list_closecell_distance)
df_closecell_distance <- left_join(df_closecell_distance, cld %>% dplyr::rename("from_cell" = "cell_id"), by = "from_cell")
df_closecell_distance <- df_closecell_distance %>% dplyr::mutate(group = case_when(distance < quantile(distance, probs = 0.25) ~ "g1",
                                                                                   distance >= quantile(distance, probs = 0.25) &
                                                                                     distance < quantile(distance, probs = 0.5) ~ "g2",
                                                                                   distance >= quantile(distance, probs = 0.5) &
                                                                                     distance < quantile(distance, probs = 0.75) ~ "g3",
                                                                                   distance >= quantile(distance, probs = 0.75) ~ "g4"))

my_comparisons <- list(c("g1", "g2"), c("g1", "g3"), c("g1", "g4"), c("g2", "g3"), c("g2", "g4"), c("g3", "g4"))

df_closecell_distance$group <- factor(df_closecell_distance$group, levels = c("g1", "g2", "g3", "g4"))
p <- ggboxplot(df_closecell_distance, x = "group", y = "PD_L1",
               color = "group",
               outlier.shape = NA,
               size = .2,
               palette = ggsci::pal_nejm("default")(4)) + ylab(glue::glue("PD_L1_intensity_of_DN_mp")) +
  stat_compare_means(aes(label = ..p.signif..),comparisons = my_comparisons, method = "wilcox.test", size = 1, label.y = c(0.5,0.55,0.6,0.65)) +
  stat_summary(fun.data = function(x) data.frame(y = 0.2, label = paste("mean=", round(mean(x),2))), geom="text", size = 1) +
  theme(axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm'))
ggsave(glue::glue("sfig4_PD_L1_intensity_of_DNmac_to_tumor.pdf"), p, width = 2, height = 2)
readr::write_rds(p, "sfig4_PD_L1_intensity_of_DNmac_to_tumor.rds")

p <- ggboxplot(df_closecell_distance, x = "group", y = "EGFR",
               color = "group",
               outlier.shape = NA,
               size = .2,
               palette = ggsci::pal_nejm("default")(4)) + ylab(glue::glue("EGFR_intensity_of_DN_mp")) +
  stat_compare_means(aes(label = ..p.signif..),comparisons = my_comparisons, method = "wilcox.test", size = 1, label.y = c(0.5,0.55,0.6,0.65)) +
  stat_summary(fun.data = function(x) data.frame(y = 0.2, label = paste("mean=", round(mean(x),2))), geom="text", size = 1) +
  theme(axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm'))
ggsave(glue::glue("sfig4_EGFR_intensity_of_DNmac_to_tumor.pdf"), p, width = 2, height = 2)
readr::write_rds(p, "sfig4_EGFR_intensity_of_DNmac_to_tumor.rds")


###########################################################################
# DN-mac and other-mac top10 pos compare ----------------------------------
sample_chemo_type_list <- readr::read_rds("/cluster/home/yjliu_jh/share/sample_chemo_type_list.rds")
pos_cells10 <- readr::read_rds("/cluster/home/yjliu_jh/share/pos_cell_list_10percent.rds")
cld <- colData(sce) %>% as_tibble()
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

#total region
Macrophage_Dn_total <- cld %>% dplyr::filter(cell_type10 == "Macrophage_HLADRn_CD163n") %>%
  dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>%
  group_by(sample_id) %>% summarise(nt = n()) %>% ungroup()

Macrophage_Other_total <- cld %>% dplyr::filter(cell_type10 %in% c("Macrophage_HLADRp",
                                                                   "Macrophage_HLADRp_CD163p",
                                                                   "Macrophage_CD163p")) %>%
  dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>%
  group_by(sample_id) %>% summarise(nt = n()) %>% ungroup()


macros <- c("Arginase_1_pos_Macrophage", "PD1_pos_Macrophage", "Caspase3_pos_Macrophage",
            "PD_L1_pos_Macrophage", "LAG_3_pos_Macrophage", "B7_H4_pos_Macrophage",
            "Vista_pos_Macrophage", "FoxP3_pos_Macrophage", "CD25_pos_Macrophage")
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
# fil_name <- "Macrophage_DN_compare_with_Others_prop_perSample_totalRegion_sfig4.pdf"
# multi_plot(fig_fn = fil_name, width = 8, height = 7,
#           p_list = plist, ncol = 5, nrow = 3, tag_levels = NULL)
readr::write_rds(plist, "sfig4_top10pos_DN_mac_compare_with_Others_prop_perSample_totalRegion_compare.rds")

###########################################################################
# function define ---------------------------------------------------------
#correlation plot function
posDNmac_posOther_correlation_func <- function(cell_pairs1, cell_pairs2){
  cld_nochemo_pos1_noBT <- cld %>% dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>%
    dplyr::filter(meta_cluster != "Bulk_tumor") %>%
    dplyr::filter(cell_type10 == "Macrophage_HLADRn_CD163n") %>%
    dplyr::filter(cell_id %in% pos_cells10[[cell_pairs1]]) %>%
    group_by(meta_merge, sample_id) %>% summarise(np = n())
  cld_nochemo_total1_noBT <- cld %>% dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>%
    dplyr::filter(meta_cluster != "Bulk_tumor") %>%
    dplyr::filter(cell_type10 == "Macrophage_HLADRn_CD163n") %>%
    group_by(meta_merge, sample_id) %>% summarise(nt = n())
  cld_nochemo_pos1_noBT <- left_join(cld_nochemo_total1_noBT, cld_nochemo_pos1_noBT, by = c("meta_merge", "sample_id")) %>%
    replace_na(list(np = 0)) %>% dplyr::mutate(prop = np/nt) %>% dplyr::select(-c(nt, np)) %>% ungroup()
  
  cld_nochemo_pos2_noBT <- cld %>% dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>%
    dplyr::filter(meta_cluster != "Bulk_tumor") %>%
    dplyr::filter(cell_id %in% pos_cells10[[cell_pairs2]]) %>%
    group_by(meta_merge, sample_id) %>% summarise(np = n())
  cld_nochemo_total2_noBT <- cld %>% dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>%
    dplyr::filter(meta_cluster != "Bulk_tumor") %>%
    dplyr::filter(cell_type11 %in% str_replace_all(cell_pairs2, "^.*pos_", "")) %>%
    group_by(meta_merge, sample_id) %>% summarise(nt = n())
  cld_nochemo_pos2_noBT <- left_join(cld_nochemo_total2_noBT, cld_nochemo_pos2_noBT, by = c("meta_merge", "sample_id")) %>%
    replace_na(list(np = 0)) %>% dplyr::mutate(prop = np/nt) %>% dplyr::select(-c(nt, np)) %>% ungroup()
  
  cld_nochemo_com <- inner_join(cld_nochemo_pos1_noBT, cld_nochemo_pos2_noBT, by = c("meta_merge", "sample_id"))
  cld_nochemo_com$meta_merge <- factor(cld_nochemo_com$meta_merge, levels = c(c("Tumor_boundary","Stroma_Macro","Other_metas")))
  
  #fig5 part3 main figure
  p <- ggscatter(cld_nochemo_com, x = "prop.x", y = "prop.y",
                 color = "meta_merge", palette = meta_colors, size = .4,
                 add = "reg.line", conf.int = TRUE, star.plot.lwd = .1,
                 cor.coef.size = .01) + xlim(0,1) +
    ylim(0,1) + xlab(glue::glue("{cell_pairs1} HLADR- CD163-")) + ylab(cell_pairs2) +
    stat_cor(aes(color = meta_merge), label.x = 0.2) +
    theme(axis.title.y = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.title.x = element_text(size = 6),
          axis.text.x = element_text(size = 6),
          legend.text = element_text(size = 2),
          legend.title = element_text(size = 2),
          legend.key.size = unit(0.1, 'cm'))
  return(p)
}

dist_to_DNmac_compare_meta3_func <- function(datadir, k, query_cell, target_cell, my_comparisons, meta_colors){
  #Tumor_boundary
  list_closecell_distance <- read_rds(glue::glue("{datadir}/TB_list_closecell_dist_k{k}_{query_cell}_to_{target_cell}.rds"))
  df_closecell_distance <- do.call("rbind", list_closecell_distance)
  TB_dist <- df_closecell_distance %>%
    dplyr::filter(to_cell %in% (cld %>% dplyr::filter(cell_type10 == "Macrophage_HLADRn_CD163n") %>% .$cell_id)) %>%
    dplyr::select(from_cell, distance) %>%
    dplyr::mutate(group = "Tumor_boundary") %>% as_tibble()
  
  #Stroma_Macro
  list_closecell_distance <- read_rds(glue::glue("{datadir}/SM_list_closecell_dist_k{k}_{query_cell}_to_{target_cell}.rds"))
  df_closecell_distance <- do.call("rbind", list_closecell_distance)
  SM_dist <- df_closecell_distance %>%
    dplyr::filter(to_cell %in% (cld %>% dplyr::filter(cell_type10 == "Macrophage_HLADRn_CD163n") %>% .$cell_id)) %>%
    dplyr::select(from_cell, distance) %>%
    dplyr::mutate(group = "Stroma_Macro") %>% as_tibble()
  
  #others no Bulk
  list_closecell_distance <- read_rds(glue::glue("{datadir}/OthersnoBT_list_closecell_dist_k{k}_{query_cell}_to_{target_cell}.rds"))
  df_closecell_distance <- do.call("rbind", list_closecell_distance)
  Other_dist <- df_closecell_distance %>%
    dplyr::filter(to_cell %in% (cld %>% dplyr::filter(cell_type10 == "Macrophage_HLADRn_CD163n") %>% .$cell_id)) %>%
    dplyr::select(from_cell, distance) %>%
    dplyr::mutate(group = "Other_metas") %>% as_tibble()
  
  df_dist <- rbind(TB_dist, SM_dist, Other_dist)
  df_dist <- df_dist %>%
    dplyr::filter(from_cell %in% (cld %>% dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>% .$cell_id))
  
  plist <- list()
  #density
  plist[[1]] <- ggdensity(df_dist, x = "distance",
                 add = "mean",
                 color = "group",
                 palette = meta_colors) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.key.size = unit(0.1, 'cm'))
  #boxplot
  plist[[2]] <- ggboxplot(df_dist, x = "group", y = "distance",
                 color = "group",
                 outlier.shape = NA,
                 size = .2,
                 palette = meta_colors) + rotate_x_text(30) +
    stat_compare_means(aes(label = ..p.signif..),comparisons = my_comparisons, method = "wilcox.test", label.y =c(800,900,1000), size = 1) +
    stat_summary(fun.data = function(x) data.frame(y = 800, label = paste("mean=", round(mean(x),2))), geom="text", size = 1) +
    theme(axis.title.y = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.title.x = element_text(size = 6),
          axis.text.x = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.key.size = unit(0.1, 'cm'))
  return(plist)
}

DNmac_to_Other_intensity_func <- function(datadir, k, query_cell, target_cell, marker){
  #total region
  list_closecell_distance <- readr::read_rds(glue::glue("{datadir}/list_closecell_dist_k{k}_{query_cell}_to_{target_cell}.rds"))
  df_closecell_distance <- do.call("rbind", list_closecell_distance)
  df_closecell_distance <- left_join(df_closecell_distance, cld %>% dplyr::rename("to_cell" = "cell_id"), by = "to_cell") %>%
    dplyr::filter(from_cell %in% (cld %>% dplyr::filter(cell_type10 == "Macrophage_HLADRn_CD163n") %>% .$cell_id)) %>%
    dplyr::filter(from_cell %in% (cld %>% dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>% .$cell_id))
  df_closecell_distance <- df_closecell_distance %>% dplyr::mutate(group = case_when(distance < quantile(distance, probs = 0.25) ~ "g1",
                                                                                     distance >= quantile(distance, probs = 0.25) &
                                                                                       distance < quantile(distance, probs = 0.5) ~ "g2",
                                                                                     distance >= quantile(distance, probs = 0.5) &
                                                                                       distance < quantile(distance, probs = 0.75) ~ "g3",
                                                                                     distance >= quantile(distance, probs = 0.75) ~ "g4"))
  
  my_comparisons <- list(c("g1", "g2"), c("g1", "g3"), c("g1", "g4"), c("g2", "g3"), c("g2", "g4"), c("g3", "g4"))
  df_closecell_distance$group <- factor(df_closecell_distance$group, levels = c("g1", "g2", "g3", "g4"))
  
  p <- ggboxplot(df_closecell_distance, x = "group", y = marker,
                 color = "group",
                 outlier.shape = NA,
                 size = .2,
                 palette = ggsci::pal_nejm("default")(4)) + ylab(glue::glue("{marker}_intensity_of_{target_cell}")) +
    stat_compare_means(aes(label = ..p.signif..),comparisons = my_comparisons, method = "wilcox.test", size = 1, label.y = c(0.5,0.55,0.6,0.65)) +
    stat_summary(fun.data = function(x) data.frame(y = 0.2, label = paste("mean=", round(mean(x),2))), geom="text", size = 1) +
    theme(axis.title.y = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.key.size = unit(0.1, 'cm'))
  return(p)
}


###########################################################################
# ki67 DN-mac and Cas3+ CD8T ----------------------------------------------
# part3 Ki67 DN-mac Cas3CD8T correlation ----------------------------------
cld <- colData(sce) %>% as_tibble()
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

meta_colors <- c(Tumor_boundary = "#BC3C29FF", Stroma_Macro = "#0072B5FF", Other_metas = "#E18727FF")

cell_pairs1 <- "Ki67_pos_Macrophage"
cell_pairs2 <- "Caspase3_pos_CD8T"
p1 <- posDNmac_posOther_correlation_func(cell_pairs1, cell_pairs2)
ggsave("sfig4_DN_mac_Ki67_cor_with_Cas3CD8T.pdf", p1, width = 4, height = 4)
readr::write_rds(p1, "sfig4_DN_mac_Ki67_cor_with_Cas3CD8T.rds")

# distance of Cas+CD8T to Ki67+ DN-mac in 3 meta --------------------------
datadir <- glue::glue("~/projects/{project}/analysis/{dataset}/human/steinbock/downAnalysis/meta3_cor_cells_distance")
meta_colors <- c(Tumor_boundary = "#BC3C29FF", Stroma_Macro = "#0072B5FF", Other_metas = "#E18727FF")
k <- 1
query_cell <- "Caspase3_pos_CD8T"
target_cell <- "Ki67_pos_Macrophage"
my_comparisons <- list(c("Stroma_Macro", "Tumor_boundary"), c("Stroma_Macro", "Other_metas"), c("Tumor_boundary", "Other_metas"))
p2 <- dist_to_DNmac_compare_meta3_func(datadir, k, query_cell, target_cell, my_comparisons, meta_colors)
ggsave("sfig4_DN_mac_Ki67_dist_with_Cas3CD8T_density.pdf", p2[[1]], width = 4, height = 4)
ggsave("sfig4_DN_mac_Ki67_dist_with_Cas3CD8T_boxplot.pdf", p2[[2]], width = 4, height = 4)
readr::write_rds(p2, glue::glue("sfig4_meta3_dist_{query_cell}_to_{target_cell}.rds"))

# marker intensity of CD8T in distance groups of ki67 DN-mac to Cas3CD8T--------
datadir <- glue::glue("~/projects/{project}/analysis/{dataset}/{species}/steinbock/downAnalysis/meta3_cor_cells_dist_marker")
cld <- colData(sce) %>% as_tibble()
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
comcc <- assay(sce, "compcounts_censored")
comcc_t <- as.data.frame(t(comcc))
comcc_t$cell_id <- rownames(comcc_t)
cld <- left_join(cld, comcc_t, by = "cell_id")

#plot
k <- 1
query_cell <- "Ki67_pos_Macrophage"
target_cell <- "CD8T"
marker <- "Caspase3"
p3 <- DNmac_to_Other_intensity_func(datadir, k, query_cell, target_cell, marker)
ggsave("sfig4_DN_mac_Ki67_Cas3CD8T_intensity.pdf", p3, width = 4, height = 4)
readr::write_rds(p3, glue::glue("sfig4_Total_{query_cell}_Dn_to_{target_cell}_distance_{marker}_intensity.rds"))

#plot
k <- 1
query_cell <- "Ki67_pos_Macrophage"
target_cell <- "CD8T"
marker <- "Granzyme_B"
p4 <- DNmac_to_Other_intensity_func(datadir, k, query_cell, target_cell, marker)
ggsave("sfig4_DN_mac_Ki67_GB_CD8T_intensity.pdf", p4, width = 4, height = 4)
readr::write_rds(p4, glue::glue("sfig4_Total_{query_cell}_Dn_to_GB_{target_cell}_distance_{marker}_intensity.rds"))

###########################################################################
# PDL1 tumor and Ki67 DN-mac ----------------------------------------------
# Ki67 DN-mac PD-L1 Tumor correlation ----------------------------------
cell_pairs1 <- "Ki67_pos_Macrophage"
cell_pairs2 <- "PD_L1_pos_Epithelial_tumor"

p5 <- posDNmac_posOther_correlation_func(cell_pairs1, cell_pairs2)
ggsave("sfig4_DN_mac_Ki67_cor_with_PDL1Tumor.pdf", p5, width = 4, height = 4)
readr::write_rds(p5, "sfig4_DN_mac_Ki67_cor_with_PDL1Tumor.rds")

# distance of PDL1 Tumor to Ki67+ DN-mac in 3 meta --------------------------
datadir <- glue::glue("~/projects/{project}/analysis/{dataset}/human/steinbock/downAnalysis/meta3_cor_cells_distance")
meta_colors <- c(Tumor_boundary = "#BC3C29FF", Stroma_Macro = "#0072B5FF", Other_metas = "#E18727FF")

k <- 1
query_cell <- "PD_L1_pos_Epithelial_tumor"
target_cell <- "Ki67_pos_Macrophage"
my_comparisons <- list(c("Stroma_Macro", "Tumor_boundary"), c("Stroma_Macro", "Other_metas"), c("Tumor_boundary", "Other_metas"))
p6 <- dist_to_DNmac_compare_meta3_func(datadir, k, query_cell, target_cell, my_comparisons, meta_colors)
ggsave("sfig4_DN_mac_Ki67_dist_with_PDL1tumor_density.pdf", p6[[1]], width = 4, height = 4)
ggsave("sfig4_DN_mac_Ki67_dist_with_PDL1tumor_boxplot.pdf", p6[[2]], width = 4, height = 4)
readr::write_rds(p6, glue::glue("sfig4_meta3_dist_{query_cell}_to_{target_cell}.rds"))

# marker intensity of Epithelial_tumor in distance groups of ki67 DN-mac to Epithelial_tumor--------
datadir <- glue::glue("~/projects/{project}/analysis/{dataset}/{species}/steinbock/downAnalysis/meta3_cor_cells_dist_marker")

#plot
k <- 1
query_cell <- "Ki67_pos_Macrophage"
target_cell <- "Epithelial_tumor"
marker <- "PD_L1"
p7 <- DNmac_to_Other_intensity_func(datadir, k, query_cell, target_cell, marker)
ggsave("sfig4_DN_mac_Ki67_PDL1tumor_intensity.pdf", p7, width = 4, height = 4)
readr::write_rds(p7, glue::glue("sfig4_Total_{query_cell}_Dn_to_{target_cell}_distance_{marker}_intensity.rds"))


# os of PD-L1 tumor in total region ---------------------------------
epi_total <- cld %>%
  dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>%
  dplyr::filter(cell_type10 == "Epithelial_tumor") %>%
  group_by(sample_id) %>% summarise(nt = n())

epi_pos <- cld %>%
  dplyr::filter(sample_id %in% sample_chemo_type_list$no_chemo_all) %>%
  dplyr::filter(cell_id %in% pos_cells10[["PD_L1_pos_Epithelial_tumor"]]) %>%
  group_by(sample_id) %>% summarise(nc = n()) %>%
  left_join(epi_total, by = "sample_id") %>%
  dplyr::mutate(prop = nc/nt) %>%
  left_join(sinfo, by = "sample_id") %>%
  dplyr::select(sample_id, prop, os_state, os_month) %>%
  dplyr::mutate(group_43 = case_when(prop >= quantile(prop, probs = 0.75) ~ "high", prop < quantile(prop, probs = 0.75) ~ "low"))
epi_pos_os <- epi_pos %>% dplyr::select(sample_id, os_state, os_month, group_43) %>% na.omit()

p8 <- ggsurvplot(surv_fit(Surv(os_month, os_state) ~ group_43, data = epi_pos_os), palette = ggsci::pal_nejm("default")(2),
                  legend.labs = levels(droplevels(as.factor(unlist(epi_pos_os[, "group_43"])))), pval=T, risk.table = F,
                  censor.size=0.3, size = 0.3, font.tickslab = c(6))
ggsave(glue::glue("sfig4_PD_L1_pos_Epithelial_tumor_prop_os_group43.pdf"),
       plot = p8, width = 3, height = 3)
readr::write_rds(p8, glue::glue("sfig4_PD_L1_pos_Epithelial_tumor_prop_os_group43.rds"))


# plot --------------------------------------------------------------------
p1 <- read_rds("cell_type_meta3_DN_mac_Stroma_Macro_os_mean_sfig4.rds")
p2 <- read_rds("cell_type_meta3_DN_mac_Other_metas_os_mean_sfig4.rds")
plist <- read_rds("sfig4_top10pos_DN_mac_compare_with_Others_prop_perSample_totalRegion_compare.rds")
p3 <- read_rds("sfig4_DN_mac_Ki67_cor_with_Cas3CD8T.rds")
p4 <- read_rds("sfig4_meta3_dist_Caspase3_pos_CD8T_to_Ki67_pos_Macrophage.rds")
p5 <- read_rds("sfig4_Total_Ki67_pos_Macrophage_Dn_to_CD8T_distance_Caspase3_intensity.rds")
p6 <- read_rds("sfig4_DN_mac_Ki67_cor_with_PDL1Tumor.rds")
p7 <- read_rds("sfig4_meta3_dist_PD_L1_pos_Epithelial_tumor_to_Ki67_pos_Macrophage.rds")
p8 <- read_rds("sfig4_Total_Ki67_pos_Macrophage_Dn_to_Epithelial_tumor_distance_PD_L1_intensity.rds")
p9 <- read_rds("sfig4_PD_L1_pos_Epithelial_tumor_prop_os_group43.rds")
p10 <- read_rds("sfig4_Total_Ki67_pos_Macrophage_Dn_to_GB_CD8T_distance_Granzyme_B_intensity.rds")
p11 <- read_rds("sfig4_PD_L1_intensity_of_DNmac_to_tumor.rds")
p12 <- read_rds("sfig4_EGFR_intensity_of_DNmac_to_tumor.rds")

p1$plot <- p1$plot + xlab("Months") +
  theme(axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 2),
        legend.title = element_text(size = 2),
        legend.key.size = unit(0.1, 'cm'))

p2$plot <- p2$plot + xlab("Months") +
  theme(axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 2),
        legend.title = element_text(size = 2),
        legend.key.size = unit(0.1, 'cm'))

p9$plot <- p16$plot + xlab("Months") +
  theme(axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 2),
        legend.title = element_text(size = 2),
        legend.key.size = unit(0.1, 'cm'))


p12$plot <- p12$plot + xlab("Months") +
  theme(axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 2),
        legend.title = element_text(size = 2),
        legend.key.size = unit(0.1, 'cm'))


pn <- (p1$plot|p2$plot|p11|p12)/(plist[[1]]|plist[[2]]|plist[[3]])/(plist[[4]]|plist[[5]]|plist[[6]])/(plist[[7]]|plist[[8]]|plist[[9]])
ggsave(glue::glue("{outdir}/sfig4_p1.pdf"), pn, width = 6.7, height = 8.5)

pn1 <- p3|p4|p5|p10
ggsave(glue::glue("{outdir}/sfig4_p2.pdf"), pn1, width = 6.7, height = 2.3)

pn2 <- p6|p7|p8|p9$plot
ggsave(glue::glue("{outdir}/sfig4_p3.pdf"), pn2, width = 7.7, height = 2.5)


