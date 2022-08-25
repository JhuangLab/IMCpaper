pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "survminer", "survival", "BiocNeighbors",
          "vroom", "jhtools", "glue", "openxlsx", "ggsci", "patchwork", "cowplot", "tidyverse", "dplyr", "SingleCellExperiment")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "hyperion"
dataset <- "qzhang"
workdir <- glue::glue("~/projects/{project}/output/{dataset}/human/figure5_v1")
setwd(workdir)

datadir <- glue::glue("~/projects/{project}/analysis/{dataset}/human/polaris_new/inform_batch_0718")

grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

#this version only use polaris data of slide1 to analyze.
#slide1 + 3_Scan1_8E + 3_Scan1_10H

# community size compare --------------------------------------------------
polaris_comsize <- readr::read_csv(glue::glue("{datadir}/com_preprocess/core_cell_community_dt_cl.csv"))
polaris_comsize_su <- polaris_comsize %>% na.omit() %>%
  dplyr::filter(slide %in% c("1_Scan2") | (slide == "3_Scan1" & Core_id %in% c("8E", "10H"))) %>%
  dplyr::select(community_id, n) %>%
  unique()

#read in data
sce <- readr::read_rds("/cluster/home/jhuang/projects/hyperion/analysis/qzhang/human/steinbock/rds/combined/all_anno.rds")
cld <- colData(sce) %>% as_tibble()
hyperion_comsize_su <- cld %>% dplyr::filter(stype2 %in% c("tumor", "after_chemo")) %>%
  na.omit() %>% group_by(community_name) %>%
  summarise(n = n())

#fig4 density plot of community size compare
comsize_su_com <- rbind(polaris_comsize_su %>% dplyr::rename(group = community_id) %>% dplyr::mutate(group = "polaris"),
                        hyperion_comsize_su %>% dplyr::rename(group = community_name) %>% dplyr::mutate(group = "hyperion"))

p1 <- ggdensity(comsize_su_com, x = "n", xlab = "community size",
               add = "mean", rug = F,
               color = "group", fill = "group",
               palette = ggsci::pal_nejm("default")(2)) +
  theme(axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.2, 'cm'))
ggsave("hyperion_polaris_community_size_compare.pdf", p1, height = 3, width = 3)
readr::write_rds(p1, "hyperion_polaris_community_size_compare.rds")


# os analysis -------------------------------------------------------------
surv_data_total <- readr::read_rds(glue::glue("{datadir}/surv_data_total.rds"))

#read in community results
meta_cluster_definition <- readr::read_csv(glue::glue("{datadir}/ana_community/meta_cluster_definition_1QuMax.csv"))
meta_cluster_definition <- janitor::clean_names(meta_cluster_definition)
meta_cluster_definition <- meta_cluster_definition %>% dplyr::mutate(roi = paste0(slide, "_", core_id))

meta_cluster_definition <- meta_cluster_definition %>%
  dplyr::mutate(meta_merge = case_when(meta_cl == "TB" ~ "TB", TRUE ~ "Others")) %>% 
  dplyr::filter(slide == "1_Scan2" | roi %in% c("3_Scan1_8E", "3_Scan1_10H"))

cell_type_meta_total <- meta_cluster_definition %>%
  group_by(roi, meta_merge) %>% summarise(nt = n()) %>% ungroup()

cell_type_meta <- meta_cluster_definition %>%
  group_by(roi, phenotype, meta_merge) %>% summarise(nc = n()) %>% ungroup()

cell_type_meta$phenotype <- factor(cell_type_meta$phenotype)
cell_type_meta <- cell_type_meta %>% tidyr::complete(roi, phenotype, meta_merge, fill = list(nc = 0))
cell_type_meta <- left_join(cell_type_meta_total, cell_type_meta, by = c("roi", "meta_merge")) %>%
  dplyr::mutate(prop = nc/nt) %>% dplyr::select(-c(nt, nc))

#add os info
cell_type_meta <- inner_join(cell_type_meta, surv_data_total %>% dplyr::filter(neo_adj == 0), by = "roi") %>%
  dplyr::select(roi, phenotype, meta_merge, prop, os_state, os_month)

#plot survival
cell_type_meta_tmp <- cell_type_meta %>% dplyr::filter(phenotype == "Macrophage")
psurvx_list <- list()
for (j in c("TB", "Others")) {
  cell_type_meta_tmp1 <- cell_type_meta_tmp %>% dplyr::filter(meta_merge == j) %>%
    dplyr::mutate(group_43 = case_when(prop >= quantile(prop, probs = 0.75) ~ "high", prop < quantile(prop, probs = 0.75) ~ "low"))
  
  cell_type_meta_tmp1_os <- cell_type_meta_tmp1 %>% dplyr::select(roi, os_state, os_month, group_43) %>% na.omit()
  
  #3/4 group
  psurvx_list[[j]] <- ggsurvplot(surv_fit(Surv(os_month, os_state) ~ group_43, data = cell_type_meta_tmp1_os), palette = ggsci::pal_nejm("default")(2),
                       legend.labs = levels(droplevels(as.factor(unlist(cell_type_meta_tmp1_os[, "group_43"])))),
                       pval=T, risk.table = F, xlim = c(0,75))
  ggsave(glue::glue("cell_type_{j}_Macrophage_os_group43.pdf"),
         plot = psurvx_list[[j]], width = 2, height = 2)
}
readr::write_rds(psurvx_list, "polaris_macrophage_os_TB_and_Others.rds")


# DN-mac and tumor distance compare in TB and others ----------------------
meta_cluster_definition <- readr::read_csv(glue::glue("{datadir}/ana_community/meta_cluster_definition_1QuMax.csv"))
meta_cluster_definition <- janitor::clean_names(meta_cluster_definition)
meta_cluster_definition <- meta_cluster_definition %>%
  dplyr::mutate(roi = paste0(slide, "_", core_id)) %>%
  dplyr::mutate(meta_merge = case_when(meta_cl == "TB" ~ "TB", TRUE ~ "Others"))

cell_ids <- meta_cluster_definition %>%
  dplyr::filter(slide == "1_Scan2" | roi %in% c("3_Scan1_8E", "3_Scan1_10H")) %>% .$cell_id
k = 1

# Macrophage to tumor -----------------------------------------------------
target_cell <- "Tumor"
query_cell <- "Macrophage"

list_closecell_distance <- readr::read_rds(glue::glue("{datadir}/distance/TB_list_closecell_dist_k{k}_{query_cell}_to_{target_cell}.rds"))
df_closecell_distance_TB <- do.call("rbind", list_closecell_distance)
df_closecell_distance_TB$group <- "Tumor_boundary"

list_closecell_distance <- readr::read_rds(glue::glue("{datadir}/distance/Others_list_closecell_dist_k{k}_{query_cell}_to_{target_cell}.rds"))
df_closecell_distance_others <- do.call("rbind", list_closecell_distance)
df_closecell_distance_others$group <- "Others"

df_closecell_distance <- rbind(df_closecell_distance_TB, df_closecell_distance_others)
df_closecell_distance <- df_closecell_distance %>% dplyr::filter(from_cell %in% cell_ids)

p2 <- ggboxplot(df_closecell_distance, x = "group", y = "distance",
                color = "group",
                outlier.shape = NA,
                size = 0.2, palette = ggsci::pal_nejm("default")(2)) + ylim(0,100) +
  stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test", label.y = c(85, 90, 95), size = 1) +
  stat_summary(fun.data = function(x) data.frame(y = 70, label = paste("median=", round(median(x),2))), geom="text", size = 1) +
  theme(axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm')) + rotate_x_text(30)
ggsave(glue::glue("{query_cell}_to_{target_cell}_dist_compareTBandOthers.pdf"), p2, height = 3, width = 2)
readr::write_rds(p2, glue::glue("{query_cell}_to_{target_cell}_dist_compareTBandOthers.rds"))
