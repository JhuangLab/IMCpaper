pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "SingleCellExperiment", 
          "RColorBrewer", "vroom", "jhtools", "glue", "jhuanglabHyperion", "openxlsx", "ggsci", "ggraph",
          "patchwork", "cowplot", "tidyverse", "dplyr", "rstatix", "magrittr", "igraph","tidygraph")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "hyperion"
dataset <- "qzhang"
species <- "human"
workdir <- glue::glue("~/projects/{project}/analysis/{dataset}/{species}/workdir")
setwd(workdir)

#dir define
sce_dir <- "/cluster/home/jhuang/projects/hyperion/analysis/qzhang/human/steinbock/rds/combined"
meta_dir <- "/cluster/home/jhuang/projects/hyperion/analysis/qzhang/human/steinbock/rds/combined/lyj"
ci_dir <- "/cluster/home/yjliu_jh/share"

#data prepare
sce <- readr::read_rds(glue::glue("{sce_dir}/all_anno.rds"))
meta_clu <- readr::read_rds(glue::glue("{meta_dir}/meta_clu3.rds"))
ci <- readr::read_rds(glue::glue("{ci_dir}/ci_list.rds"))
col_new <- colData(sce) %>% as_tibble()
col_meta <- unique(meta_clu$meta_color) %>% `names<-`(unique(meta_clu$meta_short_new))
metadata <- metadata(sce)
metadata[["differentiation_degree"]]$diff_degree <- factor(metadata[["differentiation_degree"]]$diff_degree, levels = c("low", "middle", "high"))

metadata[["stages_tme"]] <- metadata[["stages_tme"]] %>% dplyr::mutate(stage = case_when(stage == 1 ~ "RPC",
                                                                                         stage == 2 ~ "BRPC_LAPC",
                                                                                         stage == 3 ~ "MPC"))
levels(metadata[["stages_tme"]]$stage) = c("MPC", "BRPC_LAPC", "RPC")
keycols_list <- c(list("stage"), list("treatment_type"), list("os_group_24"))
compare_groups <- c("stages_tme", "neoadj_vs_direct_surgery", "os_analysis")
names(keycols_list) <- compare_groups
meta_cluster_names <- rev(c("MC-M1M2-mixed", "MC-immune-enriched", "MC-M1-like", "MC-immune-mixed",
                            "MC-stroma-CAF", "MC-stroma-macro",
                            "MC-stroma-mCAF", "MC-tumor-boundary", "MC-tumor-core"))
total_group <- col_new %>% dplyr::select(c(sample_id, sample_tiff_id, patient_id, stype2, meta_cluster_new)) %>%
  na.omit() |> dplyr::filter(meta_cluster_new %in% meta_cluster_names)
#---------------------------------------- main figures fig3a----------------------------------------
plist <- list()
# total meta cluster ----------------------------------------
total_group_pro <- total_group %>% group_by(meta_cluster_new, sample_tiff_id) %>% summarise(pro = n()/nrow(total_group))
total_group_pro$meta_cluster_new <- fct_relevel(total_group_pro$meta_cluster_new, meta_cluster_names)
p1 <- ggplot(total_group_pro, aes(x = pro, y = meta_cluster_new, fill = meta_cluster_new)) +
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_manual(values = col_meta,
                    labels = vars(meta_cluster_new)) +
  theme(strip.placement  = "outside",
        panel.spacing    = unit(3, "points"),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 5),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.line.x = element_line(size = 0.4),
        axis.line.y = element_line(size = 0.4),
        legend.position="none")
ggsave("community_total_metacluster.pdf", p1, width = 3, height = 4)
plist[["total"]] <- p1
# compare groups ----------------------------------------
# -------- all total for alteration? and for new p1  --------
all_total <- total_group_pro %>% group_by(meta_cluster_new) %>% summarise(allpro = sum(pro))
factor_list <- list("stages_tme" = c("MPC", "BRPC_LAPC", "RPC"), 
                    "neoadj_vs_direct_surgery" = c("direct_surgery", "surgery_after_chemo"),
                    "os_analysis" = c("short", "long"))
# compare groups ----------------------------------------
for(compare_group in compare_groups){
  cli::cli_h1(compare_group)
  groups <- metadata[[compare_group]]
  dt_groups <- total_group %>% dplyr::filter(patient_id %in% groups$patient_id)
  if (compare_group == "chemo_outcome_before") {
    dt_groups <- dt_groups %>% dplyr::filter(stype2 %in% c("before_chemo","puncture_pdac"))
  } else if (compare_group %notin% c("paired_after_neoadj", "mpc_primary_metastasis", "tumor_para")) {
    dt_groups <- dt_groups %>% dplyr::filter(stype2 %in% c("after_chemo","tumor"))
  }
  dat_plot <- dt_groups %>%
    group_by(meta_cluster_new, sample_tiff_id) |>
    summarise(nc = n()) |>
    group_by(sample_tiff_id) |>
    dplyr::mutate(nt = sum(nc)) |>
    dplyr::mutate(pro = nc/nt) |> ungroup()
  dat_plot$meta_cluster_new <- factor(dat_plot$meta_cluster_new, levels = meta_cluster_names)
  dat_plot <- dat_plot |> tidyr::complete(meta_cluster_new, sample_tiff_id, fill = list(pro = 0)) |>
    dplyr::mutate(pro = sd_censor(pro, range = c(-3, 3)))
  metainfo <- col_new[, c("sample_tiff_id", "patient_id")] |> distinct()
  dat_plot <- left_join(dat_plot, metainfo, by = "sample_tiff_id") |>
    left_join(groups, by = "patient_id") %>% na.omit() |> distinct()
  
  gp_key <- keycols_list[[compare_group]]
  cli::cli_alert_info(gp_key)
  dat_plot[[gp_key]] <- factor(dat_plot[[gp_key]], levels = factor_list[[compare_group]])
  dat_plot$meta_cluster_new <- fct_relevel(dat_plot$meta_cluster_new, meta_cluster_names)
  p <- ggboxplot(dat_plot, x = "meta_cluster_new", y = "pro", fill = gp_key, outlier.shape = NA,
                 palette = pal_nejm("default")(3), xlab = NULL,size = 0.2) + theme_classic() +
    theme(strip.placement  = "outside",
          panel.spacing    = unit(3, "points"),
          strip.background = element_blank(),
          strip.text       = element_text(face = "bold", size = 5),
          axis.text.x = element_text(size = 10), axis.text.y = element_blank(),    
          axis.line.x = element_line(size = 0.4),
          axis.line.y = element_line(size = 0.4),
          legend.position="right") +
    labs(x= NULL, y = NULL)
  exp1 <- expr(pro ~ !!ensym(gp_key))
  stat_test <- dat_plot |>
    group_by(meta_cluster_new) %>% rstatix::t_test(eval(exp1), p.adjust.method = "none")
  stat_test <- stat_test |> mutate(p.adj.signif = case_when(p >= 0.05 ~ "ns",
                                                            p >= 0.01 & p < 0.05 ~ "*",
                                                            p >= 0.001 & p < 0.01 ~ "**",
                                                            p >= 0.0001 & p < 0.001 ~ "***",
                                                            p < 0.0001 ~ "****",
                                                            TRUE ~ "ns"))
  stat_test <- stat_test %>%
    add_xy_position(x = "meta_cluster_new", dodge = 0.8)
  p1 <- p +
    stat_pvalue_manual(
      stat_test, tip.length = 0.01, hide.ns = T, label = "p.adj.signif",
      coord.flip = TRUE
    ) + coord_flip() +
    guides(fill = guide_legend(reverse = TRUE))
  plist[[compare_group]] <- p1
}
pc <- (plist[[1]]|plist[[2]]|plist[[3]]|plist[[4]]) + plot_layout(guides = 'collect', widths = c(2, 2.5, 2.5, 2.5))
pc
ggsave("fig3abcd.pdf", pc, width = 10, height = 4)


# network -----------------------------------------------------------------
#function
interaction_plot <- function(x, conditionCol, group, threshold = 0.1, selfDe = T){
    from <- "from_meta"
    to <- "to_meta"
    fn <- "metacluster_"
    cols <- col_meta
  
  if (selfDe == TRUE) {
    fil <- "filteredself"
    dt_cl <- x %>% dplyr::filter(.data[[from]] != .data[[to]]) %>% group_by(.data[[from]], .data[[to]], .data[[conditionCol]]) %>% summarise(N = n()) %>%
      left_join(x %>% dplyr::filter(.data[[from]] != .data[[to]]) %>% group_by(.data[[conditionCol]]) %>% summarise(NT = n()), by = conditionCol) %>% dplyr::mutate(pct = N/NT * 100)
  }else {
    fil <- ""
    dt_cl <- x %>% group_by(.data[[from]], .data[[to]], .data[[conditionCol]]) %>% summarise(N = n()) %>%
      left_join(x %>% group_by(.data[[conditionCol]]) %>% summarise(NT = n()), by = conditionCol) %>% dplyr::mutate(pct = N/NT * 100)
  }
  
  dt_cl <- dt_cl %>% group_by(.data[[conditionCol]])
  glist <- dt_cl %>% dplyr::filter(pct > 0.1) %>% group_map(~as_tbl_graph(dplyr::select(.,c(.data[[from]], .data[[to]], pct))))
  gkeys <- dt_cl %>% group_keys()
  plist <- list()
  for(i in 1:length(glist)){
    g <- glist[[i]]
    if(length(g) == 0) next
    names <- paste0(fn, group, "_", gkeys[i,][[conditionCol]], "_", fil)
    plist[[names]] <- ggraph(g, "stress") +
      geom_edge_link(aes(edge_width = pct, alpha = pct)) +
      geom_node_point(aes(color = name), size = 3) +
      #geom_node_text(aes(label = name), size = 4, repel = T, colour = "red") +
      scale_edge_width("interaction percentage", range = c(0.2, 1.5), breaks = c(0.1,1,3,5,10,15,20,25),limits=c(0.1, 25), oob = scales::squish,
                       guide = guide_legend(title.theme = element_text(size = 8), ncol = 1, byrow = FALSE, keywidth = 0.8, keyheight = 0.1)) +
      scale_edge_alpha("interaction percentage", range = c(0.1, 1.5), breaks = c(0.1,1,3,5,10,15,20,25),limits=c(0.1, 25), oob = scales::squish) +
      scale_color_manual(limits = as.factor(V(g)$name), values = cols) +
      theme_graph(base_family = 'Arial', base_size =8) +
      guides(col = guide_legend(title = "interaction", ncol = 1, byrow = FALSE, keywidth = 0.1, keyheight = 0.1,
                                override.aes = list(size=1),title.theme = element_text(size = 8)))
  }
  return(plist)
}

col_meta <- unique(meta_clu$meta_color) %>% `names<-`(unique(meta_clu$meta_short))

# clinical metadata
keycols_list <- c(list("stage"), list("treatment_type"), list("stype2"),
                  list(c("size_alteration", "outcome_group")),
                  list(c("size_alteration", "outcome_group")),
                  list("stype2"), list("pfs_group"),
                  list(c("os_group_12", "os_group_18", "os_group_24")),
                  list("stype2"), list(c("relapse_site_group")),
                  list("hyphology_type"), list("nerve_invasion"), list("vascular_invasion"),
                  list(c("lymph_metastasis_degree", "lymph_metastasis_status")), list("diff_degree"))

names(keycols_list) <- names(metadata(sce))[3:17]

keycols_list <- keycols_list[c("stages_tme", "neoadj_vs_direct_surgery", "pfs_analysis",
                               "os_analysis", "relapse_analysis", "rare_hyphology", "nerve_invasion", "vascular_invasion", "lymph_metastasis", "differentiation_degree")]

dt_list <- list()
for (i in 1:length(keycols_list)) {
  dt_list[[i]] <- metadata(sce)[[names(keycols_list[i])]][c("patient_id", keycols_list[[i]])]
  dt_list[[i]] <- col_new %>% dplyr::filter(patient_id %in% dt_list[[i]]$patient_id) %>%
    dplyr::select(sample_id, stype2, patient_id) %>% unique() %>% left_join(dt_list[[i]], by = "patient_id") %>% dplyr::filter(stype2 %in% c("after_chemo", "tumor"))
}

for (i in 1:length(dt_list)) {
  dt_list[[i]] <- left_join(dt_list[[i]] %>% dplyr::select(-patient_id), col_new %>% dplyr::select(sample_id, sample_tiff_id) %>% unique(), by = "sample_id")
  dt_list[[i]] <- left_join(dt_list[[i]], do.call("rbind", ci), by = "sample_tiff_id") %>% dplyr::select(-c(from, to))
}

p_list <- list()
for (i in 1:length(dt_list)) {
  for (j in 1:length(keycols_list[[i]])) {
    group1 <- str_c(names(keycols_list)[i],keycols_list[[i]][j], sep= "_")
    p_list[[group1]] <- interaction_plot(dt_list[[i]] %>% dplyr::select(sample_id, keycols_list[[i]][j], sample_tiff_id, from_cluster, to_cluster, from_meta, to_meta) %>% na.omit(),
                                          conditionCol = keycols_list[[i]][j], group = group1)
  }
}
p_mian <- p_list[[1]][[1]] + p_list[[1]][[2]] + p_list[[1]][[3]] +
  p_list[[2]][[1]] + p_list[[2]][[2]] + p_list[[6]][[1]] + p_list[[6]][[2]] + 
  plot_layout(design = "123#\n4567", guides='collect')
ggsave(p_mian,filename = "fig3efg_network.pdf",width = 12,height = 6)


