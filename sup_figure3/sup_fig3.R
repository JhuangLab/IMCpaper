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

keycols_list <- c(list("nerve_invasion"), list("vascular_invasion"),
                  list("lymph_metastasis_status"), list("size_alteration"),
                  list("size_alteration"), list("pfs_group"))
compare_groups <-c("nerve_invasion", "vascular_invasion", "lymph_metastasis",
                   "chemo_outcome_before", "chemo_outcome_after", "pfs_analysis")
names(keycols_list) <- compare_groups

meta_cluster_names <- rev(c("MC-M1M2-mixed", "MC-immune-enriched", "MC-M1-like", "MC-immune-mixed",
                            "MC-stroma-CAF", "MC-stroma-macro",
                            "MC-stroma-mCAF", "MC-tumor-boundary", "MC-tumor-core"))
total_group <- col_new %>% dplyr::select(c(sample_id, sample_tiff_id, patient_id, stype2, meta_cluster_new)) %>%
  na.omit() |> dplyr::filter(meta_cluster_new %in% meta_cluster_names)

factor_list <- list("nerve_invasion" = c("1", "0"),
                    "vascular_invasion" = c("1", "0"),
                    "lymph_metastasis" = c("1", "0"),
                    "chemo_outcome_before" = c("increase", "decrease"),
                    "chemo_outcome_after" = c("increase", "decrease"),
                    "pfs_analysis" = c("short","long"))

plist <- list()
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
          axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),    
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
  ggsave(glue("metacluster_{compare_group}.pdf"), p1, width = 6, height = 6)
  plist[[compare_group]] <- p1
}

# sfig3-1
pc1 <- (plist[[1]]/plist[[2]]/plist[[3]]) + plot_layout(widths = c(2, 2, 2))
ggsave("sfig3-1.pdf", pc1, width = 6, height = 10)

# sfig3-2
pc2 <- (plist[[4]]/plist[[5]]/plist[[6]]) + plot_layout(widths = c(2, 2, 2))
ggsave("sfig3-2.pdf", pc2, width = 16, height = 6)



# network -----------------------------------------------------------------
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
#mpc_primary_metastasis
dt_mpc <- metadata(sce)[["mpc_primary_metastasis"]]
dt_mpc <- col_new %>% dplyr::filter(patient_id %in% dt_mpc$patient_id) %>% dplyr::select(sample_id, stype2, patient_id) %>% unique()
dt_mpc <- dt_mpc %>% dplyr::filter(!(patient_id %in% (table(dt_mpc$patient_id)[table(dt_mpc$patient_id) == 1] %>% names())))

#paired_after_neoadj
dt_neoadj <- metadata(sce)[["paired_after_neoadj"]]
dt_neoadj <- col_new %>% dplyr::filter(patient_id %in% dt_neoadj$patient_id) %>% dplyr::select(sample_id, stype2, patient_id) %>% unique()
dt_neoadj <- dt_neoadj %>% dplyr::filter(!(patient_id %in% (table(dt_neoadj$patient_id)[table(dt_neoadj$patient_id) == 1] %>% names())))

#tumor_para
dt_para <- metadata(sce)[["tumor_para"]]
dt_para <- col_new %>% dplyr::filter(patient_id %in% dt_para$patient_id) %>% dplyr::select(sample_id, stype2, patient_id) %>% unique()
dt_para <- dt_para %>% dplyr::filter(!(patient_id %in% (table(dt_para[dt_para$stype2 != "normal",]$patient_id)[table(dt_para[dt_para$stype2 != "normal",]$patient_id) == 1] %>% names())))

#chemo_outcome_before
dt_cob <- metadata(sce)[["chemo_outcome_before"]][c("patient_id", "size_alteration")]
dt_cob <- col_new %>% dplyr::filter(patient_id %in% dt_cob$patient_id) %>%
  dplyr::select(sample_id, stype2, patient_id) %>% unique() %>% left_join(dt_cob, by = "patient_id") %>% dplyr::filter(stype2 %in% c("before_chemo", "puncture_pdac"))

#chemo_outcome_after
dt_coa <- metadata(sce)[["chemo_outcome_after"]][c("patient_id", "size_alteration")]
dt_coa <- col_new %>% dplyr::filter(patient_id %in% dt_coa$patient_id) %>%
  dplyr::select(sample_id, stype2, patient_id) %>% unique() %>% left_join(dt_coa, by = "patient_id") %>% dplyr::filter(stype2 %in% c("after_chemo", "tumor"))

dt_list <- list(dt_mpc, dt_neoadj, dt_para, dt_cob, dt_coa)
for (i in 1:length(dt_list)) {
  dt_list[[i]] <- left_join(dt_list[[i]] %>% dplyr::select(-patient_id), col_new %>% dplyr::select(sample_id, sample_tiff_id) %>% unique(), by = "sample_id")
  dt_list[[i]] <- left_join(dt_list[[i]], do.call("rbind", ci), by = "sample_tiff_id") %>% dplyr::select(-c(from, to))
}

names(dt_list) <- c("mpc_primary_metastasis", "paired_after_neoadj", "tumor_para", "chemo_outcome_before", "chemo_outcome_after")

group_list <- list(c("mpc_primary_metastasis", "stype2"), c("paired_after_neoadj", "stype2"),
                   c("tumor_para", "stype2"), c("chemo_outcome_before", "size_alteration"), c("chemo_outcome_after", "size_alteration"))

p1 <- purrr::map2(dt_list, group_list, ~interaction_plot(x = .x, conditionCol = .y[2], group = .y[1]))

p_before <- p1$chemo_outcome_before[[1]] + p1$chemo_outcome_before[[2]] + 
  plot_layout(guides='collect')
p_after <- p1$chemo_outcome_after[[1]] + p1$chemo_outcome_after[[2]] + 
  plot_layout(guides='collect')
p_group_12 <- p_list$os_analysis_os_group_12[[1]] + p_list$os_analysis_os_group_12[[2]] + 
  plot_layout(guides='collect')
ggsave(p_before, filename = "p_before.pdf",width = 6,height = 3)
ggsave(p_after, filename = "p_after.pdf",width = 6,height = 3)
ggsave(p_group_12, filename = "p_group_12.pdf",width = 6,height = 3)

p_nerve <- p_list$nerve_invasion_nerve_invasion[[1]] + p_list$nerve_invasion_nerve_invasion[[2]] + 
  plot_layout(guides='collect')
p_vascular <- p_list$vascular_invasion_vascular_invasion[[1]] + p_list$vascular_invasion_vascular_invasion[[2]] + 
  plot_layout(guides='collect')
p_lymph <- p_list$lymph_metastasis_lymph_metastasis_status[[1]] + p_list$lymph_metastasis_lymph_metastasis_status[[2]] + 
  plot_layout(guides='collect')
ggsave(p_nerve, filename = "p_nerve.pdf",width = 6,height = 3)
ggsave(p_vascular, filename = "p_vascular.pdf",width = 6,height = 3)
ggsave(p_lymph, filename = "p_lymph.pdf",width = 6,height = 3)


