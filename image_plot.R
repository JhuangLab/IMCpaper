#!/usr/bin/env Rscript
pkgs <- c("fs", "futile.logger", "optparse", "configr", "stringr", "ggpubr", "ggthemes", "BiocParallel",
          "jhtools", "glue", "ggsci", "data.table", "tidyverse", "dplyr", "imcRtools",
          "cytomapper", "ggraph", "ggplot2", "Rtsne", "jhuanglabHyperion")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-s","--sample_id"),
              help="sample_id"),
  make_option(c("-t","--stype"),
              help="stype name")
)
opt <- parse_args(OptionParser(option_list=option_list))
sample_id <- opt$sample_id
stype <- opt$stype

conflicted::conflict_prefer("mutate", "dplyr")
conflicted::conflict_prefer("channelNames<-", "cytomapper")
set.seed(2022)
project <- "hyperion"
dataset <- "qzhang"
species <- "human"
wkdir <- "~/projects/hyperion/analysis/qzhang/human/steinbock/cytomapper/pixels/"
setwd(wkdir)

message("Processing sample: ", sample_id)

#psudo-color image
dat_dir <- "/cluster/home/jhuang/projects/hyperion/analysis/qzhang/human/steinbock"
#data prepare
panel_fn <- glue("~/projects/{project}/analysis/{dataset}/{species}/steinbock/panel.csv")
panel <- fread(panel_fn)
img_path <- file.path(dat_dir, stype, "img", sample_id)

# img
img <- loadImages(img_path, pattern = "*.tiff",
                  BPPARAM = MulticoreParam(8))
channelNames(img) <- panel$name
mcols(img) <- data.frame(sample_tiff_id = names(img))

# plot 
pdf_fn <- glue("{wkdir}/{sample_id}_pixels.pdf")
pdf(pdf_fn)
  plotPixels(img, img_id = "sample_tiff_id",
             colour_by = c("CD68", "HLA_DR", "DNA1", "Pan_Keratin", "CD163", "a_smooth"),
             bcg = list(CD68 = c(0, 8, 1),
                        HLA_DR = c(0, 8, 1),
                        DNA1 = c(0, 2, 1),
                        Pan_Keratin = c(0, 5, 1),
                        CD163 = c(0, 8, 1),
                        a_smooth = c(0, 8, 1)),
             colour = list(CD68 = c("black", "green"),
                           a_smooth = c("black", "red"),
                           Pan_Keratin = c("black", "magenta"),
                           HLA_DR = c("black", "yellow"),
                           CD163 = c("black", "cyan"),
                           DNA1 = c("black", "blue")),
             display = "single",
             image_title = list(cex = 0.3,
                                colour = "red"),
             legend = list(colour_by.labels.cex = 0.7,
                           colour_by.legend.cex = 0.7))
dev.off()




#sce_path <- "/cluster/home/jhuang/projects/hyperion/analysis/qzhang/human/steinbock/rds/combined/all_anno.rds"
#cl_path <- "/cluster/home/jhuang/projects/hyperion/analysis/qzhang/human/figures/heatmap/combined/annotation/cell_types.tsv"

#sce <- readr::read_rds(sce_path)
#dat <- colData(sce) %>% as.data.frame()
#dat <- dat %>% mutate(cell_type = ifelse(cell_type != "Unknown", "Identified", cell_type))
#colData(sce)$cell_type <- dat$cell_type
