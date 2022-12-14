# Code

# Single-Cell-Based Clinically-Related Mesoscale Pathology of Pancreatic Cancer

This repository contains all code used to produce the results and figures of the publication "Single-Cell-Based Clinically-Related Mesoscale Pathology of Pancreatic Cancer".

Visulization for figures:

-   `figure1/figure1.R`

    code for plots of heatmap, barplots and tSNE in figure1.

-   `sup_figure1/sup_fig1.R`

    code for boxplots in sfigure2.

-   `figure2/figure2.R`

    code for dotplot, tSNE, circluar barplot and chord plot in figure2.

-   `figure3/figure3.R`

    code for barplots and network graphs in figure3.

-   `sup_figure3/sup_fig3.R`

    code for barplots and network graphs in supplementary figures related to figure3.

-   `figure4/figure4.R`

    code for interaction heatmaps, os analysis and boxplots in figure4.
    
-   `sup_figure4/sup_fig4.R`

    code for correlation dotplots, density plots, boxplots and os analysis in supplementary figures related to figure4.

-   `figure5/figure5.R`

    code for dotplots, density plot, os analysis and boxplots in figure5.
    
-   `image_plot.R` for all images of single ROI showed in figures.

    usage example:

        image_plot.R -s p113_normal -t normal

-   `cell_segmentation.cppipe` for cell segmentation pipeline of CellProfiler
