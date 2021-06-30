# ieCS
ieCS: Interactive Explorer of Single Cell Cluster Similarity

ieCS is a ShinyApp that aids to explore similarity of cell clusters 
    in single cell dataset interactively. ieCS implements three methods to 
    identify superclusters of similar cell clusters within a sample 
    and also across different samples based on a novel similarity metric.

## Installation
Depends: 
    R (>= 3.5.0)<br />
Imports: 
    ape,
    collapsibleTree,
    DT,
    factoextra,
    ggplot2,
    ggraph,
    grid,
    gridExtra,
    igraph,
    pheatmap,
    plotly,
    plyr,
    RColorBrewer,
    shiny,
    stats,
    tidygraph,
    methods,
    cluster,
    colorspace,
    dendextend,
    ggpubr
``` r
devtools::install_github("L-Hai/ieCS")
ieCS::run()
```
## Example
