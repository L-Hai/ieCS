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

Trying an online example of ieCS (for testing purpose only):
https://ling-hai.shinyapps.io/iecs/

## Graphical user interface of ieCS
<img src="https://user-images.githubusercontent.com/58597309/125817766-d0c29baa-1e1d-4a7d-8ffb-9a346a78b99d.png" alt="drawing" width="600"/>

## Input for ieCS
ieCS required markers of cell clusters in each sample and the rank of markers (e.g., fold change, p-value) as input. The markers need to be identified in single sample which can reduce the batch effect. ieCS provides a novel similarity metric to quantify pairwise similarity of cell clusters and identify superclusters within sample or across samples. 

The demo input can be downloaded from https://github.com/L-Hai/ieCS/tree/main/ExampleData. The demo input is generated from scRNA-seq dataset of PBMCs under stimulated and control conditions from Kang et al., 2018.

<img src="https://user-images.githubusercontent.com/58597309/125816599-0fb5807c-224c-4f38-ba92-3284e59b0f02.png" alt="drawing" width="300"/>


## Example results

Heatmap of pairwise similarity scores

<img src="https://user-images.githubusercontent.com/58597309/125818505-555b6b14-db3f-464e-a399-bd4324475afc.png" alt="drawing" width="600"/>

Identification of superclusters by hierarchical clustering

<img src="https://user-images.githubusercontent.com/58597309/125819365-3da277c9-67f5-4a83-894a-7730c7aa638e.png" alt="drawing" width="600"/>

Identification of superclusters by network partitioning

<img src="https://user-images.githubusercontent.com/58597309/125819831-8baf8423-df09-45e2-8f9b-199652fc9846.png" alt="drawing" width="600"/>

Identification of superclusters by tree aggregation

<img src="https://user-images.githubusercontent.com/58597309/125820196-2e1daf58-5b01-4d75-8997-2b8c6ce341a5.png" alt="drawing" width="600"/>
