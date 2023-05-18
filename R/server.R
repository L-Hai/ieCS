###############################################################################
### Set the maximum size of each file uploaded to ShinyApp
options(shiny.maxRequestSize = 30 * 1024 ^ 2) # maximum file size is 30MB
###############################################################################
### Shiny Server
server <- function(input, output, session) {
  ### Create working directory
  analysisID <- session$token
  workPath <- paste0('~/', analysisID, '/')
  dir.create(workPath)
  setwd(workPath)

  ###############################################################################
  ### Upload markers files
  observeEvent(input$MarkersSubmit, {
    if (!is.null(input$MarkersFile)) {
      apath <- input$MarkersFile$datapath
      df <- read.csv(apath, stringsAsFactors = FALSE)
      df <- data.frame(lapply(df, function(x) {
        gsub(" ", "_", x)
      }))
      write.csv(df, paste0(workPath, 'Que.Markers.csv'), row.names = TRUE)
      updateSelectInput(session,
                        "OrderM",
                        choices = colnames(df))
      updateSelectInput(session,
                        "ClusterM",
                        choices = colnames(df))
      updateSelectInput(session,
                        "GeneM",
                        choices = colnames(df))
    }
    if (!is.null(input$RefMarkersFile)) {
      apath <- input$RefMarkersFile$datapath
      df <- read.csv(apath, stringsAsFactors = FALSE)
      df <- data.frame(lapply(df, function(x) {
        gsub(" ", "_", x)
      }))
      write.csv(df, paste0(workPath, 'Ref.Markers.csv'), row.names = TRUE)
      updateSelectInput(session,
                        "OrderMR",
                        choices = colnames(df))
      updateSelectInput(session,
                        "ClusterMR",
                        choices = colnames(df))
      updateSelectInput(session,
                        "GeneMR",
                        choices = colnames(df))
    }
  })

  ###############################################################################
  ### Reactive values

  ### scoreMode 1,Query x Query; 2,(Query + Reference) x (Query + Reference); 3,Query x Reference
  reaValScore <-
    reactiveValues(scoreMode1 = NULL,
                   scoreMode2 = NULL,
                   scoreMode3 = NULL)

  ### Query markers table
  reaVal <- reactiveValues(treeObject = NULL, newickCS = NULL)
  reaValU <- reactiveValues(newickCS = NULL)
  reaValG <- reactiveValues(GraphObject = NULL, GraphCS = NULL)
  reaValGU <- reactiveValues(GraphCS = NULL)

  ### With Reference markers table
  reaValRef <- reactiveValues(treeObject = NULL, newickCS = NULL)
  reaValURef <- reactiveValues(newickCS = NULL)
  reaValGRef <- reactiveValues(GraphObject = NULL, GraphCS = NULL)
  reaValGURef <- reactiveValues(GraphCS = NULL)

  ###############################################################################
  ### Score cluster similarity with reference marker table
  observeEvent(input$RunScore, {
    df <-
      read.csv(
        paste0(workPath, 'Que.Markers.csv'),
        stringsAsFactors = FALSE,
        row.names = 1
      )
    df$cluster <- df[, input$ClusterM]
    df$gene <- df[, input$GeneM]
    df$order <- df[, input$OrderM]
    df <- df[order(df$order, decreasing = input$DirectionM),]
    write.csv(df, paste0(workPath, 'Que.Markers.csv'))
    reaValScore$scoreMode1 <- scoreCS(1)
    if (!is.null(input$RefMarkersFile)) {
      dfR <-
        read.csv(
          paste0(workPath, 'Ref.Markers.csv'),
          stringsAsFactors = FALSE,
          row.names = 1
        )
      dfR$cluster <- dfR[, input$ClusterMR]
      dfR$gene <- dfR[, input$GeneMR]
      dfR$order <- dfR[, input$OrderMR]
      dfR <- dfR[order(dfR$order, decreasing = input$DirectionMR),]
      write.csv(dfR, paste0(workPath, 'Ref.Markers.csv'))
      reaValScore$scoreMode2 <- scoreCS(2)
      reaValScore$scoreMode3 <- scoreCS(3)
      dfR$tag <- 'Reference'
      dfR <- dfR[, c('cluster', 'gene', 'order', 'tag')]
      df$tag <- 'Query'
      df <- df[, c('cluster', 'gene', 'order', 'tag')]
      df <- rbind(df, dfR)
    }
    output$contents <- DT::renderDataTable({
      df
    })
  })

  ###############################################################################
  ### Run cluster similarity interactively

  ### Tree aggregation method
  rea <- reactive({
    MIN <- input$minScore
    pair.score.matrix <- reaValScore$scoreMode1@pair.score.matrix
    marker.table <- reaValScore$scoreMode1@marker.table
    df <- as.matrix(pair.score.matrix)
    seqn <- seq(MIN, max(df), 1)
    reaValObj <- CSCoreFunc(marker.table, seqn, df)
    reaVal$treeObject <- reaValObj@treeObject
    reaVal$newickCS <- reaValObj@newickCS
  })

  reaU <- reactive({
    MIN <- input$minScoreU
    pair.score.matrix <- reaValScore$scoreMode1@pair.score.matrix
    marker.table <- reaValScore$scoreMode1@marker.table
    df <- as.matrix(pair.score.matrix)
    seqn <- seq(MIN, max(df), 1)
    reaValObj <- CSCoreFunc(marker.table, seqn, df)
    reaValU$newickCS <- reaValObj@newickCS
  })

  ### Network partitioning method
  reaG <- reactive({
    MIN <- input$minScoreG
    edges <- reaValScore$scoreMode1@edges
    nodes <- reaValScore$scoreMode1@nodes
    marker.table <- reaValScore$scoreMode1@marker.table
    edges <- edges[edges[, 'weight'] > MIN, ]
    GraphObject <-
      igraph::graph_from_data_frame(d = edges,
                                    vertices = nodes,
                                    directed = FALSE)
    reaValG$GraphObject <- GraphObject
    GraphGroupOri <-
      igraph::cluster_louvain(GraphObject, weights = igraph::E(GraphObject)$weight)
    GraphGroup <- list()
    for (i in 1:length(GraphGroupOri)) {
      ggroup <-
        plyr::mapvalues(GraphGroupOri[[i]], nodes$id, nodes$label, warn_missing = FALSE)
      GraphGroup[[paste0('LouvainGroup_', i)]] <-
        ggroup
    }
    GraphGene <-
      findOverlapGene(marker.table, GraphGroup)
    GraphCS <-
      CS(cluster = GraphGroup,
         gene = GraphGene)
    reaValG$GraphCS <- GraphCS
  })

  reaGU <- reactive({
    MIN <- input$minScoreGU
    edges <- reaValScore$scoreMode1@edges
    nodes <- reaValScore$scoreMode1@nodes
    marker.table <- reaValScore$scoreMode1@marker.table
    edges <- edges[edges[, 'weight'] > MIN, ]
    GraphObject <-
      igraph::graph_from_data_frame(d = edges,
                                    vertices = nodes,
                                    directed = FALSE)
    GraphGroupOri <-
      igraph::cluster_louvain(GraphObject, weights = igraph::E(GraphObject)$weight)
    GraphGroup <- list()
    for (i in 1:length(GraphGroupOri)) {
      ggroup <-
        plyr::mapvalues(GraphGroupOri[[i]], nodes$id, nodes$label, warn_missing = FALSE)
      GraphGroup[[paste0('LouvainGroup_', i)]] <-
        ggroup
    }
    GraphGene <-
      findOverlapGene(marker.table, GraphGroup)
    GraphCS <-
      CS(cluster = GraphGroup,
         gene = GraphGene)
    reaValGU$GraphCS <- GraphCS
  })

  ###############################################################################
  ### Run cluster similarity with reference markers table interactively

  ### Tree aggregation method
  reaRef <- reactive({
    MIN <- input$minScore
    pair.score.matrix <- reaValScore$scoreMode2@pair.score.matrix
    marker.table <- reaValScore$scoreMode2@marker.table
    df <- as.matrix(pair.score.matrix)
    seqn <- seq(MIN, max(df), 1)
    reaValObj <- CSCoreFunc(marker.table, seqn, df)
    reaValRef$treeObject <- reaValObj@treeObject
    reaValRef$newickCS <- reaValObj@newickCS
  })

  reaURef <- reactive({
    MIN <- input$minScoreU
    pair.score.matrix <- reaValScore$scoreMode2@pair.score.matrix
    marker.table <- reaValScore$scoreMode2@marker.table
    df <- as.matrix(pair.score.matrix)
    seqn <- seq(MIN, max(df), 1)
    reaValObj <- CSCoreFunc(marker.table, seqn, df)
    reaValURef$newickCS <- reaValObj@newickCS
  })

  ### Network partitioning method
  reaGRef <- reactive({
    MIN <- input$minScoreG
    edges <- reaValScore$scoreMode2@edges
    nodes <- reaValScore$scoreMode2@nodes
    marker.table <- reaValScore$scoreMode2@marker.table
    edges <- edges[edges[, 'weight'] > MIN, ]
    GraphObject <-
      igraph::graph_from_data_frame(d = edges,
                                    vertices = nodes,
                                    directed = FALSE)
    reaValGRef$GraphObject <- GraphObject
    GraphGroupOri <-
      igraph::cluster_louvain(GraphObject, weights = igraph::E(GraphObject)$weight)
    GraphGroup <- list()
    for (i in 1:length(GraphGroupOri)) {
      ggroup <-
        plyr::mapvalues(GraphGroupOri[[i]], nodes$id, nodes$label, warn_missing = FALSE)
      GraphGroup[[paste0('LouvainGroup_', i)]] <-
        ggroup
    }
    GraphGene <-
      findOverlapGene(marker.table, GraphGroup)
    GraphCS <-
      CS(cluster = GraphGroup,
         gene = GraphGene)
    reaValGRef$GraphCS <- GraphCS
  })

  reaGURef <- reactive({
    MIN <- input$minScoreGU
    edges <- reaValScore$scoreMode2@edges
    nodes <- reaValScore$scoreMode2@nodes
    marker.table <- reaValScore$scoreMode2@marker.table
    edges <- edges[edges[, 'weight'] > MIN, ]
    GraphObject <-
      igraph::graph_from_data_frame(d = edges,
                                    vertices = nodes,
                                    directed = FALSE)
    GraphGroupOri <-
      igraph::cluster_louvain(GraphObject, weights = igraph::E(GraphObject)$weight)
    GraphGroup <- list()
    for (i in 1:length(GraphGroupOri)) {
      ggroup <-
        plyr::mapvalues(GraphGroupOri[[i]], nodes$id, nodes$label, warn_missing = FALSE)
      GraphGroup[[paste0('LouvainGroup_', i)]] <-
        ggroup
    }
    GraphGene <-
      findOverlapGene(marker.table, GraphGroup)
    GraphCS <-
      CS(cluster = GraphGroup,
         gene = GraphGene)
    reaValGURef$GraphCS <- GraphCS
  })

  ###############################################################################
  ### Plot results from supercluster identification methods

  ###############################################################################
  ### Hierarchical clustering - global distances
  hideTab(inputId = "GeneralModeTab", target = "Detail: Query x Reference")
  hideTab(inputId = "GeneralModeTabD", target = "Dendrogram of reference")

  ### Visualize pair-wise similarity of clusters in Heatmap
  output$Heatmap <- renderPlot({
    if (input$HeatmapMode == 1) {
      p1 <- reaValScore$scoreMode1@p1
      matrix <- reaValScore$scoreMode1@pair.score.matrix
    }
    if (input$HeatmapMode == 2) {
      p1 <- reaValScore$scoreMode2@p1
      matrix <- reaValScore$scoreMode2@pair.score.matrix
    }
    if (input$HeatmapMode == 3) {
      p1 <- reaValScore$scoreMode3@p1
      matrix <- reaValScore$scoreMode3@pair.score.matrix
    }
    output$MDownMatrixT <- downloadHandler(
      filename = function() {
        'ModeMatrixTable.csv'
      },
      content = function(file) {
        write.csv(matrix, file)
      }
    )
    output$DownHeatmap <- downloadHandler(
      filename = function() {
        'Heatmap.pdf'
      },
      content = function(file) {
        pdf(file)
        print(p1)
        dev.off()
      }
    )
    return(p1)
  }, 1000, 1000)

  ### Dendrogram
  output$Hclust <- renderPlot({
    if (input$HeatmapMode == 1) {
      p1 <- reaValScore$scoreMode1@p1
    }
    if (input$HeatmapMode == 2) {
      p1 <- reaValScore$scoreMode2@p1
    }
    if (input$HeatmapMode == 3) {
      p1 <- reaValScore$scoreMode3@p1
    }
    hc <- p1$tree_row
    plot(hc , hang = -1)
    output$DownHclust <- downloadHandler(
      filename = function() {
        'Hclust.pdf'
      },
      content = function(file) {
        pdf(file)
        print(plot(hc , hang = -1))
        dev.off()
      }
    )
  }, 1000, 500)

  ### Dendrogram with user-defined number of cluster
  output$CutHclust <- renderPlot({
    if (input$HeatmapMode == 1) {
      p1 <- reaValScore$scoreMode1@p1
    }
    if (input$HeatmapMode == 2) {
      p1 <- reaValScore$scoreMode2@p1
    }
    if (input$HeatmapMode == 3) {
      p1 <- reaValScore$scoreMode3@p1
    }
    hc <- p1$tree_row
    output$DownCutHclust <- downloadHandler(
      filename = function() {
        'CutHclust.pdf'
      },
      content = function(file) {
        pdf(file)
        print(factoextra::fviz_dend(hc, rect = TRUE, k = input$NBCutHclust,labels_track_height = 100))
        dev.off()
      }
    )
    factoextra::fviz_dend(hc, rect = TRUE, k = input$NBCutHclust,labels_track_height = 100)
  },1000,500)

  ### Dendrogram with optimal number of cluster
  output$OptHclust <- renderPlot({
    if (input$HeatmapMode == 1) {
      df <- reaValScore$scoreMode1@pair.score.matrix
    }
    if (input$HeatmapMode == 2) {
      df <- reaValScore$scoreMode2@pair.score.matrix
    }
    if (input$HeatmapMode == 3) {
      df <- reaValScore$scoreMode3@pair.score.matrix
    }
    k.max <- nrow(df) - 1
    sil <-
      factoextra::fviz_nbclust(df,
                               FUN = factoextra::hcut,
                               method = "silhouette",
                               k.max = 24) + 
      ggplot2::labs(subtitle = "Silhouette method")
    n_clust <- sil$data
    max_cluster <-
      as.numeric(n_clust$clusters[which.max(n_clust$y)])
    figsize <- length(sil$data$clusters) / max(sil$data$y) / (3)
    sil <- sil + ggplot2::coord_fixed(figsize)
    wss <-
      factoextra::fviz_nbclust(df,
                               FUN = factoextra::hcut,
                               method = "wss",
                               k.max = k.max) + 
      ggplot2::geom_vline(xintercept = max_cluster , linetype = 2) + 
      ggplot2::labs(subtitle = "Elbow method")
    figsize <- length(wss$data$clusters) / max(wss$data$y) / (3)
    wss <- wss + ggplot2::coord_fixed(figsize)
    res <-
      factoextra::hcut(df, k = max_cluster, hc_method = "complete")
    # figsize <- nrow(df) / 200 / (3)
    dend <- factoextra::fviz_dend(res, rect = TRUE,labels_track_height = 100)
      # ggplot2::coord_fixed(figsize)
    gridExtra::grid.arrange(
      dend,
      sil,
      wss,
      nrow = 3,
      top = grid::textGrob(
        paste0("Optimal number of supercluster is ", max_cluster),
        gp = grid::gpar(fontsize = 20)
      )
    )
    output$DownOptHclust <- downloadHandler(
      filename = function() {
        'OptHclust.pdf'
      },
      content = function(file) {
        pdf(file)
        print(gridExtra::grid.arrange(
          dend,
          sil,
          wss,
          nrow = 3,
          top = grid::textGrob(
            paste0("Optimal number of supercluster is ", max_cluster),
            gp = grid::gpar(fontsize = 20)
          )
        ))
        dev.off()
      }
    )
  }, 1000, 1500)

  ### Query x Reference, assignment of superclusters with reference
  observeEvent(input$HeatmapMode, {
    if (input$HeatmapMode == 1) {
      hideTab(inputId = "GeneralModeTab", target = "Detail: Query x Reference")
    }
    if (!is.null(input$RefMarkersFile)) {
      if (input$HeatmapMode == 2) {
        hideTab(inputId = "GeneralModeTab", target = "Detail: Query x Reference")
      }
      if (input$HeatmapMode == 3) {
        showTab(inputId = "GeneralModeTab", target = "Detail: Query x Reference")
        p1 <- reaValScore$scoreMode3@p1
        matrix <- reaValScore$scoreMode3@pair.score.matrix
        marker.table <- reaValScore$scoreMode3@marker.table
        celltype.result <- c()
        for (i in 1:nrow(matrix)) {
          celltype.result <-
            c(celltype.result, colnames(matrix)[which.max(matrix[i, ])])
        }
        df <-
          data.frame(Cluster = rownames(matrix), Reference = celltype.result)
        Reference <- df
        Reference$Color <- as.factor(Reference$Reference)
        levels(Reference$Color) <-
          colorspace::rainbow_hcl(length(unique(Reference$Color)))
        coTree <- collapsibleTree::collapsibleTree(Reference,
                                                   c("Reference", "Cluster"),
                                                   collapsed = FALSE)
        coTreeGroup <- list()
        for (i in unique(df$Reference)) {
          coTreeGroup[[i]] <- as.character(df[df$Reference == i, 'Cluster'])
        }
        coTreeGene <-
          findOverlapGene(marker.table, coTreeGroup)
        output$HclustMode3 <- renderPlot({
          if (input$HeatmapMode == 3) {
            p1 <- reaValScore$scoreMode3@p1
            hc <- p1$tree_col
            plot(hc, hang = -1)
            output$DownHclustMode3 <- downloadHandler(
              filename = function() {
                'HclustMode3.pdf'
              },
              content = function(file) {
                pdf(file)
                print(plot(hc, hang = -1))
                dev.off()
              }
            )
          }
        }, 1000, 500)
        output$CollectMode3 <-
          collapsibleTree::renderCollapsibleTree({
            coTree
          })
        output$M3GeneList <- renderPrint({
          print(coTreeGene)
        })
        output$M3ClusterList <- renderPrint({
          print(coTreeGroup)
        })
        output$M3GeneListT <- DT::renderDataTable({
          df <- data.frame(Reference = character(), Gene = character())
          for (i in 1:length(coTreeGene)) {
            gl <- coTreeGene[i]
            adf <- data.frame(Reference = names(gl), Gene = gl[[1]])
            df <- rbind(df, adf)
          }
          output$M3DownGeneListT <- downloadHandler(
            filename = function() {
              'Mode3GeneTable.csv'
            },
            content = function(file) {
              write.csv(df, file)
            }
          )
          return(df)
        })
        output$M3ClusterListT <- DT::renderDataTable({
          df <- data.frame(SuperCluster = character(), Cluster = character())
          for (i in 1:length(coTreeGroup)) {
            gl <- coTreeGroup[i]
            adf <-
              data.frame(SuperCluster = names(gl),
                         Cluster = gl[[1]])
            df <- rbind(df, adf)
          }
          output$M3DownClusterListT <- downloadHandler(
            filename = function() {
              'Mode3ClusterTable.csv'
            },
            content = function(file) {
              write.csv(df, file)
            }
          )
          return(df)
        })
      }
    }
  })


  ###############################################################################
  ### Hierarchical clustering - direct distances
  hideTab(inputId = "GeneralModeTabD", target = "Dendrogram of reference")

  ### Dendrogram
  output$HclustD <- renderPlot({
    if (input$HeatmapModeD == 1) {
      matrix <- reaValScore$scoreMode1@pair.score.matrix
    }
    if (input$HeatmapModeD == 2) {
      matrix <- reaValScore$scoreMode2@pair.score.matrix
    }
    dis.direct <- apply(matrix, 1, function(x)
      (max(x) - x) / max(x))
    dis.direct <- t(dis.direct)
    dis.direct <- stats::as.dist(dis.direct)
    hc <- stats::hclust(dis.direct)
    plot(hc , hang = -1)
    output$DownHclustD <- downloadHandler(
      filename = function() {
        'HclustD.pdf'
      },
      content = function(file) {
        pdf(file)
        print(plot(hc , hang = -1))
        dev.off()
      }
    )
  }, 1000, 500)

  ### Dendrogram of user-defined number of clusters
  output$CutHclustD <- renderPlot({
    if (input$HeatmapModeD == 1) {
      matrix <- reaValScore$scoreMode1@pair.score.matrix
    }
    if (input$HeatmapModeD == 2) {
      matrix <- reaValScore$scoreMode2@pair.score.matrix
    }
    dis.direct <- apply(matrix, 1, function(x)
      (max(x) - x) / max(x))
    dis.direct <- t(dis.direct)
    dis.direct <- stats::as.dist(dis.direct)
    hc <- stats::hclust(dis.direct)
    d1 <- stats::as.dendrogram(hc)
    d2 <- dendextend::color_branches(d1, k = input$NBCutHclustD)
    d3 <- dendextend::color_labels(d2, k = input$NBCutHclustD)
    par(mar = c(15,3,1,1))
    plot(d3)
    stats::rect.hclust(hc, k = input$NBCutHclustD, border = "grey")
    output$DownCutHclustD <- downloadHandler(
      filename = function() {
        'CutHclustD.pdf'
      },
      content = function(file) {
        pdf(file)
        print(plot(d3))
        print(stats::rect.hclust(hc, k = input$NBCutHclustD, border = "grey"))
        dev.off()
      }
    )
  }, 1000, 500)

  ### Dendrogram of optimal number of clusters
  output$OptHclustD <- renderPlot({
    if (input$HeatmapModeD == 1) {
      matrix <- reaValScore$scoreMode1@pair.score.matrix
    }
    if (input$HeatmapModeD == 2) {
      matrix <- reaValScore$scoreMode2@pair.score.matrix
    }
    dis.direct <- apply(matrix, 1, function(x)
      (max(x) - x) / max(x))
    dis.direct <- t(dis.direct)
    dis.direct <- stats::as.dist(dis.direct)
    hc <- stats::hclust(dis.direct)
    k.max <- nrow(matrix) - 1
    avg.slihouette <- 0
    for (i in 2:k.max) {
      cutc <- stats::cutree(hc, i)
      sil <- cluster::silhouette(cutc, dist = dis.direct)
      avg.slihouette <- c(avg.slihouette, mean(sil[, 3]))
    }
    adf <-
      data.frame(clusters = as.factor(1:k.max), y = avg.slihouette)
    p <- ggpubr::ggline(
      adf,
      x = "clusters",
      y = "y",
      group = 1,
      ylab = "Average Slihouette width",
      xlab = "Number of clusters k"
    )
    n.opt <- which.max(avg.slihouette)
    p <- p + ggplot2::geom_vline(xintercept = n.opt, linetype = 2)
    dend <- factoextra::fviz_dend(hc, rect = TRUE, k = n.opt)
    gridExtra::grid.arrange(dend,
                            p,
                            nrow = 2,
                            top = grid::textGrob(
                              paste0("Optimal number of supercluster is ", n.opt),
                              gp = grid::gpar(fontsize = 20)
                            ))
    output$DownOptHclustD <- downloadHandler(
      filename = function() {
        'OptHclustD.pdf'
      },
      content = function(file) {
        pdf(file)
        print(gridExtra::grid.arrange(dend,
                                p,
                                nrow = 2,
                                top = grid::textGrob(
                                  paste0("Optimal number of supercluster is ", n.opt),
                                  gp = grid::gpar(fontsize = 20)
                                )))
        dev.off()
      }
    )
  }, 1000, 800)

  ###############################################################################
  ### Tree aggregation

  ### Plot tree with detailed scores
  output$TreePlot1 <- renderPlot({
    if (input$refTree) {
      reaRef()
      t1 <- reaValRef$treeObject
    } else{
      rea()
      t1 <- reaVal$treeObject
    }
    plot(t1)
    ape::edgelabels(t1$edge.length, font = 0.5)
    output$DownTreePlot1 <- downloadHandler(
      filename = function() {
        'TreePlot1.pdf'
      },
      content = function(file) {
        pdf(file)
        print(plot(t1))
        print(ape::edgelabels(t1$edge.length, font = 0.5))
        dev.off()
      }
    )
  }, 1000, 1200)

  ### Plot tree with selected layout
  output$TreePlot2 <- renderPlot({
    if (input$refTree) {
      reaRef()
      t1 <- reaValRef$treeObject
      newickCS <- reaValRef$newickCS
    } else{
      rea()
      t1 <- reaVal$treeObject
      newickCS <- reaVal$newickCS
    }
    TreeGroup <- newickCS@cluster
    nTreeColor <- length(TreeGroup)
    if(nTreeColor <= 8){
      TreeColor <- RColorBrewer::brewer.pal(nTreeColor,'Dark2')
    }else{
      TreeColor <- sample(grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,'Dark2'))(nTreeColor),nTreeColor)
    }
    i <- 1
    LabelColor <- c()
    EdgeColor <- rep("black", ape::Nedge(t1))
    for(i in 1:nTreeColor){
      nATree <- length(TreeGroup[[i]])
      AEdge <- ape::which.edge(t1, TreeGroup[[i]])
      EdgeColor[AEdge] <- TreeColor[i]
      LabelColor <- c(LabelColor,rep(TreeColor[i],nATree))
    }
    plot(t1,
         type = input$TreeLayout, # Layout
         use.edge.length = F,
         tip.color = LabelColor,
         edge.color = EdgeColor)
    hcp <- legend('bottomleft',
           title = 'SuperCluster',
           names(TreeGroup), col=TreeColor,
           lty=1, cex=1, lwd = 2,box.lty=0)
    hcp
    output$DownTreePlot2 <- downloadHandler(
      filename = function() {
        'TreePlot2.pdf'
      },
      content = function(file) {
        pdf(file)
        print(plot(t1,
             type = input$TreeLayout, # Layout
             use.edge.length = F,
             tip.color = LabelColor,
             edge.color = EdgeColor))
        print(hcp)
        dev.off()
      }
    )
    
  }, 1000, 1000)

  ### Plot subtree
  output$multiTreePlots <- renderUI({
    if (input$refTree) {
      reaRef()
      newickCS <- reaValRef$newickCS
    } else{
      rea()
      newickCS <- reaVal$newickCS
    }
    subtext <-
      newickCS@text[grep('\\(', newickCS@text)]
    getPlotList(subtext)
  })

  getPlotList <- function(subtext) {
    plotList <- lapply(1:length(subtext), function(x) {
      plotObject <- renderPlot({
        i <- subtext[[x]]
        i <- paste0('(', i, ');')
        t <- ape::read.tree(text = i)
        plot(t)
        ape::edgelabels(t$edge.length)
      })
    })
    tagList(plotList)
    return(plotList)
  }

  ###############################################################################
  ### Network partitioning
  output$GraphPlot <- renderPlot({
    if (input$refGraph) {
      reaGRef()
      GraphObject <- reaValGRef$GraphObject
    } else{
      reaG()
      GraphObject <- reaValG$GraphObject
    }
    
    tidyg <- tidygraph::as_tbl_graph(GraphObject)
    tidygf <-
      tidygraph::mutate(tidyg, community = as.factor(tidygraph::group_louvain(weights = weight)))
    
    tidyplot <- ggraph::ggraph(tidygf, layout = input$NetworkLayout) +
      ggraph::geom_node_point(ggplot2::aes(colour = community), size = 5) +
      ggraph::geom_edge_link(ggplot2::aes(width = weight), alpha = 0.2) +
      ggraph::scale_edge_width(range = c(0.2, 2)) +
      ggraph::geom_node_text(ggplot2::aes(label = label), repel = TRUE) +
      ggraph::theme_graph(base_family="sans")
    
    output$DownGraphPlot <- downloadHandler(
      filename = function() {
        'GraphPlot.pdf'
      },
      content = function(file) {
        pdf(file)
        print(tidyplot)
        dev.off()
      }
    )
    tidyplot

  }, 1000, 1000)

  ###############################################################################
  ### Upload cell embedding files
  observeEvent(input$embeddingsubmit, {
    req(input$umap)
    req(input$meta)
    tryCatch({
      apath <- input$umap$datapath
      umap <-
        read.csv(apath,
                 stringsAsFactors = FALSE,
                 row.names = 1)
      umap <- data.frame(lapply(umap, function(x) {
        gsub(" ", "_", x)
      }))
      write.csv(umap, paste0(workPath, 'umap.csv'))
      apath <- input$meta$datapath
      meta <-
        read.csv(apath,
                 stringsAsFactors = FALSE,
                 row.names = 1)
      meta <- data.frame(lapply(meta, function(x) {
        gsub(" ", "_", x)
      }))
      write.csv(meta, paste0(workPath, 'meta.csv'))
    },
    error = function(e) {
      stop(safeError(e))
    })
    updateSelectInput(session,
                      "color",
                      choices = colnames(meta))
    updateSelectInput(session,
                      "text",
                      choices = colnames(meta))
    updateSelectInput(session,
                      "indicluster",
                      choices = colnames(meta))
  })

  ###############################################################################
  ### Visualize cell in cell embedding plot
  observeEvent(input$metasubmit, {
    umap <-
      read.csv(
        paste0(workPath, 'umap.csv'),
        stringsAsFactors = FALSE,
        row.names = 1
      )
    meta <-
      read.csv(
        paste0(workPath, 'meta.csv'),
        stringsAsFactors = FALSE,
        row.names = 1
      )
    meta <- meta[rownames(umap),]
    umap1 <- cbind(umap, meta)

    output$UMAPPlotly <- plotly::renderPlotly({
      inputx <- colnames(umap)[1]
      inputy <- colnames(umap)[2]
      text <- ''
      for (atext in input$text) {
        text <- paste0(text, '"<br>', atext, ': ",', atext, ',')
      }
      text <-
        paste0('paste0("', substring(text, 6, nchar(text) - 1), ')')
      p <-
        ggplot2::ggplot(
          umap1,
          ggplot2::aes_string(
            x = inputx,
            y = inputy,
            color = input$color,
            text = text
          ),
          colors = colorRampPalette(RColorBrewer::brewer.pal(7, "Set1"))(100)
        ) + ggplot2::geom_point(size = 0.1, alpha = 0.5) + ggplot2::theme_bw()
      plotly::ggplotly(p)
    })
  })

  ### Visualize hierarchical clustering - global distances in cell embedding plot
  observeEvent(input$CSHierClustMetaSubmit, {
    withProgress(message = 'Create plot', value = 0.1, {
      umap <-
        read.csv(
          paste0(workPath, 'umap.csv'),
          stringsAsFactors = FALSE,
          row.names = 1
        )
      meta <-
        read.csv(
          paste0(workPath, 'meta.csv'),
          stringsAsFactors = FALSE,
          row.names = 1
        )
      meta <- meta[rownames(umap),]
      umap1 <- cbind(umap, meta)
      incProgress(0.1)
      indigrp <- meta[, input$indicluster]
      if (input$refUMAP) {
        p1 <- reaValScore$scoreMode2@p1
      } else{
        p1 <- reaValScore$scoreMode1@p1
      }
      hc <- p1$tree_row
      sub_grp.ori <- stats::cutree(hc, k = input$NBCutHclustU)
      sub_grp <- sub_grp.ori
      sub_grp <- paste0('Cluster', sub_grp)
      names(sub_grp) <- names(sub_grp.ori)
      coTreeGroup <- list()
      for (i in unique(sub_grp)) {
        coTreeGroup[[i]] <- names(sub_grp[sub_grp == i])
      }
      cs <- coTreeGroup
      CSCluster <- indigrp
      incProgress(0.2)
      for (i in 1:length(cs)) {
        acs <- cs[i]
        if (length(acs[[1]]) > 1) {
          for (aacsv in acs[[1]]) {
            CSCluster[CSCluster == aacsv] <- names(acs)
          }
        } else{
          CSCluster[CSCluster == acs[[1]]] <- acs[[1]]
        }
      }
      incProgress(0.1)
      grp <- CSCluster
      umap1$SuperCluster <- grp
      output$CSHierClustUMAPPlotly <- plotly::renderPlotly({
        inputx <- colnames(umap)[1]
        inputy <- colnames(umap)[2]
        text <- ''
        for (atext in input$text) {
          text <- paste0(text, '"<br>', atext, ': ",', atext, ',')
        }
        text <-
          paste0('paste0("', substring(text, 6, nchar(text) - 1), ')')
        p <-
          ggplot2::ggplot(
            umap1,
            ggplot2::aes_string(
              x = inputx,
              y = inputy,
              color = 'SuperCluster',
              text = text
            ),
            colors = colorRampPalette(RColorBrewer::brewer.pal(7, "Set1"))(100)
          ) + ggplot2::geom_point(size = 0.1, alpha = 0.5) + ggplot2::theme_bw()
        plotly::ggplotly(p)
      })
    })
  })

  ### Visualize hierarchical clustering - direct distances in cell embedding plot
  observeEvent(input$CSHierClustMetaSubmitD, {
    withProgress(message = 'Create plot', value = 0.1, {
      umap <-
        read.csv(
          paste0(workPath, 'umap.csv'),
          stringsAsFactors = FALSE,
          row.names = 1
        )
      meta <-
        read.csv(
          paste0(workPath, 'meta.csv'),
          stringsAsFactors = FALSE,
          row.names = 1
        )
      meta <- meta[rownames(umap),]
      umap1 <- cbind(umap, meta)
      incProgress(0.1)
      indigrp <- meta[, input$indicluster]
      if (input$refUMAP) {
        matrix <- reaValScore$scoreMode2@pair.score.matrix
      } else{
        matrix <- reaValScore$scoreMode1@pair.score.matrix
      }
      dis.direct <- apply(matrix, 1, function(x)
        (max(x) - x) / max(x))
      dis.direct <- t(dis.direct)
      dis.direct <- stats::as.dist(dis.direct)
      hc <- stats::hclust(dis.direct)
      sub_grp.ori <- stats::cutree(hc, k = input$NBCutHclustUD)
      sub_grp <- sub_grp.ori
      sub_grp <- paste0('Cluster', sub_grp)
      names(sub_grp) <- names(sub_grp.ori)
      coTreeGroup <- list()
      for (i in unique(sub_grp)) {
        coTreeGroup[[i]] <- names(sub_grp[sub_grp == i])
      }
      cs <- coTreeGroup
      CSCluster <- indigrp
      incProgress(0.2)
      for (i in 1:length(cs)) {
        acs <- cs[i]
        if (length(acs[[1]]) > 1) {
          for (aacsv in acs[[1]]) {
            CSCluster[CSCluster == aacsv] <- names(acs)
          }
        } else{
          CSCluster[CSCluster == acs[[1]]] <- acs[[1]]
        }
      }
      incProgress(0.1)
      grp <- CSCluster
      umap1$SuperCluster <- grp
      output$CSHierClustUMAPPlotlyD <- plotly::renderPlotly({
        inputx <- colnames(umap)[1]
        inputy <- colnames(umap)[2]
        text <- ''
        for (atext in input$text) {
          text <- paste0(text, '"<br>', atext, ': ",', atext, ',')
        }
        text <-
          paste0('paste0("', substring(text, 6, nchar(text) - 1), ')')
        p <-
          ggplot2::ggplot(
            umap1,
            ggplot2::aes_string(
              x = inputx,
              y = inputy,
              color = 'SuperCluster',
              text = text
            ),
            colors = colorRampPalette(RColorBrewer::brewer.pal(7, "Set1"))(100)
          ) + ggplot2::geom_point(size = 0.1, alpha = 0.5) + ggplot2::theme_bw()
        plotly::ggplotly(p)
      })
    })
  })

  ### Visualize tree aggregation in cell embedding plot
  observeEvent(input$CSTreeMetaSubmit, {
    withProgress(message = 'Create plot', value = 0.1, {
      umap <-
        read.csv(
          paste0(workPath, 'umap.csv'),
          stringsAsFactors = FALSE,
          row.names = 1
        )
      meta <-
        read.csv(
          paste0(workPath, 'meta.csv'),
          stringsAsFactors = FALSE,
          row.names = 1
        )
      meta <- meta[rownames(umap),]
      umap1 <- cbind(umap, meta)
      incProgress(0.1)
      indigrp <- meta[, input$indicluster]

      if (input$refUMAP) {
        reaURef()
        newickCS <- reaValURef$newickCS
      } else{
        reaU()
        newickCS <- reaValU$newickCS
      }

      cs <- newickCS@cluster
      CSCluster <- indigrp
      incProgress(0.2)
      for (i in 1:length(cs)) {
        acs <- cs[i]
        if (length(acs[[1]]) > 1) {
          for (aacsv in acs[[1]]) {
            CSCluster[CSCluster == aacsv] <- names(acs)
          }
        } else{
          CSCluster[CSCluster == acs[[1]]] <- acs[[1]]
        }
      }

      incProgress(0.1)
      grp <- CSCluster
      umap1$SuperCluster <- grp
      output$CSTreeUMAPPlotly <- plotly::renderPlotly({
        inputx <- colnames(umap)[1]
        inputy <- colnames(umap)[2]
        text <- ''
        for (atext in input$text) {
          text <- paste0(text, '"<br>', atext, ': ",', atext, ',')
        }
        text <-
          paste0('paste0("', substring(text, 6, nchar(text) - 1), ')')
        p <-
          ggplot2::ggplot(
            umap1,
            ggplot2::aes_string(
              x = inputx,
              y = inputy,
              color = 'SuperCluster',
              text = text
            ),
            colors = colorRampPalette(RColorBrewer::brewer.pal(7, "Set1"))(100)
          ) + ggplot2::geom_point(size = 0.1, alpha = 0.5) + ggplot2::theme_bw()
        plotly::ggplotly(p)
      })
    })
  })

  ### Visualize network partitioning in cell embedding plot
  observeEvent(input$CSGraphMetaSubmit, {
    withProgress(message = 'Create plot', value = 0.1, {
      umap <-
        read.csv(
          paste0(workPath, 'umap.csv'),
          stringsAsFactors = FALSE,
          row.names = 1
        )
      meta <-
        read.csv(
          paste0(workPath, 'meta.csv'),
          stringsAsFactors = FALSE,
          row.names = 1
        )
      meta <- meta[rownames(umap),]
      umap1 <- cbind(umap, meta)
      incProgress(0.1)
      indigrp <- meta[, input$indicluster]

      if (input$refUMAP) {
        reaGURef()
        newickCS <- reaValGURef$GraphCS
      } else{
        reaGU()
        newickCS <- reaValGU$GraphCS
      }
      cs <- newickCS@cluster
      CSCluster <- indigrp
      incProgress(0.2)
      for (i in 1:length(cs)) {
        acs <- cs[i]
        if (length(acs[[1]]) > 1) {
          for (aacsv in acs[[1]]) {
            CSCluster[CSCluster == aacsv] <- names(acs)
          }
        } else{
          CSCluster[CSCluster == acs[[1]]] <- acs[[1]]
        }
      }
      incProgress(0.1)
      grp <- CSCluster
      umap1$SuperCluster <- grp
      output$CSGraphUMAPPlotly <- plotly::renderPlotly({
        inputx <- colnames(umap)[1]
        inputy <- colnames(umap)[2]
        text <- ''
        for (atext in input$text) {
          text <- paste0(text, '"<br>', atext, ': ",', atext, ',')
        }
        text <-
          paste0('paste0("', substring(text, 6, nchar(text) - 1), ')')
        p <-
          ggplot2::ggplot(
            umap1,
            ggplot2::aes_string(
              x = inputx,
              y = inputy,
              color = 'SuperCluster',
              text = text
            ),
            colors = colorRampPalette(RColorBrewer::brewer.pal(7, "Set1"))(100)
          ) + ggplot2::geom_point(size = 0.1, alpha = 0.5) + ggplot2::theme_bw()
        plotly::ggplotly(p)
      })
    })
  })



  ###############################################################################
  ### Print overlapping genes - Hierarchical clustering - global distances
  output$MGeneList <- renderPrint({
    if (input$HeatmapMode == 1) {
      p1 <- reaValScore$scoreMode1@p1
      matrix <- reaValScore$scoreMode1@pair.score.matrix
      marker.table <- reaValScore$scoreMode1@marker.table
    }
    if (input$HeatmapMode == 2) {
      p1 <- reaValScore$scoreMode2@p1
      matrix <- reaValScore$scoreMode2@pair.score.matrix
      marker.table <- reaValScore$scoreMode2@marker.table
    }
    if (input$HeatmapMode == 3) {
      p1 <- reaValScore$scoreMode3@p1
      matrix <- reaValScore$scoreMode3@pair.score.matrix
      marker.table <- reaValScore$scoreMode3@marker.table
    }
    
    hc <- p1$tree_row
    sub_grp.ori <- stats::cutree(hc, k = input$NBCutHclust)
    sub_grp <- sub_grp.ori
    sub_grp <- paste0('Cluster', sub_grp)
    names(sub_grp) <- names(sub_grp.ori)
    coTreeGroup <- list()
    for (i in unique(sub_grp)) {
      coTreeGroup[[i]] <- names(sub_grp[sub_grp == i])
    }
    coTreeGene <- findOverlapGene(marker.table, coTreeGroup)
    print(coTreeGene)
  })

  ### Print supercluster composition - Hierarchical clustering - global distances
  output$MClusterList <- renderPrint({
    if (input$HeatmapMode == 1) {
      p1 <- reaValScore$scoreMode1@p1
    }
    if (input$HeatmapMode == 2) {
      p1 <- reaValScore$scoreMode2@p1
    }
    if (input$HeatmapMode == 3) {
      p1 <- reaValScore$scoreMode3@p1
    }
    hc <- p1$tree_row
    sub_grp.ori <- stats::cutree(hc, k = input$NBCutHclust)
    sub_grp <- sub_grp.ori
    sub_grp <- paste0('Cluster', sub_grp)
    names(sub_grp) <- names(sub_grp.ori)
    coTreeGroup <- list()
    for (i in unique(sub_grp)) {
      coTreeGroup[[i]] <- names(sub_grp[sub_grp == i])
    }
    print(coTreeGroup)
  })

  ### Download table of overlapping genes - Hierarchical clustering - global distances
  output$MGeneListT <- DT::renderDataTable({
    if (input$HeatmapMode == 1) {
      p1 <- reaValScore$scoreMode1@p1
      matrix <- reaValScore$scoreMode1@pair.score.matrix
      marker.table <- reaValScore$scoreMode1@marker.table
    }
    if (input$HeatmapMode == 2) {
      p1 <- reaValScore$scoreMode2@p1
      matrix <- reaValScore$scoreMode2@pair.score.matrix
      marker.table <- reaValScore$scoreMode2@marker.table
    }
    if (input$HeatmapMode == 3) {
      p1 <- reaValScore$scoreMode3@p1
      matrix <- reaValScore$scoreMode3@pair.score.matrix
      marker.table <- reaValScore$scoreMode3@marker.table
    }
    hc <- p1$tree_row
    sub_grp.ori <- stats::cutree(hc, k = input$NBCutHclust)
    sub_grp <- sub_grp.ori
    sub_grp <- paste0('Cluster', sub_grp)
    names(sub_grp) <- names(sub_grp.ori)
    coTreeGroup <- list()
    for (i in unique(sub_grp)) {
      coTreeGroup[[i]] <- names(sub_grp[sub_grp == i])
    }
    coTreeGene <- findOverlapGene(marker.table, coTreeGroup)

    df <- data.frame(Reference = character(), Gene = character())
    for (i in 1:length(coTreeGene)) {
      gl <- coTreeGene[i]
      adf <- data.frame(SuperCluster = names(gl), Gene = gl[[1]])
      df <- rbind(df, adf)
    }
    output$MDownGeneListT <- downloadHandler(
      filename = function() {
        'ModeGeneTable.csv'
      },
      content = function(file) {
        write.csv(df, file)
      }
    )
    return(df)
  })

  ### Download table of supercluster compostion - Hierarchical clustering - global distances
  output$MClusterListT <- DT::renderDataTable({
    if (input$HeatmapMode == 1) {
      p1 <- reaValScore$scoreMode1@p1
    }
    if (input$HeatmapMode == 2) {
      p1 <- reaValScore$scoreMode2@p1
    }
    if (input$HeatmapMode == 3) {
      p1 <- reaValScore$scoreMode3@p1
    }
    hc <- p1$tree_row
    sub_grp.ori <- stats::cutree(hc, k = input$NBCutHclust)
    sub_grp <- sub_grp.ori
    sub_grp <- paste0('Cluster', sub_grp)
    names(sub_grp) <- names(sub_grp.ori)
    coTreeGroup <- list()
    for (i in unique(sub_grp)) {
      coTreeGroup[[i]] <- names(sub_grp[sub_grp == i])
    }

    df <-
      data.frame(SuperCluster = character(), Cluster = character())
    for (i in 1:length(coTreeGroup)) {
      gl <- coTreeGroup[i]
      adf <- data.frame(SuperCluster = names(gl), Cluster = gl[[1]])
      df <- rbind(df, adf)
    }

    output$MDownClusterListT <- downloadHandler(
      filename = function() {
        'ModeClusterTable.csv'
      },
      content = function(file) {
        write.csv(df, file)
      }
    )
    return(df)
  })

  ### Print overlapping genes - Hierarchical clustering - direct distances
  output$MGeneListD <- renderPrint({
    if (input$HeatmapModeD == 1) {
      matrix <- reaValScore$scoreMode1@pair.score.matrix
      marker.table <- reaValScore$scoreMode1@marker.table
    }
    if (input$HeatmapModeD == 2) {
      matrix <- reaValScore$scoreMode2@pair.score.matrix
      marker.table <- reaValScore$scoreMode2@marker.table
    }
    if (input$HeatmapModeD == 3) {
      matrix <- reaValScore$scoreMode3@pair.score.matrix
      matrix <- t(matrix)
      marker.table <- reaValScore$scoreMode3@marker.table
    }
    dis.direct <- apply(matrix, 1, function(x)
      (max(x) - x) / max(x))
    dis.direct <- t(dis.direct)
    dis.direct <- stats::as.dist(dis.direct)
    hc <- hclust(dis.direct)
    sub_grp.ori <- stats::cutree(hc, k = input$NBCutHclustD)
    sub_grp <- sub_grp.ori
    sub_grp <- paste0('Cluster', sub_grp)
    names(sub_grp) <- names(sub_grp.ori)
    coTreeGroup <- list()
    for (i in unique(sub_grp)) {
      coTreeGroup[[i]] <- names(sub_grp[sub_grp == i])
    }
    coTreeGene <- findOverlapGene(marker.table, coTreeGroup)
    print(coTreeGene)
  })

  ### Print supercluster composition - Hierarchical clustering - direct distances
  output$MClusterListD <- renderPrint({
    if (input$HeatmapModeD == 1) {
      matrix <- reaValScore$scoreMode1@pair.score.matrix
    }
    if (input$HeatmapModeD == 2) {
      matrix <- reaValScore$scoreMode2@pair.score.matrix
    }
    if (input$HeatmapModeD == 3) {
      matrix <- reaValScore$scoreMode3@pair.score.matrix
      matrix <- t(matrix)
    }
    dis.direct <- apply(matrix, 1, function(x)
      (max(x) - x) / max(x))
    dis.direct <- t(dis.direct)
    dis.direct <- stats::as.dist(dis.direct)
    hc <- stats::hclust(dis.direct)
    sub_grp.ori <- stats::cutree(hc, k = input$NBCutHclustD)
    sub_grp <- sub_grp.ori
    sub_grp <- paste0('Cluster', sub_grp)
    names(sub_grp) <- names(sub_grp.ori)
    coTreeGroup <- list()
    for (i in unique(sub_grp)) {
      coTreeGroup[[i]] <- names(sub_grp[sub_grp == i])
    }
    print(coTreeGroup)
  })

  ### Download overlapping genes - Hierarchical clustering - direct distances
  output$MGeneListTD <- DT::renderDataTable({
    if (input$HeatmapModeD == 1) {
      matrix <- reaValScore$scoreMode1@pair.score.matrix
      marker.table <- reaValScore$scoreMode1@marker.table
    }
    if (input$HeatmapModeD == 2) {
      matrix <- reaValScore$scoreMode2@pair.score.matrix
      marker.table <- reaValScore$scoreMode2@marker.table
    }
    if (input$HeatmapModeD == 3) {
      matrix <- reaValScore$scoreMode3@pair.score.matrix
      matrix <- t(matrix)
      marker.table <- reaValScore$scoreMode3@marker.table
    }
    dis.direct <- apply(matrix, 1, function(x)
      (max(x) - x) / max(x))
    dis.direct <- t(dis.direct)
    dis.direct <- stats::as.dist(dis.direct)
    hc <- stats::hclust(dis.direct)
    sub_grp.ori <- stats::cutree(hc, k = input$NBCutHclustD)
    sub_grp <- sub_grp.ori
    sub_grp <- paste0('Cluster', sub_grp)
    names(sub_grp) <- names(sub_grp.ori)
    coTreeGroup <- list()
    for (i in unique(sub_grp)) {
      coTreeGroup[[i]] <- names(sub_grp[sub_grp == i])
    }
    coTreeGene <- findOverlapGene(marker.table, coTreeGroup)
    df <- data.frame(Reference = character(), Gene = character())
    for (i in 1:length(coTreeGene)) {
      gl <- coTreeGene[i]
      adf <- data.frame(SuperCluster = names(gl), Gene = gl[[1]])
      df <- rbind(df, adf)
    }
    output$MDownGeneListTD <- downloadHandler(
      filename = function() {
        'ModeGeneTableDirect.csv'
      },
      content = function(file) {
        write.csv(df, file)
      }
    )
    return(df)
  })

  ### Download supercluster composition  - Hierarchical clustering - direct distances
  output$MClusterListTD <- DT::renderDataTable({
    if (input$HeatmapModeD == 1) {
      matrix <- reaValScore$scoreMode1@pair.score.matrix
    }
    if (input$HeatmapModeD == 2) {
      matrix <- reaValScore$scoreMode2@pair.score.matrix
    }
    if (input$HeatmapModeD == 3) {
      matrix <- reaValScore$scoreMode3@pair.score.matrix
      matrix <- t(matrix)
    }
    dis.direct <- apply(matrix, 1, function(x)
      (max(x) - x) / max(x))
    dis.direct <- t(dis.direct)
    dis.direct <- stats::as.dist(dis.direct)
    hc <- stats::hclust(dis.direct)
    sub_grp.ori <- stats::cutree(hc, k = input$NBCutHclustD)
    sub_grp <- sub_grp.ori
    sub_grp <- paste0('Cluster', sub_grp)
    names(sub_grp) <- names(sub_grp.ori)
    coTreeGroup <- list()
    for (i in unique(sub_grp)) {
      coTreeGroup[[i]] <- names(sub_grp[sub_grp == i])
    }

    df <-
      data.frame(SuperCluster = character(), Cluster = character())
    for (i in 1:length(coTreeGroup)) {
      gl <- coTreeGroup[i]
      adf <- data.frame(SuperCluster = names(gl), Cluster = gl[[1]])
      df <- rbind(df, adf)
    }

    output$MDownClusterListTD <- downloadHandler(
      filename = function() {
        'ModeClusterTableDirect.csv'
      },
      content = function(file) {
        write.csv(df, file)
      }
    )
    return(df)
  })


  ### Print overlapping genes - tree aggregation
  output$GeneList <- renderPrint({
    if (input$refTree) {
      reaRef()
      newickCS <- reaValRef$newickCS
    } else{
      rea()
      newickCS <- reaVal$newickCS
    }
    print(newickCS@gene)
  })

  ### Print supercluster composition - tree aggregation
  output$ClusterList <- renderPrint({
    if (input$refTree) {
      reaRef()
      newickCS <- reaValRef$newickCS
    } else{
      rea()
      newickCS <- reaVal$newickCS
    }
    print(newickCS@cluster)
  })

  ### Print overlapping genes - network partitioning
  output$GGeneList <- renderPrint({
    if (input$refGraph) {
      reaGRef()
      GraphCS <- reaValGRef$GraphCS
    } else{
      reaG()
      GraphCS <- reaValG$GraphCS
    }
    print(GraphCS@gene)
  })

  ### Print supercluster composition - network partitioning
  output$GClusterList <- renderPrint({
    if (input$refGraph) {
      reaGRef()
      GraphCS <- reaValGRef$GraphCS
    } else{
      reaG()
      GraphCS <- reaValG$GraphCS
    }
    print(GraphCS@cluster)
  })

  ### Download table of overlapping genes - tree aggregation
  output$GeneListT <- DT::renderDataTable({
    if (input$refTree) {
      reaRef()
      newickCS <- reaValRef$newickCS
    } else{
      rea()
      newickCS <- reaVal$newickCS
    }
    df <- data.frame(SuperCluster = character(), Gene = character())
    for (i in 1:length(newickCS@gene)) {
      gl <- newickCS@gene[i]
      adf <- data.frame(SuperCluster = names(gl), Gene = gl[[1]])
      df <- rbind(df, adf)
    }
    output$DownGeneListT <- downloadHandler(
      filename = function() {
        'CSTreeGeneTable.csv'
      },
      content = function(file) {
        write.csv(df, file)
      }
    )
    return(df)
  })

  ### Download table of supercluster composition - tree aggregation
  output$ClusterListT <- DT::renderDataTable({
    if (input$refTree) {
      reaRef()
      newickCS <- reaValRef$newickCS
    } else{
      rea()
      newickCS <- reaVal$newickCS
    }
    df <-
      data.frame(SuperCluster = character(), Cluster = character())
    for (i in 1:length(newickCS@cluster)) {
      gl <- newickCS@cluster[i]
      adf <- data.frame(SuperCluster = names(gl), Cluster = gl[[1]])
      df <- rbind(df, adf)
    }
    output$DownClusterListT <- downloadHandler(
      filename = function() {
        'CSTreeClusterTable.csv'
      },
      content = function(file) {
        write.csv(df, file)
      }
    )
    return(df)
  })

  ### Download table of overlapping genes - network partitioning
  output$GGeneListT <- DT::renderDataTable({
    if (input$refGraph) {
      reaGRef()
      GraphCS <- reaValGRef$GraphCS
    } else{
      reaG()
      GraphCS <- reaValG$GraphCS
    }
    df <- data.frame(SuperCluster = character(), Gene = character())
    for (i in 1:length(GraphCS@gene)) {
      gl <- GraphCS@gene[i]
      adf <- data.frame(SuperCluster = names(gl), Gene = gl[[1]])
      df <- rbind(df, adf)
    }
    output$DownGGeneListT <- downloadHandler(
      filename = function() {
        'CSGraphGeneTable.csv'
      },
      content = function(file) {
        write.csv(df, file)
      }
    )
    return(df)
  })

  ### Download table of supercluster composition - network partitioning
  output$GClusterListT <- DT::renderDataTable({
    if (input$refGraph) {
      reaGRef()
      GraphCS <- reaValGRef$GraphCS
    } else{
      reaG()
      GraphCS <- reaValG$GraphCS
    }
    df <-
      data.frame(SuperCluster = character(), Cluster = character())
    for (i in 1:length(GraphCS@cluster)) {
      gl <- GraphCS@cluster[i]
      adf <- data.frame(SuperCluster = names(gl), Cluster = gl[[1]])
      df <- rbind(df, adf)
    }
    output$DownGClusterListT <- downloadHandler(
      filename = function() {
        'CSGraphClusterTable.csv'
      },
      content = function(file) {
        write.csv(df, file)
      }
    )
    return(df)
  })


  ###############################################################################
  ### Remove uploaded files
  session$onSessionEnded(function() {
    system(paste("rm -r", workPath))
    rm(list = ls())
  })
}
