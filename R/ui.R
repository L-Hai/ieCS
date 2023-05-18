#' ShinyApp UI
#'
#' @import shiny
ui <- navbarPage(
  "ieCS",
  tabPanel("UploadData",
           sidebarLayout(
             sidebarPanel(
               fileInput(
                 "MarkersFile",
                 "Choose Individual Cluster Markers CSV File",
                 multiple = FALSE,
                 accept = c("text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv")
               ),
               fileInput(
                 "RefMarkersFile",
                 "Choose Reference Cell Type Markers CSV File (Optional)",
                 multiple = FALSE,
                 accept = c("text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv")
               ),
               actionButton('MarkersSubmit', ' Submit ')
             ),
             mainPanel(
               column(
                 6,
                 h5("Individual Cluster Markers"),
                 selectInput('OrderM', 'Markers Order by', ''),
                 radioButtons('DirectionM', 'The Marker With Bigger Value is', 
                              list('More Important' = TRUE, 'Less Important' = FALSE), TRUE),
                 selectInput('GeneM', 'Markers', ''),
                 selectInput('ClusterM', 'Cell Cluster Information', '')
               ),
               column(
                 6,
                 h5("Reference Cell Type Markers"),
                 selectInput('OrderMR', 'Markers Order by', ''),
                 radioButtons('DirectionMR', 'The Marker With Bigger Value is', 
                              list('More Important' = TRUE, 'Less Important' = FALSE), TRUE),
                 selectInput('GeneMR', 'Markers', ''),
                 selectInput('ClusterMR', 'Cell Type Information', '')
               ),
               hr(),
               actionButton('RunScore', ' Run '),
               hr(),
               DT::dataTableOutput("contents")
             )
           )),
  tabPanel(
    "CSHierClust",
    radioButtons(
      "HeatmapMode",
      "Show result from:",
      c(
        "Query x Query" = 1,
        "(Query + Reference) x (Query + Reference)" = 2,
        "Query x Reference" = 3
      ),
      selected = 1
    ),
    tabsetPanel(
      id = 'GeneralModeTab',
      tabPanel(
        'Heatmap',
        downloadButton('MDownMatrixT', 'Download: CS Matrix'),
        downloadButton('DownHeatmap', 'Download: Heatmap'),
        plotOutput("Heatmap", 1000, 1000)
      ),
      tabPanel(
        'HierarchicalClustering_Global',
        tabsetPanel(
          tabPanel('Dendrogram',
                   downloadButton('DownHclust', 'Download'),
                   plotOutput("Hclust", 1000, 500)),
          tabPanel(
            'OptimalClusterAssignment',
            downloadButton('DownOptHclust', 'Download'),
            plotOutput("OptHclust", 1000, 1500)
          ),
          tabPanel(
            'CustomClusterAssignment',
            numericInput("NBCutHclust", "NumberOfSuperCluster:", 5, min = 2),
            downloadButton('DownCutHclust', 'Download'),
            plotOutput("CutHclust",1000,500)
          ),
          tabPanel("Supercluster", tabsetPanel(
            tabPanel('Text', verbatimTextOutput("MClusterList")),
            tabPanel(
              'Table',
              downloadButton('MDownClusterListT', 'Download'),
              DT::dataTableOutput("MClusterListT")
            )
          )),
          tabPanel("OverlapGene", tabsetPanel(
            tabPanel('Text', verbatimTextOutput("MGeneList")),
            tabPanel(
              'Table',
              downloadButton('MDownGeneListT', 'Download'),
              DT::dataTableOutput("MGeneListT")
            )
          ))
        )
      ),
      tabPanel(
        'Detail: Query x Reference',
        tabsetPanel(
          tabPanel('Dendrogram',
                   downloadButton('DownHclustMode3', 'Download'),
                   plotOutput("HclustMode3", 1000, 500)),
          tabPanel(
            'OptimalClusterAssignment',
            collapsibleTree::collapsibleTreeOutput("CollectMode3", height = "500px")
          ),
          tabPanel("Supercluster", tabsetPanel(
            tabPanel('Text', verbatimTextOutput("M3ClusterList")),
            tabPanel(
              'Table',
              downloadButton('M3DownClusterListT', 'Download'),
              DT::dataTableOutput("M3ClusterListT")
            )
          )),
          tabPanel("OverlapGene", tabsetPanel(
            tabPanel('Text', verbatimTextOutput("M3GeneList")),
            tabPanel(
              'Table',
              downloadButton('M3DownGeneListT', 'Download'),
              DT::dataTableOutput("M3GeneListT")
            )
          ))
        )
      )
    )
  ),
  tabPanel(
    "CSHierClustDirect",
    radioButtons(
      "HeatmapModeD",
      "Show result from:",
      c(
        "Query x Query" = 1,
        "(Query + Reference) x (Query + Reference)" = 2
      ),
      selected = 1
    ),
    tabsetPanel(id = 'GeneralModeTabD',
                tabPanel(
                  'HierarchicalClustering_Direct',
                  tabsetPanel(
                    tabPanel('Dendrogram',
                             downloadButton('DownHclustD', 'Download'),
                             plotOutput("HclustD", 1000, 500)),
                    tabPanel(
                      'OptimalClusterAssignment',
                      downloadButton('DownOptHclustD', 'Download'),
                      plotOutput("OptHclustD", 1000, 800)
                    ),
                    tabPanel(
                      'CustomClusterAssignment',
                      numericInput("NBCutHclustD", "NumberOfSuperCluster:", 5, min = 2),
                      downloadButton('DownCutHclustD', 'Download'),
                      plotOutput("CutHclustD", 1000, 500)
                    ),
                    tabPanel("Supercluster", tabsetPanel(
                      tabPanel('Text', verbatimTextOutput("MClusterListD")),
                      tabPanel(
                        'Table',
                        downloadButton('MDownClusterListTD', 'Download'),
                        DT::dataTableOutput("MClusterListTD")
                      )
                    )),
                    tabPanel("OverlapGene", tabsetPanel(
                      tabPanel('Text', verbatimTextOutput("MGeneListD")),
                      tabPanel(
                        'Table',
                        downloadButton('MDownGeneListTD', 'Download'),
                        DT::dataTableOutput("MGeneListTD")
                      )
                    ))
                  )
                ))
  ),
  tabPanel("CSTree",
           sidebarLayout(
             sidebarPanel(numericInput("minScore", "Minimum Score:", 5, min = 1)),
             mainPanel(
               radioButtons(
                 "refTree",
                 "Show result from:",
                 c("Individual Clusters" = FALSE,
                   "With Reference" = TRUE),
                 selected = FALSE
               ),
               tabsetPanel(
                 tabPanel("TreePlot", tabsetPanel(
                   tabPanel('Scored', 
                            downloadButton('DownTreePlot1', 'Download'),
                            plotOutput("TreePlot1", 1000, 1200)),
                   tabPanel('Layout', 
                            selectInput('TreeLayout', 'Change Layout', 
                                        c("fan", "radial", 'cladogram', "phylogram", "unrooted"),
                                       'fan'),
                            downloadButton('DownTreePlot2', 'Download'),
                            plotOutput("TreePlot2", 1000, 1000))
                 )),
                 tabPanel("Supercluster", tabsetPanel(
                   tabPanel('Text', verbatimTextOutput("ClusterList")),
                   tabPanel(
                     'Table',
                     downloadButton('DownClusterListT', 'Download'),
                     DT::dataTableOutput("ClusterListT")
                   )
                 )),
                 tabPanel("OverlapGene", tabsetPanel(
                   tabPanel('Text', verbatimTextOutput("GeneList")),
                   tabPanel(
                     'Table',
                     downloadButton('DownGeneListT', 'Download'),
                     DT::dataTableOutput("GeneListT")
                   )
                 )),
                 tabPanel("SubTreePlot", uiOutput("multiTreePlots"))
               )
             )
           )),
  tabPanel("CSNetwork",
           sidebarLayout(
             sidebarPanel(numericInput(
               "minScoreG",
               "Minimum Score:",
               5,
               min = 1,
               max = 112
             )),
             mainPanel(
               radioButtons(
                 "refGraph",
                 "Show result from:",
                 c("Individual Clusters" = FALSE,
                   "With Reference" = TRUE),
                 selected = FALSE
               ),
               tabsetPanel(
                 tabPanel("NetworkPlot",
                          selectInput('NetworkLayout', 'Change Layout', 
                                      c('auto','fr','kk','gem', 'dh', 'graphopt','mds', 'drl', 'lgl'),
                                      'auto'),
                          downloadButton('DownGraphPlot', 'Download'),
                          plotOutput("GraphPlot")),
                 tabPanel("Supercluster", tabsetPanel(
                   tabPanel('Text', verbatimTextOutput("GClusterList")),
                   tabPanel(
                     'Table',
                     downloadButton('DownGClusterListT', 'Download'),
                     DT::dataTableOutput("GClusterListT")
                   )
                 )),
                 tabPanel("OverlapGene", tabsetPanel(
                   tabPanel('Text', verbatimTextOutput("GGeneList")),
                   tabPanel(
                     'Table',
                     downloadButton('DownGGeneListT', 'Download'),
                     DT::dataTableOutput("GGeneListT")
                   )
                 ))
               )
             )
           )),
  tabPanel(
    "CellEmbeddingPlot",
    sidebarLayout(
      sidebarPanel(
        fileInput(
          "umap",
          "Choose Cell Embedding CSV File",
          multiple = FALSE,
          accept = c("text/csv",
                     "text/comma-separated-values,text/plain",
                     ".csv")

        ),
        fileInput(
          "meta",
          "Choose Cell Annotation CSV File",
          multiple = FALSE,
          accept = c("text/csv",
                     "text/comma-separated-values,text/plain",
                     ".csv")
        ),
        actionButton('embeddingsubmit', ' Submit ')
      ),

      mainPanel(
        selectInput('indicluster', 'Individual Cluster Information', ''),
        selectInput('text', 'Hover Annotation', '', multiple = TRUE),
        radioButtons(
          "refUMAP",
          "Show result from:",
          c("Individual Clusters" = FALSE,
            "With Reference" = TRUE),
          selected = FALSE
        ),
        tabsetPanel(
          tabPanel(
            "Meta",
            selectInput('color', 'Color by', ''),
            actionButton('metasubmit', ' Plot '),
            plotly::plotlyOutput("UMAPPlotly")
          ),
          tabPanel(
            "CSHierClust",
            numericInput("NBCutHclustU", "NumberOfSuperCluster:", 5, min = 1),
            actionButton('CSHierClustMetaSubmit', ' Plot '),
            plotly::plotlyOutput("CSHierClustUMAPPlotly")
          ),
          tabPanel(
            "CSHierClustDirect",
            numericInput("NBCutHclustUD", "NumberOfSuperCluster:", 5, min = 1),
            actionButton('CSHierClustMetaSubmitD', ' Plot '),
            plotly::plotlyOutput("CSHierClustUMAPPlotlyD")
          ),
          tabPanel(
            "CSTree",
            numericInput("minScoreU", "Minimum Score:", 5, min = 1),
            actionButton('CSTreeMetaSubmit', ' Plot '),
            plotly::plotlyOutput("CSTreeUMAPPlotly")
          ),
          tabPanel(
            "CSNetwork",
            numericInput("minScoreGU", "Minimum Score:", 5, min = 1),
            actionButton('CSGraphMetaSubmit', ' Plot '),
            plotly::plotlyOutput("CSGraphUMAPPlotly")
          )
        )
      )
    )
  )
)
