############################
# MeSHDD Shiny Application #
# Written by Adam Brown    #
# Updated 2/22/2016        #
############################


## Load ##
library(shiny)
require(DT)
require(stats)
require(dendextend)
load('data/shiny.RData')

## UI Definition ##
ui <- fluidPage(
  # Header
  titlePanel("MeSHDD: MeSH-based Drug-Drug Similarity and Repositioning"),
  p(
    'MeSHDD uses MeSH-term enrichment to discover literature-based similarities between FDA approved drugs. 
    Drug-drug similarities can then be used to robustly cluster drugs and discover clusters thare are overrepresented
    for the treatment of a given disease. MeSHDD can be used in two ways:'),
  tags$ol(
    tags$li('In a drug-centric way, guiding your search using a drug you would like to reposition'),
    tags$li('In a disease centric way, guiding your search by selecting first a disease of interest and
       then a cluster containing drugs that might treat that disease')
  ),
  
  p('MeSHDD is published in TBA, and source code
      will be made available upon publication on', a("GitHub",href='http://github.com/adam-sam-brown/')),
  p('Visit the', a('Patel Group Homepage', href='http://www.chiragjpgroup.org/'), 'for other projects.'),
  p('Application by Adam Brown. This work is licensed under a', 
    a('Creative Commons Attribution-NonCommercial 4.0 International License',
      href="http://creativecommons.org/licenses/by-nc/4.0/")),
  
  # Big tabset
  tabsetPanel(
    # Drug-centric
    tabPanel(
      title = 'Drug-Centric MeSHDD',
      sidebarLayout(
        sidebarPanel(
          p('One way to use MeSHDD is to search for a drug you would like to reposition.
            Start by selecting a drug in the panel below and customize your output as desired'),
          selectInput(
            inputId = 'drugName', 
            label = 'Choose a drug from the dropdown menu:',
            choices = drugNames,
            selected = 'Metformin'
          )
          
        ),
        
        mainPanel(
          tabsetPanel(
            tabPanel(paste('Cluster Dendrogram'), 
                     plotOutput(outputId = 'cladogram'),
                     dataTableOutput(outputId = 'indication')),
            tabPanel('Similar Drugs', dataTableOutput(outputId = 'simTable')),
            tabPanel('Enriched MeSH Terms for Selected Drug', dataTableOutput(outputId = 'meshTable'))
          )
        )
      )
    ),
    
    tabPanel(
      title = 'Disease-Centric MeSHDD',
      sidebarLayout(
        sidebarPanel(
          p('Another way to use MeshDD is to search for clusters that contain drugs overrepresented for treatment of a specific disease
        Start by selecting a disease in the dropdown menu below.'),
          uiOutput('diseaseSelect'),
          p('Then, select a cluster that is enriched for treatment of that disease:'),
          uiOutput('clusterSelect')
          
        ),
        mainPanel(
          tabsetPanel(
            tabPanel(title = 'Cluster Cladogram',
                     plotOutput(outputId = 'clustergram'))#,
            #         tabPanel(title = 'Drug List',
            #                  dataTableOutput('drugTable'))
          )
        )
      )
    )
  )
)

## Server Definition ##
server <- function(input, output) {
  # MeSH term display
  output$meshTable <- renderDataTable({
    DT::datatable(enrList[[input$drugName]], options = list(pageLength = 10))
  })
  
  # Similar drugs display
  output$simTable <- renderDataTable({
    DT::datatable(distList[[input$drugName]][-1,], options = list(pageLength = 10))
  })
  
  # Cladogram plotting
  output$cladogram <- renderPlot({
    clustNum <- unname(partition[input$drugName])
    dend <- as.dendrogram(hclustList[[clustNum]])
    labels_colors(dend) <- ifelse(labels(dend) == input$drugName, 'red', 'black')
    plot(
      dend,
      main=paste('Plotting Cluster Containing', input$drugName), 
      xlab=NA,
      sub=NA,
      ylab='Normalized Binary Distance')
  })
  output$indication <- renderDataTable({
    clustNum <- unname(partition[input$drugName])
    DT::datatable(indicationList[[clustNum]], options = list(pageLength = 5))
  })
  
  # Dynamic disease selection
  output$diseaseSelect <- renderUI({
    selectInput(
      inputId = 'diseaseName',
      label = 'Choose a disease from the dropdown menu:',
      choices = diseaseNames,
      selected = 'Diabetes mellitus'
    )
  })
  
  # Dynamic disease selection
  output$clusterSelect <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(input$diseaseName)) return()
    
    disease <- input$diseaseName
    clusters <- which(lapply(indicationList, function(x) disease %in% as.character(x$Indication)) == T)
    clusters <- clusters[!(clusters == 9 | clusters == 10 | clusters == 11)]
    clusters <- as.numeric(clusters)
    
    selectInput(
      inputId = 'clusterNumber',
      label = 'Choose a cluster from the dropdown menu:',
      choices = clusters,
      selected='29'
    )
  })
  
  # Cluster drugs display
  
  output$clustergram <- renderPlot({
    clustNum <- as.integer(input$clusterNumber)
    if (is.na(clustNum)) return()
    plot(
      hclustList[[clustNum]],
      main=paste('Plotting Cluster', clustNum), 
      xlab=NA,
      sub=NA,
      ylab='Normalized Binary Distance')
  })
}

## Compile ##
shinyApp(ui=ui, server=server)