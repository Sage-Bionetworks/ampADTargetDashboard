
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(DT)
library(visNetwork)
shinyServer(function(input, output) {

  edges <- reactive({
    network %>% filter(target == input$gene | feature == input$gene) %>% 
      mutate(from=feature, to=target)
  })
  
  output$network <- renderVisNetwork({
    
    edges <- edges()
    
    validate(need(nrow(edges) <= 20, sprintf("Network too large (%s edges) for the gene '%s'; maximum number of edges to show is 20.", nrow(edges), input$gene)))
    
    nodes <- genes %>% 
      dplyr::filter(gene %in% edges$from | gene %in% edges$to) %>% 
      dplyr::mutate(color=ifelse(gene == input$gene, "97C1FC", "FFD58F"))
    
    n <- visNetwork(nodes, edges) %>% visEdges(arrows='to')
    if (nrow(edges) <=10) {
      n <- n %>% visHierarchicalLayout()
    }
    
    n
  })

  output$edgeTable <- DT::renderDataTable(edges() %>% select(feature, target, coexpression, feature.fdr, feature.lfc, target.fdr, target.lfc),
                                          options=list(lengthChange=FALSE, pageLength=5, dom="tp"))
  
  output$status <- renderValueBox({
    whichGene <- genes %>% filter(gene==input$gene)
    
    valueBox(subtitle="Votes", value=whichGene$votes, 
             color = whichGene$votesColor)
  })
  
  output$video <- renderUI({
    whichGene <- genes %>% filter(gene==input$gene)
    
    tags$iframe(src=whichGene$vid, height=300, width=534)
  })
  
  

})
