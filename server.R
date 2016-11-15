
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(DT)
library(visNetwork)
library(igraph)

shinyServer(function(input, output, session) {

  selectedGene <- eventReactive(input$getdetails, {
    targetManifest[as.numeric(input$targetlist_rows_selected), ]$GENE_SYMBOL
  })
  
  output$targetlist <- DT::renderDataTable(targetManifest,
                                          options=list(lengthChange=FALSE, 
                                                       pageLength=15, dom="tp"), 
                                          selection = 'single')
  
  observeEvent(input$getdetails, {
    print(row)
    updateTabItems(session, "tabs", selected = "targetdetails")
  })
  
  edges <- reactive({
    ensGene <- c(filter(lillyData, GENE_SYMBOL== selectedGene())$ensembl.gene,
                 filter(ddiData, GENE_SYMBOL== selectedGene())$ensembl.gene)

    
    gg.neighbors <- ego(gg, 1, V(gg)[V(gg)$name %in% ensGene])
    
    if (length(gg.neighbors) < 1) {
      gg.neighbors <- c()
    } else {
      gg.neighbors <- gg.neighbors[[1]]
    }
    
    foo <- induced_subgraph(gg, vids = gg.neighbors) %>%
      toVisNetworkData()
    foo
  })
  
  output$network <- renderVisNetwork({
    
    gg2 <- edges()    
    validate(need(nrow(gg2$edges) > 0, sprintf("No edges for the gene '%s'.", selectedGene())))
    validate(need(nrow(gg2$edges) <= 50, sprintf("Network too large (%s edges) for the gene '%s'; maximum number of edges to show is 50.", nrow(edges), selectedGene())))
    
    nodes <- gg2$nodes %>% 
      select(id) %>% 
      left_join(genesForNetwork, by='id') %>% 
      select(gene, id, label) %>% 
      dplyr::mutate(group=ifelse((label %in% geneTargetList) & (label != selectedGene()),
                                 "target", "other")) %>% 
      dplyr::mutate(group=ifelse(label == selectedGene(), "selected", group))
      
    print(nodes)
    
    n <- visNetwork(nodes, gg2$edges) %>% 
      visPhysics(solver="forceAtlas2Based", stabilization = TRUE) %>% 
      visEdges(color='black') %>% 
      visLegend() %>% 
      visGroups(groupname='selected', color='green') %>% 
      visGroups(groupname='target', color='#97C1FC') %>% 
      visGroups(groupname='other', color='#FFD58F')
      
    # if (nrow(edges) <=10) {
    #   n <- n %>% visIgraphLayout()
    # }
    
    n
  })

  output$edgeTable <- DT::renderDataTable(edges()$edges,
                                          options=list(lengthChange=FALSE, pageLength=5, dom="tp"))
  
  output$targetInfo <- renderInfoBox({
    infoBox("Selected Target", value=selectedGene(), color = 'green')
  })
  
  output$status <- renderValueBox({
    valueBox(subtitle="Votes", value=100, 
             color = 'red')
  })
  
  output$video <- renderUI({
    tags$iframe(src=vids[1], height=300, width=534)
  })
  
  

})
