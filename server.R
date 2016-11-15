
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
    ensGene <- c(filter(lillyData, GENE_SYMBOL== input$gene)$ensembl.gene,
                 filter(ddiData, GENE_SYMBOL== input$gene)$ensembl.gene)
    
    network %>% filter(var1 %in% ensGene | var2 %in% ensGene) %>% 
      mutate(from=var1, to=var2)
  })
  
  output$network <- renderVisNetwork({
    
    edges <- edges()
    validate(need(nrow(edges) > 0, sprintf("No edges for the gene '%s'.", input$gene)))
    
    validate(need(nrow(edges) <= 50, sprintf("Network too large (%s edges) for the gene '%s'; maximum number of edges to show is 50.", nrow(edges), input$gene)))
    
    nodes <- genesForNetwork %>% 
      dplyr::filter(gene %in% edges$from | gene %in% edges$to) %>% 
      dplyr::mutate(group=ifelse(label %in% geneTargetList, "target", "other")) %>% 
      dplyr::mutate(group=ifelse(label == input$gene, "selected", "other"))
    
    n <- visNetwork(nodes, edges) %>% 
      visEdges(color='black') %>% 
      visGroups(groupname='selected', color='green') %>% 
      visGroups(groupname='target', color='#97C1FC') %>% 
      visGroups(groupname='other', color='#FFD58F')
      
    # if (nrow(edges) <=10) {
    #   n <- n %>% visIgraphLayout()
    # }
    
    n
  })

  output$edgeTable <- DT::renderDataTable(edges() %>% select(var1, var2),
                                          options=list(lengthChange=FALSE, pageLength=5, dom="tp"))
  
  output$status <- renderValueBox({
    valueBox(subtitle="Votes", value=100, 
             color = 'red')
  })
  
  output$video <- renderUI({
    tags$iframe(src=vids[1], height=300, width=534)
  })
  
  

})
