
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyServer(function(input, output) {

  output$network <- renderVisNetwork({
    
    edges <- network %>% filter(target == input$gene | feature == input$gene) %>% 
      mutate(from=feature, to=target)

    nodes <- genes %>% filter(gene %in% edges$from | gene %in% edges$to) %>% 
      mutate(color=ifelse(gene == input$gene, "97C1FC", "FFD58F"))
    
    visNetwork(nodes, edges) %>% visEdges(arrows='to')
  })

  output$status <- renderValueBox({
    whichGene <- genes %>% filter(gene==input$gene)
    
    valueBox(subtitle="Votes", value=whichGene$votes, 
             color = whichGene$votesColor)
  })
  
  output$video <- renderUI({
    whichGene <- genes %>% filter(gene==input$gene)
    
    tags$iframe(src=whichGene$vid, height=350, width=623)
  })
  
  

})
