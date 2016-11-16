
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(visNetwork)
library(igraph)

shinyServer(function(input, output, session) {

  selectedGene <- eventReactive(input$getdetails, {
    targetManifest[as.numeric(input$targetlist_rows_selected), ]$GENE_SYMBOL
  })
  
  output$targetlist <- DT::renderDataTable(targetManifest %>% rename(`Gene`=GENE_SYMBOL),
                                          options=list(lengthChange=FALSE, 
                                                       pageLength=15, dom="tp"),
                                          selection = list(mode='single', target='row', selected=1),
                                          server = TRUE,
                                          rownames = FALSE)
  
  observeEvent(input$getdetails, {
    updateTabItems(session, "tabs", selected = "targetdetails")
  })
  
  edges <- reactive({
    ensGene <- filter(ddiData, GENE_SYMBOL== selectedGene())$ensembl.gene
    
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

  output$gtex <- renderPlot({
    
    fpkmGenes <- filter(ddiData, GENE_SYMBOL == selectedGene())$ensembl.gene
    
    tmp <- gtex %>% dplyr::filter(ensembl.gene %in% fpkmGenes)
    
    p <- ggplot(tmp, aes(x=tissue, y=medianFPKM))
    p <- p + geom_col(aes(fill=tissue))
    p <- p + theme_bw()
    p <- p + theme(axis.text.x=element_text(angle=270, vjust = 0),
                   legend.position = "none", 
                   axis.title.x=element_blank())
    p
  })
  
  output$expression <- renderPlot({
    
    gg2 <- edges()

    if (nrow(gg2$nodes) == 0) {
      fpkmGenes <- filter(ddiData, GENE_SYMBOL == selectedGene())$ensembl.gene
    }
    else {
      fpkmGenes <- gg2$nodes$id
    }
    
    print(fpkmGenes)
    validate(need(length(fpkmGenes) <= 50, "Too many related genes to display - maximum is 50."))
    
    tmp <- geneFPKMLong %>% dplyr::filter(ensembl_gene_id %in% fpkmGenes)
    
    medianTmp <- tmp %>% group_by(hgnc_symbol) %>% 
      summarize(median=median(fpkm)) %>% 
      arrange(median)
    
    tmp$hgnc_symbol <- factor(tmp$hgnc_symbol, 
                              levels=medianTmp$hgnc_symbol,
                              ordered=TRUE)

    p <- ggplot(tmp, aes(x=hgnc_symbol, y=fpkm))
    p <- p + geom_boxplot(aes(fill=cogdx))
    p <- p + scale_fill_manual(values=c("1"="blue", "4"="orange"))
    p <- p + theme_bw()
    p
  })
  
  output$network <- renderVisNetwork({
    
    gg2 <- edges()    
    validate(need(nrow(gg2$edges) > 0, sprintf("No edges for the gene '%s'.", selectedGene())))
    validate(need(nrow(gg2$edges) <= 50, sprintf("Network too large (%s edges) for the gene '%s'; maximum number of edges to show is 50.", nrow(edges), selectedGene())))
    
    nodes <- gg2$nodes %>% 
      select(id) %>% 
      left_join(genesForNetwork, by='id') %>% 
      select(gene, id, label) %>% 
      dplyr::mutate(group=ifelse(label %in% ddiData$GENE_SYMBOL, "target", "other")) %>% 
      dplyr::mutate(group=ifelse(label == selectedGene(), "selected", group))
      
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

  output$lilly <- renderInfoBox({
    tmp <- lillyData %>% 
      filter(GENE_SYMBOL == selectedGene()) %>% 
      select(Score=Lilly_DrugEBIlity_Consensus_Score)
    
    valueBox("Score",value = tmp$Score, 
             color=lillyStatusColors[[as.character(tmp$Score)]])
  })
  
  output$targetInfo <- renderInfoBox({
    geneName <- selectedGene()
    geneList <- ddiData %>% filter(GENE_SYMBOL == geneName)
    ens <- paste(geneList$ensembl.gene, collapse=",")
    
    infoBox("Selected Target", value=HTML(sprintf("%s<br/>%s", geneName, ens)), color = 'green')
  })
  
  output$status <- renderPlot({
    tmp <- ddiData %>%
      filter(GENE_SYMBOL == selectedGene()) %>%
      select(Center, starts_with("status")) %>% 
      tidyr::gather(key = 'type', value = 'status', starts_with("status")) %>% 
      mutate(status=factor(status, levels=c("good", "medium", "bad", "unknown")))
    
    tmp$type <- forcats::fct_recode(tmp$type, `Known Ligands`="status_known_ligands", 
                                    `Crystal Structures`="status_crystal_structure", 
                                    Pocket="status_pocket", Assays="status_assays", 
                                    `In vivo`="status_in_vivo_work")
    
    ggplot(tmp, aes(x=type, y=Center)) + 
      facet_grid(Center ~ type, scales="free") +
      geom_tile(aes(fill=status)) + 
      scale_fill_manual(values=oddiStatusColors) + 
      theme_bw() + 
      theme(axis.text=element_blank(), axis.title=element_blank(),
            axis.ticks=element_blank(), strip.text.y=element_text(angle=360),
            strip.background=element_rect(fill="white"),
            legend.position="bottom")
  })
  
  output$video <- renderUI({
    tags$iframe(src=vids[1], height=300, width=534)
  })
  
  

})
