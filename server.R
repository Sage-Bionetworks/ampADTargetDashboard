
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(visNetwork)
library(igraph)
library(wesanderson)

shinyServer(function(input, output, session) {
  session$sendCustomMessage(type="readCookie",
                            message=list(name='org.sagebionetworks.security.user.login.token'))

  foo <- observeEvent(input$cookie, {
    
    synapseLogin(sessionToken=input$cookie)

    withProgress(message = 'Loading data...',
                 {source("load.R")})
    
    selectedGene <- eventReactive(input$getdetails, {
      targetManifest[as.numeric(input$targetlist_rows_selected), ]$Gene
    })
    
    output$targetlist <- DT::renderDataTable(targetManifest,
                                             options=list(lengthChange=FALSE, 
                                                          pageLength=10, dom="ftp"),
                                             selection = list(mode='single', target='row', selected=1),
                                             server = TRUE,
                                             rownames = FALSE,
                                             container=targetManifsetSketch)
    
    observeEvent(input$getdetails, {
      updateTabItems(session, "tabs", selected = "targetdetails")
    })
    
    output$notificationMenu <- renderMenu({
      prof <- synGetUserProfile()
      id <- prof@ownerId
      name <- prof@userName
      
      url <- sprintf("https://www.synapse.org/#!Profile:%s", id)
      msgs <- list(notificationItem(text=sprintf("Logged in as %s", name),
                                    icon=icon('user-circle'), status = "info",
                                    href=url))
    })
    edges <- reactive({
      ensGene <- filter(druggabilityData, GENE_SYMBOL== selectedGene())$ensembl.gene
      
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
      
      fpkmGenes <- filter(druggabilityData, GENE_SYMBOL == selectedGene())$ensembl.gene
      validate(need(length(fpkmGenes) > 0, "No expression data to display."))
      
      tmp <- gtex %>% dplyr::filter(ensembl.gene %in% fpkmGenes)
      
      p <- ggplot(tmp, aes(x=tissue, y=medianFPKM))
      p <- p + geom_col(aes(fill=tissue))
      p <- p + geom_hline(yintercept = medianGTEx, color='red')
      p <- p + theme_bw()
      p <- p + theme(axis.text.x=element_text(angle=270, vjust = 0),
                     legend.position = "none", 
                     axis.title.x=element_blank())
      p
    })
    
    output$expression <- renderPlot({
      
      gg2 <- edges()
      
      if (nrow(gg2$nodes) == 0) {
        fpkmGenes <- filter(druggabilityData, GENE_SYMBOL == selectedGene())$ensembl.gene
      }
      else {
        fpkmGenes <- gg2$nodes$id
      }
      
      print(fpkmGenes)
      validate(need(length(fpkmGenes) > 0, "No expression data to display."))
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
      p <- p + scale_fill_manual(values=wes_palette("Chevalier"))
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
        dplyr::mutate(group=ifelse(label %in% druggabilityData$GENE_SYMBOL, "target", "other")) %>% 
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
    
    output$lillyConsensus <- renderInfoBox({
      tmp <- druggabilityData %>% 
        filter(GENE_SYMBOL == selectedGene()) %>% 
        slice(1)
      
      validate(need(nrow(tmp) > 0, "No data to display."))
      
      valueBox("Consensus", value=tmp$Lilly_DrugEBIlity_Consensus_Score, 
               color=lillyStatusColors[[as.character(tmp$Lilly_DrugEBIlity_Consensus_Score)]])
    })
    
    output$lillyStructureBased <- renderInfoBox({
      tmp <- druggabilityData %>%
        filter(GENE_SYMBOL == selectedGene()) %>% 
        slice(1)
      
      validate(need(nrow(tmp) > 0, "No data to display."))
      
      valueBox("Structure", value=tmp$`Lilly_GW_Druggability_Structure-based`,
               color=lillyStatusColors[[as.character(tmp$`Lilly_GW_Druggability_Structure-based`
               )]])
    })
    
    output$targetInfo <- renderInfoBox({
      geneName <- selectedGene()
      geneList <- targetList %>% filter(Gene == geneName)
      ens <- paste(unique(geneList$ensembl.gene), collapse=",")
      
      infoBox("Selected Target", value=HTML(sprintf("<a href='http://www.genenames.org/cgi-bin/gene_search?search=%s' target='_blank'>%s</a> (HGNC) <a href='https://www.targetvalidation.org/target/%s' target='_blank'>%s</a> (Open Targets Platform)<br/>Nominated by: %s", 
                                                    geneName, geneName, ens, ens, paste(geneList$Center, collapse=","))), 
              icon=icon('info-circle'))
    })
    
    output$status <- renderPlot({
      tmp <- druggabilityData %>%
        filter(GENE_SYMBOL == selectedGene()) %>%
        select(starts_with("status")) %>% 
        top_n(1) %>% 
        tidyr::gather(key = 'type', value = 'status', starts_with("status")) %>% 
        mutate(status=factor(status, levels=c("good", "medium", "bad", "unknown")))
      
      validate(need(nrow(tmp) > 0, "No data to display."))
      
      tmp$type <- forcats::fct_recode(tmp$type, `Known Ligands`="status_known_ligands", 
                                      `Crystal Structures`="status_crystal_structure", 
                                      Pocket="status_pocket", Assays="status_assays", 
                                      `In vivo`="status_in_vivo_work")
      tmp$Center <- NA
      ggplot(tmp, aes(x=type, y=Center)) + 
        facet_wrap( ~ type, scales="free", ncol=5) +
        geom_tile(aes(fill=status)) + 
        scale_fill_manual(values=oddiStatusColors) + 
        theme_bw() + 
        theme(axis.text=element_blank(), axis.title=element_blank(),
              axis.ticks=element_blank(), strip.text.y=element_text(angle=360),
              strip.background=element_rect(fill="white"),
              legend.position="bottom")
    })
    
    output$video <- renderUI({
      geneName <- selectedGene()
      geneList <- targetList %>% filter(Gene == geneName) %>% select(Center) %>% top_n(1)
      
      validate(need(geneList$Center %in% names(vids), sprintf("No video from %s.", geneList$Center)))
      
      HTML(sprintf('<video height="250" controls><source src="%s" type="video/mp4"></video>', 
                   vids[[geneList$Center]]))
    })
  })
})