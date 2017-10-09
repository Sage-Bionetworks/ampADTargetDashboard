
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
  addClass(selector = "body", class = "sidebar-collapse")
  
    # session$sendCustomMessage(type="readCookie",
  #                           message=list(name='org.sagebionetworks.security.user.login.token'))
  # 
  # foo <- observeEvent(input$cookie, {
  # 
  #   synapseLogin(sessionToken=input$cookie)
    synapseLogin()
    withProgress(message = 'Loading data...',
                 {source("load.R")})
    
    
    selectedGene <- eventReactive(input$targetlist_rows_selected, {
      targetManifest[as.numeric(input$targetlist_rows_selected), ]$Gene
    })
    
    observeEvent(input$targetlist, {
      updateTabItems(session, "tabs", "targetmanifest")
    })
    
    output$targetlist <- DT::renderDataTable(targetManifestTable,
                                             server = TRUE,
                                             rownames = FALSE,
                                             container=targetManifsetSketch)

    output$gomf <- DT::renderDataTable({
      
      geneName <- selectedGene()
      ensGene <- filter(druggabilityData, GENE_SYMBOL== geneName)$ensembl.gene[1]
      
      res <- mygene::query(ensGene, fields = c('go.MF'))$hits$go$MF[[1]] %>%
        dplyr::select(term, id)
      
      DT::datatable(res, options = list(lengthChange = FALSE, dom="tp", pageLength=10),
                    rownames = FALSE)
      
      # # DT::datatable(data.frame(a=1, b=2 c=3))#,
      #               # options=list(lengthChange=FALSE, 
      #               #              pageLength=50, dom="ftp"),
      #               # rownames = FALSE)
    })

    output$reactome <- DT::renderDataTable({
      
      geneName <- selectedGene()
      ensGene <- filter(druggabilityData, GENE_SYMBOL== geneName)$ensembl.gene[1]
      
      res <- mygene::query(ensGene, fields = c('pathway.reactome'))
      
      validate(need(!is.null(res$hits$pathway$reactome), 
                    "No Reactome pathways found."))
      
      res <- res$hits$pathway$reactome[[1]] %>% 
        select(name, id)
      
      DT::datatable(res, options = list(lengthChange = FALSE, dom="tp", pageLength=10),
                    rownames = FALSE)
      
    })

    output$evidence <- DT::renderDataTable({
      
      geneName <- selectedGene()
      
      res <- targetListOrig %>% filter(gene_symbol == geneName) %>% select(group,rank,evidence)

      DT::datatable(res, options = list(lengthChange = FALSE, dom="tp", pageLength=10),
                    rownames = FALSE)
      
    })

    output$ISMR <- DT::renderDataTable({
      
      geneName <- selectedGene()

      res <- ISMR %>% filter(`Gene Symbol` == geneName) %>% select(-`Gene Symbol`) %>% 
        dplyr::select(Repository, dplyr::everything())
      
      DT::datatable(res, options = list(lengthChange = FALSE, dom="tp", pageLength=10),
                    rownames = FALSE)
      
      # # DT::datatable(data.frame(a=1, b=2 c=3))#,
      #               # options=list(lengthChange=FALSE, 
      #               #              pageLength=50, dom="ftp"),
      #               # rownames = FALSE)
    })
    
    observeEvent(input$targetlist_rows_selected, {
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
        # gg2 <- make_empty_graph(1)
        # V(gg2)$name <- ensGene
        # foo <- gg2
        foo <- induced_subgraph(gg, vids = c())
      } else {
        gg.neighbors <- gg.neighbors[[1]]
        foo <- induced_subgraph(gg, vids = gg.neighbors)
      }
      
      foo %>%
        toVisNetworkData()
      
    })
    
    # output$gtex <- renderPlot({
    #   
    #   fpkmGenes <- filter(gtex, hgnc_symbol == selectedGene())$ensembl.gene
    #   message(sprintf('fpkmGenes = %s', fpkmGenes))
    #   validate(need(length(fpkmGenes) > 0, "No expression data to display."))
    #   
    #   tmp <- gtex %>% dplyr::filter(ensembl.gene %in% fpkmGenes)
    #   
    #   p <- ggplot(tmp, aes(x=tissue, y=medianFPKM))
    #   p <- p + geom_col(fill="black")
    #   p <- p + geom_hline(yintercept = medianGTEx, color='red')
    #   p <- p + theme_bw()
    #   p <- p + theme(axis.text.x=element_text(angle=270, vjust = 0),
    #                  legend.position = "none", 
    #                  axis.title.x=element_blank())
    #   p
    # })

    output$gtex <- renderImage({
        # When input$n is 3, filename is ./images/image3.jpeg
        filename <- normalizePath(file.path('./gxa_static.png'))
        
        # Return a list containing the filename and alt text
        list(src = filename,
             height = 400,
             alt = 'GXA Static Image')
        
      }, deleteFile = FALSE)    
    
    output$gtexText <- renderText({"See more expression data at Expression Atlas. This expression view is provided by Expression Atlas. Please send any queries or feedback to arrayexpress-atlas@ebi.ac.uk."})
    
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
      tmp <- tmp %>% rename(`Cognitive Diagnosis`=cogdx)
      p <- ggplot(tmp, aes(x=hgnc_symbol, y=fpkm))
      p <- p + geom_boxplot(aes(fill=`Cognitive Diagnosis`))
      p <- p + scale_fill_manual(values=wes_palette("Chevalier"))
      p <- p + theme_bw() + theme(legend.position="bottom")
      p
    })
    
    output$network <- renderVisNetwork({
      
      gg2 <- edges()    
      validate(need(nrow(gg2$edges) > 0, sprintf("No nodes for the gene '%s'.", selectedGene())))
      validate(need(nrow(gg2$edges) <= 50, sprintf("Network too large (%s edges) for the gene '%s'; maximum number of edges to show is 50.", nrow(edges), selectedGene())))
      
      nodes <- gg2$nodes %>% 
        select(id) %>% 
        left_join(genesForNetwork, by='id') %>% 
        select(gene, id, label) %>% 
        dplyr::mutate(group=ifelse(label %in% targetManifest$Gene, "target", "other")) %>% 
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

    output$lillyDrugability <- renderPlot({
      selGene <- selectedGene()
      
      tmp <- druggabilityData %>%
        filter(GENE_SYMBOL == selGene)
      
      if (nrow(tmp) == 0) {
        tmp <- data_frame(GENE_SYMBOL=selGene, 
                          Lilly_DrugEBIlity_Consensus_Score='unk',
                          `Lilly_GW_Druggability_Structure-based`='unk')
      }
      tmp <- tmp %>%
        select(GENE_SYMBOL, 
               `Consensus Score`=Lilly_DrugEBIlity_Consensus_Score,
               `Structure-based Score`=`Lilly_GW_Druggability_Structure-based`) %>% 
        top_n(1) %>% 
        tidyr::gather(key = 'type', value = 'score',
                      c(`Consensus Score`,
                      `Structure-based Score`)) %>% 
        mutate(score=factor(score, levels=c("unk", "0", "1", "2", "3")))
      
      validate(need(nrow(tmp) > 0, "No data to display."))
      
      tmp$Center <- NA
      ggplot(tmp, aes(x=type, y=Center)) + 
        facet_wrap( ~ type, scales="free", ncol=5) +
        geom_tile(aes(fill=score)) + 
        scale_fill_manual(values=lillyStatusColors, drop = FALSE) + 
        theme_bw() + 
        theme(axis.text=element_blank(), axis.title=element_blank(),
              axis.ticks=element_blank(), strip.text.y=element_text(angle=360),
              strip.background=element_rect(fill="white"),
              legend.position="bottom")
    })

    output$targetInfo <- renderUI({
      geneName <- selectedGene()
      
      ensGene <- filter(druggabilityData, GENE_SYMBOL== geneName)$ensembl.gene[1]
      
      geneList <- targetList %>% filter(Gene == geneName)
      ens <- paste(unique(geneList$ensembl.gene), collapse=",")
      
      res <- mygene::getGene(ensGene, fields = c('name', 'summary'))[[1]]
      
      tagList(tags$h3(tags$a(href=sprintf("http://www.genenames.org/cgi-bin/gene_search?search=%s", geneName),
                             target="blank",geneName),
                      sprintf("(%s)", res$name)),
              tags$h4(tags$a(href=sprintf('https://www.targetvalidation.org/target/%s', ens), target="blank", ens),
                      "(Open Targets Platform)"),
              tags$h4(sprintf("Nominated by: %s", paste(geneList$Center, collapse=", "))),
              tags$p(sprintf("Summary: %s", res$summary))
              )
    })
    
    output$status <- renderPlot({
      
      geneName <- selectedGene()
      tmp <- druggabilityData %>%
        filter(GENE_SYMBOL == geneName)
      
      if (nrow(tmp) == 0) {
        y <- rbind(tmp, rep('unknown', ncol(tmp)))
        colnames(y) <- colnames(tmp)
        y$GENE_SYMBOL <- geneName
        tmp <- y
      }

      tmp <- tmp %>%
        select(starts_with("status")) %>% 
        top_n(1) %>% 
        tidyr::gather(key = 'type', value = 'status', starts_with("status")) %>% 
        mutate(status=factor(status, levels=c("good", "medium", "bad", "unknown"), ordered=TRUE))
      
      validate(need(nrow(tmp) > 0, "No data to display."))
      
      tmp$type <- forcats::fct_recode(tmp$type, `Known Ligands`="status_known_ligands", 
                                      `Crystal Structures`="status_crystal_structure", 
                                      Pocket="status_pocket", Assays="status_assays", 
                                      `In vivo`="status_in_vivo_work")
      tmp$Center <- NA
      ggplot(tmp, aes(x=type, y=Center)) + 
        facet_wrap( ~ type, scales="free", ncol=5) +
        geom_tile(aes(fill=status)) + 
        scale_fill_manual(values=oddiStatusColors, drop = FALSE) + 
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

      HTML(sprintf('<a href="%s">Video</a>', 
                   vids[[geneList$Center]]))
      
      # HTML(sprintf('<video height="250" controls><source src="%s" type="video/mp4"></video>', 
      #              vids[[geneList$Center]]))
    })
   ## })
})