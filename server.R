library(shiny)
library(visNetwork)
library(igraph)
library(wesanderson)

shinyServer(function(input, output, session) {
  addClass(selector = "body", class = "sidebar-collapse")
  
  session$sendCustomMessage(type="readCookie",
                            message=list(name='org.sagebionetworks.security.user.login.token'))

  foo <- observeEvent(input$cookie, {

    #synLogin(sessionToken=input$cookie)
    synLogin(silent=TRUE)
    withProgress(message = 'Loading data...',
                 {source("load.R")})
    
    gene <- reactiveValues(geneName=NULL)

    observeEvent(input$targetlist_rows_selected, {
      gene$geneName <- targetManifest[as.numeric(input$targetlist_rows_selected), ]$Gene
    })
    
    observeEvent(input$selectGeneBoxButton, {
      gene$geneName <- input$inputSelectedGene

    })
    
    selectedGene <- reactive({
      gene$geneName
    })
    
    observeEvent(input$targetlist, {
      updateTabItems(session, "tabs", "targetmanifest")
    })

    observeEvent(input$selectGeneBoxButton, {
      updateTabItems(session, "tabs", "targetdetails")
    })
    
    output$targetlist <- DT::renderDataTable(targetManifestTable,
                                             server = TRUE,
                                             container=targetManifsetSketch)

    output$gomf <- DT::renderDataTable({
      
      geneName <- selectedGene()
      ensGene <- filter(geneDF, Gene == geneName)$ensembl.gene[1]
      
      res <- mygene::query(ensGene, fields = c('go.MF'))$hits$go$MF[[1]] %>%
        dplyr::select(term, id)
      
      DT::datatable(res, options = list(lengthChange = FALSE, dom="tp", pageLength=5),
                    rownames = FALSE)
      
    })

    output$reactome <- DT::renderDataTable({
      
      geneName <- selectedGene()
      ensGene <- filter(geneDF, Gene == geneName)$ensembl.gene[1]
      
      res <- mygene::query(ensGene, fields = c('pathway.reactome'))
      
      validate(need(!is.null(res$hits$pathway$reactome), 
                    "No Reactome pathways found."))
      
      res <- res$hits$pathway$reactome[[1]] %>% 
        dplyr::select(name, id)
      
      DT::datatable(res, options = list(lengthChange = FALSE, dom="tp", pageLength=5),
                    rownames = FALSE)
      
    })

    output$evidence <- DT::renderDataTable({
      
      geneName <- selectedGene()
      
      res <- targetListOrig %>% 
        dplyr::filter(gene_symbol == geneName) %>% 
        dplyr::select(group,rank,evidence)

      DT::datatable(res, options = list(lengthChange = FALSE, dom="tp", pageLength=10),
                    rownames = FALSE)
      
    })

    output$IMSR <- DT::renderDataTable({
      
      geneName <- selectedGene()

      res <- IMSR %>% filter(`Gene Symbol` == geneName) %>% 
        dplyr::select(-`Gene Symbol`, -URL) %>% 
        dplyr::select(Repository, dplyr::everything())
      
      DT::datatable(res, options = list(lengthChange = FALSE, 
                                        dom="tp", pageLength=5),
                    rownames = FALSE, escape=2)
    })
    
    observeEvent(input$targetlist_rows_selected, {
      updateTabItems(session, "tabs", selected = "targetdetails")
    })
    
    output$notificationMenu <- renderMenu({
      prof <- synGetUserProfile()
      id <- prof$ownerId
      name <- prof$userName
      
      url <- sprintf("https://www.synapse.org/#!Profile:%s", id)
      msgs <- list(notificationItem(text=sprintf("Logged in as %s", name),
                                    icon=icon('user-circle'), status = "info",
                                    href=url))
    })

    subnetwork <- reactive({
      ensGene <- filter(nodesForNetwork, label == selectedGene())$id

      gg.neighbors <- ego(gg, 1, V(gg)[V(gg)$name %in% ensGene])

      if (length(gg.neighbors) < 1) {
        foo <- induced_subgraph(gg, vids = c())
      } else {
        gg.neighbors <- gg.neighbors[[1]]
        foo <- induced_subgraph(gg, vids = gg.neighbors)
      }

      foo

    })
    
    output$expression <- renderPlot({

      ensGene <- filter(geneDF, Gene == selectedGene())$ensembl.gene[1]
      
      validate(need(length(ensGene) > 0, "No expression data to display."))
      validate(need(length(ensGene) <= 50, "Too many related genes to display - maximum is 50."))

      tmp <- geneFPKMLong %>% dplyr::filter(ensembl_gene_id %in% ensGene)

      medianTmp <- tmp %>% group_by(hgnc_symbol) %>%
        summarize(median=median(fpkm)) %>%
        arrange(median)

      tmp$hgnc_symbol <- factor(tmp$hgnc_symbol,
                                levels=medianTmp$hgnc_symbol,
                                ordered=TRUE)
      
      p <- ggplot(tmp, aes(x=hgnc_symbol, y=fpkm))
      p <- p + geom_boxplot(aes(fill=`Cognitive Diagnosis`))
      p <- p + scale_fill_manual(values=wes_palette("Chevalier"))
      p <- p + scale_color_manual(values=wes_palette("Chevalier"))
      p <- p + theme_bw() + theme(legend.position="bottom")
      p
    })
    
    output$network <- renderVisNetwork({

      network <- subnetwork()
      gg2 <- network %>% toVisNetworkData()
      myEdges <- gg2$edges %>% dplyr::top_n(50, value)
      
      validate(need(nrow(myEdges) > 0, sprintf("No nodes for the gene '%s'.", selectedGene())))
      
      if (nrow(myEdges) > 50) {
        tmp <- myEdges %>% filter(from == selectedGene() | to == selectedGene())
        maxRowsLeft <- max(0, 50 - nrow(tmp))
        tmp2 <- myEdges %>% anti_join(tmp) %>% arrange(-value) %>% slice(1:maxRowsLeft)
        
        myEdges <- rbind(tmp, tmp2)
        gg2$nodes <- gg2$nodes %>% filter(id %in% myEdges$from | id %in% myEdges$to)
        message(sprintf("Network too large (%s edges) for the gene '%s'; maximum number of edges to show is 50.", nrow(edges), selectedGene()))
      }
      
      nodes <- gg2$nodes %>%
        select(id) %>%
        left_join(nodesForNetwork, by='id') %>%
        select(id, label, title, color) %>%
        dplyr::mutate(group=ifelse(label %in% targetManifest$Gene, "target", "other")) %>%
        dplyr::mutate(group=ifelse(label == selectedGene(), "selected", group))

      n <- visNetwork(nodes, myEdges) %>%
        visPhysics(solver="forceAtlas2Based", stabilization = TRUE) %>%
        visEdges(color='black') %>%
        visLegend(width=0.1) %>%
        visGroups(groupname='selected', shape='star', size=40) %>%
        visGroups(groupname='target', size=25) %>%
        visGroups(groupname='other', size=15)


      n
    })
    
    targetInfoValues <- reactiveValues(summary=NULL)
    
    output$targetInfoSummary <- reactive({
      targetInfoValues$summary
    })
    
    output$targetInfo <- renderUI({
      geneName <- selectedGene()
      
      ensGene <- filter(geneDF, Gene == geneName)$ensembl.gene[1]
      
      geneList <- targetListOrig %>% filter(gene_symbol == geneName)
      
      if (nrow(geneList) > 0) {
        ens <- paste(unique(geneList$ensembl_id), collapse=",")
        centers <- paste(geneList$group, collapse=", ")
      }
      else {
        ens <- ensGene
        centers = ""
      }
      
      res <- mygene::getGene(ensGene, fields = c('name', 'summary'))[[1]]
      
      targetInfoValues$summary <- res$summary
      
      tagList(tags$h3(tags$a(href=sprintf("http://www.genenames.org/cgi-bin/gene_search?search=%s", geneName),
                             target="blank",geneName),
                      sprintf("(%s)", res$name)),
              tags$h4(tags$a(href=sprintf('https://www.targetvalidation.org/target/%s', ens), target="blank", ens),
                      "(Open Targets Platform)"),
              tags$h4(sprintf("Nominated by: %s", centers)),
              wellPanel(tags$p(sprintf("Summary: %s", res$summary),
                               style="white-space:pre-wrap;"))
              )
    })
    
    output$selectGeneBox <- renderUI({
      selectizeInput('inputSelectedGene', label='Gene Symbol', multiple=FALSE,
                     choices=c("Select a gene"="", unique(geneExprData$hgnc_symbol)), 
                     selected=NULL, width="50%", 
                     options=list(highlight=TRUE, openOnFocus=FALSE, 
                                  closeAfterSelect=TRUE, selectOnTab=TRUE))
    })
  
    output$de <- renderText({
      
      params <- data.frame(hgnc_symbol=selectedGene(),
                           tissue_study=input$Tissue, 
                           comparison_model_sex=input$Model) %>% 
        separate(tissue_study, into=c("Tissue", "Study"), sep=", ") %>% 
        separate(comparison_model_sex, into=c("Comparison", "Model", "Sex"), sep=", ")
      
      dForPlot <- inner_join(geneExprData, params)
      isDE <- ifelse(dForPlot$adj.P.Val < 0.05, "is", "is not")
      
      glue::glue("{gene} {isDE} differentially expressed (log fold change = {lfc}; adjusted p-value = {pval}) in the {tissue} from {study} when comparing {model} in {sex}.", 
                 gene=params$hgnc_symbol, isDE=isDE, lfc=round(dForPlot$logFC, 3), pval=round(dForPlot$adj.P.Val, 3),
                 tissue=params$Tissue, model=params$Model, sex=params$Sex, study=params$Study)

    })
    output$forest <- renderPlot({
      params <- data.frame(hgnc_symbol=selectedGene(), 
                           comparison_model_sex=input$Model) %>% 
        separate(comparison_model_sex, into=c("Comparison", "Model", "Sex"), sep=", ") %>% 
        select(-Sex)
      
      dForPlot <- inner_join(geneExprData, params) %>% 
        mutate(study_tissue_sex=paste(Study, Comparison, Tissue, Sex))
      
      p <- ggplot(dForPlot)
      p <- p + geom_point(aes(y=logFC, x=study_tissue_sex))
      p <- p + geom_pointrange(aes(ymax = CI.R, ymin = CI.L, y=logFC, 
                                   x=study_tissue_sex, color=Study))
      p <- p + geom_hline(yintercept = 0, linetype = 2)
      p <- p + coord_flip()
      p <- p + theme_minimal()
      p <- p + theme(axis.text=element_text(size=12), 
                     axis.title=element_text(size=14))
      p <- p + labs(x="Log Fold Change", y=NULL)
      p
    }) 
    
    output$volcanoSelect <- renderUI({
      tagList(
        div(style="display: inline-block;vertical-align:top; width: 150px;",
            selectInput(inputId="Tissue", label="Tissue",
                        choices=tissueStudySelections, multiple = FALSE,
                        selected="TCX, MAYO")),
        div(style="display: inline-block;vertical-align:top; width: 25px;",HTML("<br>")),
        div(style="display: inline-block;vertical-align:top; width: 250px;",
            selectInput(inputId="Model", label="Model",
                        choices=modelSexSelections, multiple = FALSE,
                        selected="AD-CONTROL, Diagnosis, Males and Females"))
      )
    
    })
    
    output$volcano <- renderPlotly({
      
      geneName <- selectedGene()
      
      params <- data.frame(tissue_study=input$Tissue, 
                           model_sex=input$Model) %>% 
        separate(tissue_study, into=c("Tissue", "Study"), sep=", ") %>% 
        separate(model_sex, into=c("Comparison", "Model", "Sex"), sep=", ")
        
      dForPlot <- inner_join(geneExprData, params)
      dIsSelected <- dForPlot %>% filter(hgnc_symbol == geneName)
      
      p <- ggplot(dForPlot)
      p <- p + geom_point(aes(x=logFC, y=neg.log10.adj.P.Val, label=hgnc_symbol), alpha=(1/3))
      p <- p + geom_point(aes(x=logFC, y=neg.log10.adj.P.Val, label=hgnc_symbol), 
                          data=dIsSelected, shape=21, color="red", fill="white", 
                          size=1, stroke=2)
      
      p <- p + theme_minimal()
      p <- p + labs(x="Log Fold Change", y="-log10(Adjusted p-value)")
      
      ggplotly(p) %>% config(displayModeBar = F)
    })
    
    output$video <- renderUI({
      geneName <- selectedGene()
      geneList <- targetListOrig %>% filter(Gene == geneName) %>% 
        dplyr::select(group) %>% top_n(1)
      
      validate(need(geneList$group %in% names(vids), sprintf("No video from %s.", geneList$Center)))

      HTML(sprintf('<a href="%s">Video</a>', 
                   vids[[geneList$group]]))
      
      # HTML(sprintf('<video height="250" controls><source src="%s" type="video/mp4"></video>', 
      #              vids[[geneList$Center]]))
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
    
    # output$gtex <- renderImage({
    #     # When input$n is 3, filename is ./images/image3.jpeg
    #     filename <- normalizePath(file.path('./gxa_static.png'))
    #     
    #     # Return a list containing the filename and alt text
    #     list(src = filename,
    #          height = 400,
    #          alt = 'GXA Static Image')
    #     
    #   }, deleteFile = FALSE)    
    
    # output$gtexText <- renderText({"See more expression data at Expression Atlas. This expression view is provided by Expression Atlas. Please send any queries or feedback to arrayexpress-atlas@ebi.ac.uk."})
    
    # output$edgeTable <- DT::renderDataTable(edges()$edges,
    #                                         options=list(lengthChange=FALSE, pageLength=5, dom="tp"))
    
    # output$lillyDrugability <- renderPlot({
    #   selGene <- selectedGene()
    #   
    #   tmp <- druggabilityData %>%
    #     filter(GENE_SYMBOL == selGene)
    #   
    #   if (nrow(tmp) == 0) {
    #     tmp <- data_frame(GENE_SYMBOL=selGene, 
    #                       Lilly_DrugEBIlity_Consensus_Score='unk',
    #                       `Lilly_GW_Druggability_Structure-based`='unk')
    #   }
    #   tmp <- tmp %>%
    #     select(GENE_SYMBOL, 
    #            `Consensus Score`=Lilly_DrugEBIlity_Consensus_Score,
    #            `Structure-based Score`=`Lilly_GW_Druggability_Structure-based`) %>% 
    #     top_n(1) %>% 
    #     tidyr::gather(key = 'type', value = 'score',
    #                   c(`Consensus Score`,
    #                   `Structure-based Score`)) %>% 
    #     mutate(score=factor(score, levels=c("unk", "0", "1", "2", "3")))
    #   
    #   validate(need(nrow(tmp) > 0, "No data to display."))
    #   
    #   tmp$Center <- NA
    #   ggplot(tmp, aes(x=type, y=Center)) + 
    #     facet_wrap( ~ type, scales="free", ncol=5) +
    #     geom_tile(aes(fill=score)) + 
    #     scale_fill_manual(values=lillyStatusColors, drop = FALSE) + 
    #     theme_bw() + 
    #     theme(axis.text=element_blank(), axis.title=element_blank(),
    #           axis.ticks=element_blank(), strip.text.y=element_text(angle=360),
    #           strip.background=element_rect(fill="white"),
    #           legend.position="bottom")
    # })
    
    # output$status <- renderPlot({
    #   
    #   geneName <- selectedGene()
    #   tmp <- druggabilityData %>%
    #     filter(GENE_SYMBOL == geneName)
    #   
    #   if (nrow(tmp) == 0) {
    #     y <- rbind(tmp, rep('unknown', ncol(tmp)))
    #     colnames(y) <- colnames(tmp)
    #     y$GENE_SYMBOL <- geneName
    #     tmp <- y
    #   }
    # 
    #   tmp <- tmp %>%
    #     select(starts_with("status")) %>% 
    #     top_n(1) %>% 
    #     tidyr::gather(key = 'type', value = 'status', starts_with("status")) %>% 
    #     mutate(status=factor(status, levels=c("good", "medium", "bad", "unknown"), ordered=TRUE))
    #   
    #   validate(need(nrow(tmp) > 0, "No data to display."))
    #   
    #   tmp$type <- forcats::fct_recode(tmp$type, `Known Ligands`="status_known_ligands", 
    #                                   `Crystal Structures`="status_crystal_structure", 
    #                                   Pocket="status_pocket", Assays="status_assays", 
    #                                   `In vivo`="status_in_vivo_work")
    #   tmp$Center <- NA
    #   ggplot(tmp, aes(x=type, y=Center)) + 
    #     facet_wrap( ~ type, scales="free", ncol=5) +
    #     geom_tile(aes(fill=status)) + 
    #     scale_fill_manual(values=oddiStatusColors, drop = FALSE) + 
    #     theme_bw() + 
    #     theme(axis.text=element_blank(), axis.title=element_blank(),
    #           axis.ticks=element_blank(), strip.text.y=element_text(angle=360),
    #           strip.background=element_rect(fill="white"),
    #           legend.position="bottom")
    # })
    
  })
})