
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shinydashboard)
library(visNetwork)
library(shinyBS)

sidebar <- dashboardSidebar(width = 150,
                            sidebarMenu(id = 'tabs',
                                        menuItem("Targets",
                                                 tabName = "targetmanifest", 
                                                 icon = icon("list")),
                                        menuItem("Details",
                                                 icon = icon("info-sign", lib="glyphicon"),
                                                 tabName = "targetdetails"),
                                        menuItem("Help",
                                                 tabName = "help", 
                                                 icon = icon("question-sign", lib="glyphicon")),
                                        menuItem("About", 
                                                 tabName = "feedback", 
                                                 icon = icon("flag", lib="glyphicon"))
                            )
)

dashboardPage(
  skin = "blue",
  
  
  dashboardHeader(title = "Wall Of Targets",
                  dropdownMenuOutput("notificationMenu")),
  sidebar,
  # body

  dashboardBody(
    
    useShinyjs(),
    
    tags$head(
      singleton(
        includeScript("www/readCookie.js")
      )
    ),
    
    tabItems(
      tabItem(tabName = "targetmanifest",
              fluidRow(
                column(width=2),
                column(width=8,
                       includeMarkdown('welcome.md'),
                       tabBox(id='inputgene', width=NULL,
                              tabPanel("Gene Search", id="selectAGene",
                                       uiOutput("selectGeneBox"),
                                       actionButton('selectGeneBoxButton', 'Go')),
                              tabPanel("Nominated Target Genes", id="nominatedTargets",
                                       includeMarkdown('info.md'),
                                       DT::dataTableOutput('targetlist', width='100%'))
                              )
                ),
                column(width=2))
      ),
      tabItem(tabName = "targetdetails",
              
              # Boxes need to be put in a row (or column)
              fluidRow(
                column(width=12,
                       box(actionButton('targetlist', 'Back to Target List',
                                        icon("paper-plane"), 
                                        style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                           width=NULL))
              ),
              fluidRow(
                column(width=12,
                       box(title=tagList("Target overview",
                                         tipify(icon("question-sign", lib="glyphicon"), "Annotations associated with this gene from Gene Ontology")),
                           solidHeader=TRUE, status="primary",
                           width=NULL, collapsible = TRUE, collapsed = FALSE,
                           uiOutput('targetInfo')
                       ),
                       box(title=tagList("Gene Ontology annotations",
                                         tipify(icon("question-sign", lib="glyphicon"),
                                                "Annotations associated with this gene from Gene Ontology")),
                           solidHeader=TRUE, status="primary",
                           width=NULL, collapsible = TRUE, collapsed = TRUE,
                           DT::dataTableOutput('gomf')
                       ),
                       box(title=tagList("Biological pathways",
                                         tipify(icon("question-sign", lib="glyphicon"), "Curated biological pathways that this gene is involved with.")),
                           solidHeader=TRUE, status="primary",
                           width=NULL, collapsible = TRUE, collapsed = TRUE,
                           DT::dataTableOutput('reactome')),
                       # box(title=tagList("Baseline RNA Expression",
                       #                   tipify(icon("question-sign", lib="glyphicon"), "RNA-Seq median expression values in different brain regions (red line indicates median of all genes)")),
                       #     solidHeader=TRUE, status="primary",
                       #     width=NULL, height = NULL, collapsible = TRUE, collapsed = TRUE,
                       #     # plotOutput("gtex"),
                       #     imageOutput("gtex"),
                       #     textOutput("gtexText")
                       # ),
                       # box(title=tagList("Gene SNPs and structural variants",
                       #                   tipify(icon("question-sign",
                       #                               lib="glyphicon"),
                       #                          title="SNPs and structural variation affecting the genome region containing the gene.")),
                       #     solidHeader=TRUE, collapsible = TRUE, collapsed = TRUE,
                       #     status="primary", width=NULL),
                       # box(title=tagList("Protein information",
                       #                   tipify(icon("question-sign", lib="glyphicon"), "Details on the protein product of this gene.")),
                       #     solidHeader=TRUE, status="primary",
                       #     width=NULL, collapsible = TRUE, collapsed = TRUE),
                       box(title=tagList("Mouse model systems",
                                         tipify(icon("question-sign",
                                                     lib="glyphicon"),
                                                title="Availability of mouse models for the selected gene from the International Mouse Strain Resource (IMSR).")),
                           solidHeader=TRUE, collapsible = TRUE, collapsed = TRUE,
                           status="success", width=NULL,
                           DT::dataTableOutput("ISMR")),
                       # box(title=tagList("Assays",
                       #                   tipify(icon("question-sign",
                       #                               lib="glyphicon"),
                       #                          title="Availability of in vivo or in vitro assays.")),
                       #     solidHeader=TRUE, collapsible = TRUE, collapsed = TRUE,
                       #     status="success", width=NULL),
                       # box(title=tagList("Druggability",
                       #                   tipify(icon("question-sign",
                       #                               lib="glyphicon"),
                       #                          title="Evidence of target drugability from the Oxford Drug Discovery Institute and Lilly/EBI DrugEBIlity")),
                       #     solidHeader=TRUE, collapsible = TRUE, collapsed = TRUE,
                       #     status="success", width=NULL,
                       #     tagList(h3("Oxford Drug Discovery Institute Scores",
                       #                tipify(icon("question-sign",
                       #                            lib="glyphicon"),
                       #                       title="Evidence of target drugability from the Oxford Drug Discovery Institute and Lilly/EBI DrugEBIlity"))),
                       #     plotOutput("status", height=150, width=450),
                       #     h3("Lilly / EBI DrugEBIlity Scores"),
                       #     plotOutput('lillyDrugability', height=150, width=300)
                       # ),
                       # box(title=tagList("Nomination Evidence",
                       #                   tipify(icon("question-sign", lib="glyphicon"),
                       #                          title="A presentation from the group who nominated the target.")), 
                       #     solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                       #     status="danger", width=NULL, 
                       #     DT::dataTableOutput("evidence"),
                       #     htmlOutput('video')),
                       # box(title=tagList("Gene Co-Expression Network",
                       #                   tipify(icon("question-sign", lib="glyphicon"), 
                       #                          title="Co-expression networks from the ROSMAP study.")),
                       #     solidHeader=TRUE, collapsible = TRUE, collapsed = TRUE,
                       #     status="danger", width=NULL, visNetworkOutput("network", height = "350px")),
                       box(title=tagList("Case-Control Differential Expression", 
                                         tipify(icon("question-sign", lib="glyphicon"),
                                                title="Differential gene expression of target genes between individuals with AD and no cognitive impairment (NCI).")),
                           solidHeader=TRUE, collapsible = TRUE, collapsed = TRUE,
                           status="danger", width=NULL, 
                           plotOutput("expression")),
                       box(title=tagList("Expression Forest Plot", 
                                         tipify(icon("question-sign", lib="glyphicon"),
                                                title="Forest plot of 95% confidence intervals around the log fold change across studies and comparison models.")),
                           solidHeader=TRUE, collapsible = TRUE, collapsed = TRUE,
                           status="danger", width=NULL,
                           uiOutput('selectForestPlot'),
                           actionButton('reset','Reset Filters'),
                           actionButton('plotForest', 'Plot'),
                           plotOutput("forest"))
                )
              )
      ),
      tabItem(tabName = "help",
              includeMarkdown("overview.md")
      ),
      tabItem(tabName = "feedback",
              includeMarkdown("feedback.md")
      )
      
    )
  )
)
