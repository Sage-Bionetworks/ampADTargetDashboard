
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
  
  
  dashboardHeader(title = "wall Of Targets",
                  dropdownMenuOutput("notificationMenu")),
  sidebar,
  # body

  dashboardBody(
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
                       includeMarkdown('info.md'),
                       DT::dataTableOutput('targetlist', width='100%')
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
                           htmlOutput('targetInfo')
                       ),
                       box(title=tagList("Protein information",
                                         tipify(icon("question-sign", lib="glyphicon"), "Details on the protein product of this gene.")),
                           solidHeader=TRUE, status="primary",
                           width=NULL, collapsible = TRUE, collapsed = TRUE),

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
                           width=NULL, collapsible = TRUE, collapsed = TRUE),
                       box(title=tagList("Baseline RNA Expression",
                                         tipify(icon("question-sign", lib="glyphicon"), "RNA-Seq median expression values in different brain regions (red line indicates median of all genes)")),
                           solidHeader=TRUE, status="primary",
                           width=NULL, collapsible = TRUE, collapsed = TRUE,
                           plotOutput("gtex")),
                       box(title=tagList("Mouse model systems",
                                         tipify(icon("question-sign",
                                                     lib="glyphicon"),
                                                title="Availability of mouse model systems for AD.")),
                           solidHeader=TRUE, collapsible = TRUE, collapsed = TRUE,
                           status="success", width=NULL),
                       box(title=tagList("Assays",
                                         tipify(icon("question-sign",
                                                     lib="glyphicon"),
                                                title="Availability of in vivo or in vitro assays.")),
                           solidHeader=TRUE, collapsible = TRUE, collapsed = TRUE,
                           status="success", width=NULL),
                       box(title=tagList("ODDI Druggability",
                                         tipify(icon("question-sign",
                                                     lib="glyphicon"),
                                                title="Evidence of target drugability from ODDI.")),
                           solidHeader=TRUE, collapsible = TRUE, collapsed = TRUE,
                           status="success", width=NULL, 
                           plotOutput("status", height=150)),
                       box(title=tagList("Lilly DrugEBIlity",
                                         tipify(icon("question-sign", 
                                                     lib="glyphicon"),
                                                title="Evidence of target drugability from Lilly.")),
                           solidHeader = TRUE, status="success",
                           width=NULL, collapsible = TRUE, collapsed = TRUE,
                           plotOutput('lillyDrugability', height=150)
                       ),
                       box(title=tagList("Nomination Video",
                                         tipify(icon("question-sign", lib="glyphicon"),
                                                title="A presentation from the group who nominated the target.")), 
                           solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                           status="danger", width=NULL, htmlOutput('video')),
                       box(title=tagList("Gene Co-Expression Network",
                                         tipify(icon("question-sign", lib="glyphicon"), 
                                                title="Co-expression networks from the ROSMAP study.")),
                           solidHeader=TRUE, collapsible = TRUE, collapsed = TRUE,
                           status="danger", width=NULL, visNetworkOutput("network", height = "350px")),
                       box(title=tagList("Case-Control Differential Expression", 
                                         tipify(icon("question-sign", lib="glyphicon"),
                                                title="Differential gene expression of target genes between individuals with AD and no cognitive impairment (NCI).")),
                           solidHeader=TRUE, collapsible = TRUE, collapsed = TRUE,
                           status="danger", width=NULL, 
                           plotOutput("expression"))
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
