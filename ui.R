
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shinydashboard)
library(visNetwork)

sidebar <- dashboardSidebar(width = 125,
                            sidebarMenu(id = 'tabs',
                                        menuItem("Targets",
                                                 tabName = "targetmanifest", 
                                                 icon = icon("list")),
                                        menuItem("Details",
                                                 icon = icon("info-sign", lib="glyphicon"),
                                                 tabName = "targetdetails"),
                                        menuItem("Help",
                                                 tabName = "help", 
                                                 icon = icon("question-sign", lib="glyphicon")),                                        menuItem("Feedback", 
                                                 tabName = "feedback", 
                                                 icon = icon("flag", lib="glyphicon"))
                                        
                            )
)

dashboardPage(skin = "blue",
  dashboardHeader(title = "AMP-AD Targets"),
  sidebar,
  # body

  dashboardBody(
    tabItems(
      tabItem(tabName = "targetmanifest",
              includeMarkdown('info.md'),
              DT::dataTableOutput('targetlist'),
              actionButton('getdetails', 'Get Target Details')),
      tabItem(tabName = "targetdetails",
              
              # Boxes need to be put in a row (or column)
              fluidRow(
                column(width=6,
                       infoBoxOutput('targetInfo', width = NULL),
                       box(title="ODDI Druggability", solidHeader=TRUE, 
                           status="info", height=200, width=NULL, 
                           plotOutput("status", height=150)),
                       box(title="Lilly DrugEBIlity", solidHeader = TRUE, status="info",
                           width=NULL, valueBoxOutput('lilly')),
                       box(title="GTEx", solidHeader=TRUE, status="info",
                           width=NULL,
                           plotOutput("gtex"))),
                column(width=6,
                       box(title="Nomination Video", solidHeader = FALSE, 
                           status="info", width=NULL, htmlOutput('video')),
                       box(title="Expression", solidHeader=TRUE,
                           status="info", width=NULL, 
                           plotOutput("expression")),
                       box(width=NULL, visNetworkOutput("network", height = "350px"))
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
