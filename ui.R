
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shinydashboard)
library(visNetwork)

sidebar <- dashboardSidebar(
  sidebarMenu(id = 'tabs',
    menuItem("Target Manifest", tabName = "targetmanifest", icon = icon("dashboard")),
    menuItem("Target Details", icon = icon("th"), tabName = "targetdetails")
  )
)

dashboardPage(
  dashboardHeader(title = "AMP-AD Targets"),
  sidebar,
  # body

  dashboardBody(
    tabItems(
      tabItem(tabName = "targetmanifest",
              dataTableOutput('targetlist'),
              actionButton('getdetails', 'Get Target Details')),
      tabItem(tabName = "targetdetails",
              
              # Boxes need to be put in a row (or column)
              fluidRow(
                column(width=6,
                       # box(width=NULL, selectInput("gene", "Gene", choices=geneTargetList)),
                       infoBoxOutput('targetInfo', width = NULL),
                       valueBoxOutput("status", width=NULL),
                       box(title="Nomination Video", solidHeader = TRUE, status="info", 
                           width=NULL, htmlOutput('video'))
                ),
                column(width=6,
                       box(width=NULL, visNetworkOutput("network", height = "350px")),
                       DT::dataTableOutput('edgeTable')
                )
              )
      )
    )
  )
)