
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shinydashboard)
library(visNetwork)

dashboardPage(
  dashboardHeader(title = "AMP-AD Targets"),
  dashboardSidebar(disable = TRUE),
  # body
  dashboardBody(
    # Boxes need to be put in a row (or column)
    fluidRow(
      column(width=4,
             box(width=NULL, selectInput("gene", "Gene", choices=genes$gene)),
             valueBoxOutput("status", width=NULL),
             box(width=NULL, htmlOutput('video'))
      ),
      column(width=8,
             box(width=NULL, visNetworkOutput("network", height = "350px")),
             DT::dataTableOutput('edgeTable')
      )
    )
  )
)
