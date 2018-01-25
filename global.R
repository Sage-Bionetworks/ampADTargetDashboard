library(synapseClient)
library(data.table)
library(mygene)
library(igraph)
# library(plyr)
# library(dplyr)
# library(stringr)
# library(forcats)
library(tidyverse)
library(ggplot2)
library(rjson)
library(shinyjs)
library(plotly)
library(feather)

usePublic <- FALSE

oddiStatusColors <- c("good"="#5e933f", "medium"="#ef7d0b", "bad"="#a30b0d", "unknown"="#a3a3a3")
lillyStatusColors <- c("3"="#5e933f", "2"="yellow", "1"="#ef7d0b", "0"="#a30b0d", "NA"="#a3a3a3", "unk"="#a3a3a3")

targetManifsetSketch <- htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th('Gene', title = 'The target gene'),
      th('Center(s)', title = 'The center(s) nominating this target'),
      th('Nominations', title= 'Number of centers nominating this target')
    )
  )
))
vids <- c("UFL-ISB-Mayo"="https://s3.amazonaws.com/static.synapse.org/kdaily/AMP-AD/AMP-AD_ExperimentalValidationWGWebinar_092916.mp4",
          "Broad-Rush"="https://s3.amazonaws.com/static.synapse.org/kdaily/AMP-AD/AMP-AD_ExperimentalValidationWGWebinar_101316.mp4",
          "Emory"="https://s3.amazonaws.com/static.synapse.org/kdaily/AMP-AD/AMP-AD_ExperimentalValidationWGWebinar_10616.mp4",
          "MSSM"="https://s3.amazonaws.com/static.synapse.org/kdaily/AMP-AD/AMP-AD_ExperimentalValidationWGWebinar_1473951893.mp4")
