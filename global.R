library(synapseClient)
library(data.table)
library(mygene)
library(plyr)
library(dplyr)
library(igraph)
library(forcats)
library(ggplot2)
library(stringr)
library(rjson)

oddiStatusColors <- c("good"="#5e933f", "medium"="#ef7d0b", "bad"="#a30b0d", "unknown"="#a3a3a3")
lillyStatusColors <- c("3"="green", "2"="yellow", "1"="orange", "0"="red", "NA"="red", "unk"="red")

targetManifsetSketch <- htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th('Gene', title = 'The target gene'),
      # th('Center', title = 'The center nominating this target'),
      th('ODDI Druggability Score', title = 'Sum of individual ODDI Druggability scores'),
      th('Lilly Druggability Consensus', title = 'A consensus score for druggability from Lilly')
    )
  )
))
vids <- c("UFL-ISB-Mayo"="https://s3.amazonaws.com/static.synapse.org/kdaily/AMP-AD/AMP-AD_ExperimentalValidationWGWebinar_092916.mp4",
          "Broad-Rush"="https://s3.amazonaws.com/static.synapse.org/kdaily/AMP-AD/AMP-AD_ExperimentalValidationWGWebinar_101316.mp4",
          "Emory"="https://s3.amazonaws.com/static.synapse.org/kdaily/AMP-AD/AMP-AD_ExperimentalValidationWGWebinar_10616.mp4",
          "MSSM"="https://s3.amazonaws.com/static.synapse.org/kdaily/AMP-AD/AMP-AD_ExperimentalValidationWGWebinar_1473951893.mp4")
