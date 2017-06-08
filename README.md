# AMP-AD Wall of Targets

A R Shiny dashboard of nominated gene targets with information about their genetics, interactions, and druggability.

## Requirements

- Bioconductor 3.4
- `shiny`
- `shinydashboard`
- `synapseClient`
- `data.table`
- `mygene==1.10`
- `plyr`
- `dplyr`
- `igraph`
- `forcats`
- `ggplot2`
- `stringr`
- `rstudio/DT`
- `visNetwork`
- `wesanderson`

## Run

Currently requires to be logged into Synapse to run in the browser.

```
library(shiny)
shiny::runGitHub('Sage-Bionetworks/WallOfTargets')
```
