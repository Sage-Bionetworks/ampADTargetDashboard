library(synapseClient)
library(data.table)
library(mygene)
library(plyr)
library(dplyr)
library(igraph)

synapseLogin()

vids <- c("https://www.youtube.com/embed/gc3Kd4ez1iY",
          "https://www.youtube.com/embed/JHCMnr1bU2Q",
          "https://www.youtube.com/embed/Ul0HXq5Tmok",
          "https://www.youtube.com/embed/fZ-N5qqZTGU",
          "https://www.youtube.com/embed/DO5f5R4qW1s"
)

ddiData <- fread(getFileLocation(synGet("syn7537835")), 
                 data.table=FALSE)

lillyData <- fread(getFileLocation(synGet('syn7525109')),
                   data.table=FALSE)

geneTargetList <- sort(unique(c(ddiData$GENE_SYMBOL, lillyData$GENE_SYMBOL)))

out <- queryMany(unique(c(ddiData$GENE_SYMBOL, lillyData$GENE_SYMBOL)),
                        scopes="symbol", fields="ensembl.gene", species="human",
                 returnall=TRUE, size=1)

res <- as.data.frame(out$response)

# Two genes - CHRH1 (in DDI and Lilly) and MOAP1 (in Lilly) - have multiple
# matches to some other chromosomes in Ensembl. Get only the first one.
res <- res %>% mutate(ensembl=laply(ensembl, function(x) ifelse(is.null(x), '', x[[1]][1])))
res[is.na(res$ensembl.gene), 'ensembl.gene'] <- res[is.na(res$ensembl.gene), 'ensembl']
res <- res %>% dplyr::select(-ensembl, -X_id)

ddiData <- ddiData %>% left_join(res, by=c('GENE_SYMBOL'='query'))
lillyData <- lillyData %>% left_join(res, by=c('GENE_SYMBOL'='query'))

targetManifest <- ddiData %>%
  arrange(GENE_SYMBOL) %>%
  select(GENE_SYMBOL, Center, activity_direction)

network <- fread(getFileLocation(synGet("syn7537683")), 
                 data.table=FALSE)

genesForNetwork <- network %>%
  tidyr::gather(source, gene, var1, var2) %>% 
  dplyr::select(gene) %>% 
  mutate(id=gene) %>% 
  distinct() %>% 
  arrange(gene)

out <- queryMany(unique(genesForNetwork$gene),
                 scopes="ensembl.gene", fields="symbol", species="human",
                 returnall=TRUE, size=1)

res <- as.data.frame(out$response) %>% select(symbol, query)

genesForNetwork <- genesForNetwork %>% 
  left_join(res, by=c('id'='query')) %>% 
  mutate(label=ifelse(is.na(symbol), gene, symbol))

gg <- graph_from_data_frame(network)
