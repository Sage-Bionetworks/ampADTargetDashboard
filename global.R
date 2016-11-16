library(synapseClient)
library(data.table)
library(mygene)
library(plyr)
library(dplyr)
library(igraph)
library(forcats)
library(ggplot2)
library(stringr)

synapseLogin()

oddiStatusColors <- c("good"="#5e933f", "medium"="#ef7d0b", "bad"="#a30b0d", "unknown"="#a3a3a3")
lillyStatusColors <- c("3"="green", "2"="yellow", "1"="orange", "0"="red")

vids <- c("https://www.youtube.com/embed/gc3Kd4ez1iY",
          "https://www.youtube.com/embed/JHCMnr1bU2Q",
          "https://www.youtube.com/embed/Ul0HXq5Tmok",
          "https://www.youtube.com/embed/fZ-N5qqZTGU",
          "https://www.youtube.com/embed/DO5f5R4qW1s"
)

ddiData <- fread(getFileLocation(synGet("syn7537835")), 
                 data.table=FALSE) %>% 
  mutate(status_assays=fct_recode(status_assays, unknown=""),
         status_crystal_structure=fct_recode(status_crystal_structure, unknown=""),
         status_pocket=fct_recode(status_pocket, unknown=""),
         status_in_vivo_work=fct_recode(status_in_vivo_work, unknown=""),
         status_known_ligands=fct_recode(status_known_ligands, unknown="")
  )

lillyData <- fread(getFileLocation(synGet('syn7525109')),
                   data.table=FALSE)

out <- queryMany(unique(c(ddiData$GENE_SYMBOL, lillyData$GENE_SYMBOL)),
                        scopes="symbol", fields="ensembl.gene", species="human",
                 returnall=TRUE, size=1)

res <- as.data.frame(out$response)

# Two genes - CHRH1 (in DDI and Lilly) and MOAP1 (in Lilly) - have multiple
# matches to some other chromosomes in Ensembl. Get only the first one.
res <- res %>% mutate(ensembl.gene=laply(ensembl, 
                                         function(x) ifelse(is.null(x), 
                                                            '', x[[1]][1])))

# res[is.na(res$ensembl), 'ensembl.gene'] <- res[is.na(res$ensembl), 'ensembl']
res <- res %>% dplyr::select(ensembl.gene, query)

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

geneFPKM <- fread(getFileLocation(synGet("syn5581268")), 
                  data.table=FALSE) %>% 
  filter(ensembl_gene_id %in% c(genesForNetwork$gene, 
                                ddiData$ensembl.gene)) 

geneCovariates <- fread(getFileLocation(synGet("syn5581227")),
                        data.table=FALSE) %>% 
  filter(cogdx %in% c(1, 4)) %>%
  mutate(cogdx=factor(cogdx, ordered=TRUE))

geneFPKMLong <- geneFPKM %>% 
  tidyr::gather(sample, fpkm, 3:ncol(geneFPKM)) %>% 
  left_join(geneCovariates %>% select(Sampleid_batch, cogdx), 
            by=c("sample"="Sampleid_batch")) %>% 
  filter(!is.na(cogdx))

gtexObj <- synGet('syn7542283')

gtex <- fread(getFileLocation(gtexObj), data.table=FALSE) %>% 
  mutate(ensembl.gene=str_replace(Name, "\\..*", "")) %>% 
  dplyr::filter(ensembl.gene %in% c(genesForNetwork$gene, 
                                    ddiData$ensembl.gene)) %>% 
  select(ensembl.gene, hgnc_symbol=Description, starts_with('Brain'))

gtex <- gtex %>% 
  tidyr::gather(tissue, medianFPKM, 3:ncol(gtex)) %>% 
  mutate(tissue=str_replace(tissue, "Brain - ", ""))

