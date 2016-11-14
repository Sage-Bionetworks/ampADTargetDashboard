library(synapseClient)
library(data.table)
library(dplyr)

synapseLogin()

vids <- c("https://www.youtube.com/embed/gc3Kd4ez1iY",
          "https://www.youtube.com/embed/JHCMnr1bU2Q",
          "https://www.youtube.com/embed/Ul0HXq5Tmok",
          "https://www.youtube.com/embed/fZ-N5qqZTGU",
          "https://www.youtube.com/embed/DO5f5R4qW1s"
)

ddiData <- fread("~/Projects/AMP-AD/DDI_medium_list_processed.csv", data.table=FALSE) # synGet()
out <- queryMany(ddiData$GENE_SYMBOL, scopes="symbol", fields="ensembl.gene", species="human",
                 returnall=TRUE)

network <- fread(getFileLocation(synGet("syn7346460", version=11)), 
                 data.table=FALSE) %>% 
  dplyr::filter(feature.assay=='TF', target.assay=='mrna',
                from.state=='SC', to.state=='DE')

highDegree <- network %>% count(target, sort=TRUE) %>% filter(n > 5)

network <- network %>% 
  filter(target %in% highDegree$target)

genes <- network %>%
  select(feature, target) %>% 
  tidyr::gather(source, gene, feature, target) %>% 
  select(gene) %>% 
  mutate(id=gene, label=gene) %>% 
  distinct() %>% 
  arrange(gene)

genes <- genes %>% 
  mutate(vid=sample(vids, nrow(genes), replace = T),
         votes=sample(1:100, nrow(genes), replace=T),
         votesColor=cut(votes, breaks=c(0, 33, 66, 100), 
                   labels=c("red", "yellow", "green")))
