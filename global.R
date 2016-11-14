library(synapseClient)
library(data.table)
library(dplyr)
library(mygene)
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


out <- queryMany(unique(c(ddiData$GENE_SYMBOL, lillyData$GENE_SYMBOL)),
                        scopes="symbol", fields="ensembl.gene", species="human",
                 returnall=TRUE, size=1)

res <- as.data.frame(out$response)

# Two genes - CHRH1 (in DDI and Lilly) and MOAP1 (in Lilly) - have multiple
# matches to some other chromosomes in Ensembl. Get only the first one.
res <- res %>% mutate(ensembl=laply(ensembl, function(x) ifelse(is.null(x), '', x[[1]][1])))
res[is.na(res$ensembl.gene), 'ensembl.gene'] <- res[is.na(res$ensembl.gene), 'ensembl']
res <- res %>% select(-ensembl, -X_id)

ddiData <- ddiData %>% left_join(res, by=c('GENE_SYMBOL'='query'))
lillyData <- lillyData %>% left_join(res, by=c('GENE_SYMBOL'='query'))

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
