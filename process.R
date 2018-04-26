library(tidyverse)
library(synapser)
library(feather)

synLogin()

wotFolderId <- 'syn7525089'
scoreDataId <- 'syn11688680'
targetListOrigId <- "syn8656625"

targetListOutputFile <- "./targetList.csv"
targetListDistinctOutputFile <- "./targetListDistinct.csv"
targetManifestOutputFile <- "./targetManifest.csv"

scoreData <- synGet(scoreDataId)$path %>% 
  read_csv() %>% 
  rename(ensembl.gene=gene, Score=adDriverScore, Gene=external_gene_name)

####### Process target list
targetListOrig <- synGet(targetListOrigId)$path %>% 
  read_csv()

targetList <- targetListOrig %>%
  select(Center=group, Gene=gene_symbol, ensembl.gene=ensembl_id)

write_csv(targetList, targetListOutputFile)

fTargetList <- synStore(File(targetListOutputFile, 
                             parent=wotFolderId), 
                        used=targetListOrigId, forceVersion=FALSE)

targetListDistinct <- targetList %>%
  group_by(Gene, ensembl.gene) %>%
  mutate(Centers=paste(Center, collapse=", "),
         nominations=length(unique(Center))) %>%
  ungroup() %>%
  select(-Center) %>%
  distinct() %>%
  arrange(-nominations)

write_csv(targetList, targetListDistinctOutputFile)

fTargetListDistinct <- synStore(File(targetListDistinctOutputFile, 
                                     parent=wotFolderId), 
                        used=fTargetList, forceVersion=FALSE)

targetManifest <- scoreData %>%
  filter(Gene %in% targetListDistinct$Gene) %>% 
  select(Gene, Score) %>% 
  mutate(Score=round(Score, 1)) %>% 
  arrange(-Score)

write_csv(targetManifest, targetManifestOutputFile)

fTargetManifest <- synStore(File(targetManifestOutputFile, 
                                     parent=wotFolderId), 
                                used=c(fTargetListDistinct, scoreDataId), 
                            forceVersion=FALSE)
