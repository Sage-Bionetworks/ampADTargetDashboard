library(tidyverse)
library(synapseClient)
library(feather)

synapseLogin()

wotFolderId <- 'syn7525089'
druggabilityDataId <- 'syn7555804'
targetListOrigId <- "syn8656625"


druggabilityDataOutputFile <- "./druggabilityData.csv"
targetListOutputFile <- "./targetList.csv"
targetListDistinctOutputFile <- "./targetListDistinct.csv"
targetManifestOutputFile <- "./targetManifest.csv"

druggabilityData <- synGet(druggabilityDataId) %>% 
  getFileLocation() %>% 
  read_csv() %>%
  select(-Center, -DDI_Interest, -some_number, -GENE_DESCRIPTION) %>%
  rename(ensembl.gene=ENSG) %>%
  mutate(status_assays=forcats::fct_recode(status_assays, unknown=""),
         status_crystal_structure=forcats::fct_recode(status_crystal_structure, unknown=""),
         status_pocket=forcats::fct_recode(status_pocket, unknown=""),
         status_in_vivo_work=forcats::fct_recode(status_in_vivo_work, unknown=""),
         status_known_ligands=forcats::fct_recode(status_known_ligands, unknown=""),
         Lilly_DrugEBIlity_Consensus_Score=forcats::fct_explicit_na(as.character(Lilly_DrugEBIlity_Consensus_Score),
                                                           na_level = "unk"),
         `Lilly_GW_Druggability_Structure-based`=forcats::fct_explicit_na(as.character(`Lilly_GW_Druggability_Structure-based`),
                                                                 na_level = "unk")) %>% 
  group_by(GENE_SYMBOL) %>% 
  slice(1) %>% 
  ungroup()

druggabilityData <- druggabilityData %>%
  select(GENE_SYMBOL, starts_with('status')) %>%
  distinct() %>%
  tidyr::gather(category, status, starts_with('status')) %>%
  mutate(status_numeric=forcats::fct_recode(status, `0`="unknown", `0`='bad',
                                   `1`="medium", `2`="good"),
         status_numeric=levels(status_numeric)[as.numeric(status_numeric)]) %>%
  select(GENE_SYMBOL, status_numeric) %>%
  group_by(GENE_SYMBOL) %>%
  summarize(sum_status=sum(as.numeric(status_numeric))) %>%
  ungroup() %>%
  right_join(druggabilityData, by=c('GENE_SYMBOL')) %>%
  arrange(GENE_SYMBOL)

write_csv(druggabilityData, druggabilityDataOutputFile)

fDrugData <- synStore(File(druggabilityDataOutputFile, 
                           parentId=wotFolderId), 
                      used=druggabilityDataId, forceVersion=FALSE)

####### Process target list
targetListOrig <- synGet(targetListOrigId) %>% 
  getFileLocation %>% 
  read_csv()

targetList <- targetListOrig %>%
  select(Center=group, Gene=gene_symbol, ensembl.gene=ensembl_id)

write_csv(targetList, targetListOutputFile)

fTargetList <- synStore(File(targetListOutputFile, 
                             parentId=wotFolderId), 
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
                             parentId=wotFolderId), 
                        used=fTargetList, forceVersion=FALSE)

targetManifest <- targetListDistinct %>%
  left_join(druggabilityData, by=c('Gene'='GENE_SYMBOL',
                                   'ensembl.gene'='ensembl.gene')) %>%
  select(Gene,
         Centers,
         nominations) %>%
  distinct() %>%
  arrange(-nominations)

write_csv(targetManifest, targetManifestOutputFile)

fTargetManifest <- synStore(File(targetManifestOutputFile, 
                                     parentId=wotFolderId), 
                                used=c(fTargetListDistinct, fDrugData), 
                            forceVersion=FALSE)
