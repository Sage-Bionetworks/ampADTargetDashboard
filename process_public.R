library(tidyverse)
library(synapseClient)
library(feather)

synapseLogin()

wotFolderId <- 'syn7525089'
druggabilityDataId <- 'syn11420932'
targetListOrigId <- "syn11420934"

geneExprDataId <- 'syn11180450'
IMSRId <- 'syn11149859'
geneFPKMId <- "syn5581268"
geneCovariatesId <- "syn5581227"

druggabilityDataOutputFile <- "./druggabilityData_IGAP.csv"
targetListOutputFile <- "./targetList_IGAP.csv"
targetListDistinctOutputFile <- "./targetListDistinct_IGAP.csv"
targetManifestOutputFile <- "./targetManifest_GAP.csv"

fGeneExprDataOutputFile <- "./geneExprData.feather"
IMSROutputFile <- "./IMSR_processed.feather"
geneFPKMLongOutputFile <- "./geneFPKMLong.feather"

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

### Process gene expression (logfc and CI) data
geneExprData <- synGet(geneExprDataId) %>%
  getFileLocation() %>%
  read_tsv() %>%
  filter(stringr::str_detect(Model, "Diagnosis")) %>%
  rename(Sex=Gender) %>%
  mutate(hgnc_symbol=ifelse(is.na(hgnc_symbol), ensembl_gene_id, hgnc_symbol),
         Sex=forcats::fct_recode(Sex, Males="MALE",
                                 Females="FEMALE",
                                 `Males and Females`="ALL"),
         Model=stringr::str_replace(Model, "\\.", " x "))

write_feather(geneExprData, fGeneExprDataOutputFile)
fGeneExprData <- synStore(File(fGeneExprDataOutputFile,
                               parentId=wotFolderId),
                          used=c(geneExprDataId),
                          forceVersion=FALSE)

IMSR <- synGet(IMSRId) %>%
  getFileLocation() %>%
  readr::read_csv() %>%
  dplyr::select(`Strain ID`, `Strain/Stock`, Repository, `Gene Symbol`, URL) %>%
  mutate(`Gene Symbol`=toupper(`Gene Symbol`)) %>%
  mutate(`Strain/Stock`=stringr::str_c("<a href='", URL, "'>", `Strain/Stock`, "</a>"))

write_feather(IMSR, IMSROutputFile)
fIMSR <- synStore(File(IMSROutputFile,
                       parentId=wotFolderId),
                  used=c(IMSRId),
                  forceVersion=FALSE)

geneFPKM <- synGet(geneFPKMId) %>%
  getFileLocation() %>%
  read_tsv()

geneCovariates <- synGet(geneCovariatesId) %>%
  getFileLocation() %>%
  read_tsv() %>%
  filter(cogdx %in% c(1, 4)) %>%
  mutate(cogdx=factor(cogdx, ordered=TRUE))

geneFPKMLong <- geneFPKM %>%
  tidyr::gather(sample, fpkm, 3:ncol(geneFPKM)) %>%
  left_join(geneCovariates %>% dplyr::select(Sampleid_batch, cogdx),
            by=c("sample"="Sampleid_batch")) %>%
  dplyr::filter(!is.na(cogdx), !is.na(hgnc_symbol))

geneFPKMLong <- geneFPKMLong %>%
  mutate(cogdx=forcats::fct_recode(geneFPKMLong$cogdx, NCI='1', AD='4')) %>%
  dplyr::rename(`Cognitive Diagnosis`=cogdx)

write_feather(geneFPKMLong, geneFPKMLongOutputFile)
fGeneFPKMLong <- synStore(File(geneFPKMLongOutputFile,
                               parentId=wotFolderId),
                          used=c(geneFPKMId, geneCovariatesId),
                          forceVersion=FALSE)

# network <- fread(getFileLocation(synGet("syn7770770")),
#                  data.table=FALSE)
# 
# genesForNetwork <- network %>%
#   tidyr::gather(source, gene, var1, var2) %>%
#   dplyr::select(gene) %>%
#   mutate(id=gene) %>%
#   distinct() %>%
#   arrange(gene)
# 
# out <- queryMany(unique(genesForNetwork$gene),
#                  scopes="ensembl.gene", fields="symbol", species="human",
#                  returnall=TRUE, size=1)
# 
# res <- as.data.frame(out$response) %>% select(symbol, query)
# 
# genesForNetwork <- genesForNetwork %>%
#   left_join(res, by=c('id'='query')) %>%
#   mutate(label=ifelse(is.na(symbol), gene, symbol))
# 
# gg <- graph_from_data_frame(network %>% distinct())
# 
# gtexObj <- synGet('syn7542283')
# 
# gtex <- fread(getFileLocation(gtexObj), data.table=FALSE) %>%
#   mutate(ensembl.gene=str_replace(Name, "\\..*", "")) %>%
#   # dplyr::filter(ensembl.gene %in% c(genesForNetwork$gene,
#   #                                   targetList$ensembl.gene)) %>%
#   select(ensembl.gene, hgnc_symbol=Description, starts_with('Brain'))
# 
# gtex <- gtex %>%
#   tidyr::gather(tissue, medianFPKM, 3:ncol(gtex)) %>%
#   mutate(tissue=str_replace(tissue, "Brain - ", ""))
# 
# medianGTEx <- median(gtex$medianFPKM)
# # 
# nTargets <- targetListOrig %>% count(group)
# 
# save(ISMR, nTargets, network, targetList, targetListOrig,
#      druggabilityData, targetManifest, targetListDistinct, genesForNetwork,
#      gg, geneFPKMLong, gtex, medianGTEx, geneExprData, geneDF, dForFilter,
#      file="./data2load.RData")
# synStore(File("./data2load.RData", parentId="syn7525089"))
# diffExprData <- synapseClient::synGet("syn11180450") %>% 
#   synapseClient::getFileLocation() %>% 
#   readr::read_csv()
