druggabilityData <- fread(getFileLocation(synGet("syn7555804")), 
                          data.table=FALSE) %>% 
  rename(ensembl.gene=ENSG) %>% 
  mutate(status_assays=fct_recode(status_assays, unknown=""),
         status_crystal_structure=fct_recode(status_crystal_structure, unknown=""),
         status_pocket=fct_recode(status_pocket, unknown=""),
         status_in_vivo_work=fct_recode(status_in_vivo_work, unknown=""),
         status_known_ligands=fct_recode(status_known_ligands, unknown=""),
         Lilly_DrugEBIlity_Consensus_Score=fct_explicit_na(as.character(Lilly_DrugEBIlity_Consensus_Score), 
                                                           na_level = "unk"),
         `Lilly_GW_Druggability_Structure-based`=fct_explicit_na(as.character(`Lilly_GW_Druggability_Structure-based`), 
                                                                 na_level = "unk")
  )

druggabilityData <- druggabilityData %>% 
  select(GENE_SYMBOL, starts_with('status')) %>% 
  distinct() %>% 
  tidyr::gather(category, status, starts_with('status')) %>% 
  mutate(status_numeric=fct_recode(status, `0`="unknown", `0`='bad',
                                   `1`="medium", `2`="good"),
         status_numeric=levels(status_numeric)[as.numeric(status_numeric)]) %>% 
  select(GENE_SYMBOL, status_numeric) %>% 
  group_by(GENE_SYMBOL) %>% 
  summarize(sum_status=sum(as.numeric(status_numeric))) %>% 
  ungroup() %>% 
  right_join(druggabilityData, by=c('GENE_SYMBOL')) %>% 
  arrange(GENE_SYMBOL) %>%
  select(Gene=GENE_SYMBOL,
         ensembl.gene,
         `ODDI Druggability Score`=sum_status,
         `Lilly DrugEBIlity Consensus`=Lilly_DrugEBIlity_Consensus_Score) %>% 
  distinct()


targetList <- synGet("syn8656625") %>% getFileLocation %>% fread(data.table=FALSE) %>% 
  select(Center=group, Gene=gene_symbol, ensembl.gene=ensembl_id)

targetManifest <- targetList %>% left_join(druggabilityData)

network <- fread(getFileLocation(synGet("syn7770770")), 
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

# geneFPKM <- fread(getFileLocation(synGet("syn5581268")), 
#                   data.table=FALSE) %>% 
#   filter(ensembl_gene_id %in% c(genesForNetwork$gene, 
#                                 druggabilityData$ensembl.gene)) 
# 
# geneCovariates <- fread(getFileLocation(synGet("syn5581227")),
#                         data.table=FALSE) %>% 
#   filter(cogdx %in% c(1, 4)) %>%
#   mutate(cogdx=factor(cogdx, ordered=TRUE))
# 
# geneFPKMLong <- geneFPKM %>% 
#   tidyr::gather(sample, fpkm, 3:ncol(geneFPKM)) %>% 
#   left_join(geneCovariates %>% select(Sampleid_batch, cogdx), 
#             by=c("sample"="Sampleid_batch")) %>% 
#   filter(!is.na(cogdx))

geneFPKMLong <- fread(getFileLocation(synGet("syn7555798")), 
                      data.table=FALSE)

geneFPKMLong$cogdx <- forcats::fct_recode(geneFPKMLong$cogdx, 
                                          NCI='1', AD='4')

gtexObj <- synGet('syn7542283')

gtex <- fread(getFileLocation(gtexObj), data.table=FALSE) %>% 
  mutate(ensembl.gene=str_replace(Name, "\\..*", "")) %>% 
  dplyr::filter(ensembl.gene %in% c(genesForNetwork$gene, 
                                    druggabilityData$ensembl.gene)) %>% 
  select(ensembl.gene, hgnc_symbol=Description, starts_with('Brain'))

gtex <- gtex %>% 
  tidyr::gather(tissue, medianFPKM, 3:ncol(gtex)) %>% 
  mutate(tissue=str_replace(tissue, "Brain - ", ""))

medianGTEx <- median(gtex$medianFPKM)