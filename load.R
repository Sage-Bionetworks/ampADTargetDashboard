
if (usePublic) {
  # Public data
  targetListOrigId <- "syn8656625"
  targetManifestId <- "syn11421445"
} else {
  # Unreleased data
  scoreDataId <- "syn11688680"
  targetListOrigId <- "syn11421406"
  targetListDistinctId <- "syn11421426"
  targetManifestId <- "syn11318664"
}

geneExprDataId <- "syn11326321"
fGeneFPKMLongId <- 'syn11327106'
IMSRId <- "syn11420915"

scoreData <- synGet(scoreDataId) %>% 
  getFileLocation() %>% 
  read_csv() %>% 
  rename(ensembl.gene=gene, Score=adDriverScore, Gene=external_gene_name)

targetListOrig <- synGet(targetListOrigId) %>% 
  getFileLocation %>% 
  read_csv()

targetManifest <- synGet(targetManifestId) %>% 
  getFileLocation %>% 
  read_csv()

targetManifestTable <- targetManifest %>% 
  DT::datatable(options=list(lengthChange=FALSE,
                             autoWidth=TRUE, scrollX=FALSE,
                             pageLength=50, dom="ftp"),
                rownames = FALSE,
                selection = list(mode='single', target='row'))

geneExprData <- synGet(geneExprDataId) %>% 
  getFileLocation %>% 
  read_feather()

geneDF <- geneExprData %>%
  dplyr::select(Gene=hgnc_symbol, `ensembl.gene`=ensembl_gene_id) %>%
  dplyr::distinct()

# create list of filters
dForFilter <- geneExprData %>% 
  dplyr::distinct(Study, Tissue, Model, Sex)

dForFilter <- dForFilter %>% 
  tidyr::unite("tissue_study", Tissue, Study, sep=", ", remove=FALSE) %>% 
  tidyr::unite("model_sex", Model, Sex, sep=", ", remove=FALSE) %>% 
  mutate(tissue_study_pretty=glue::glue("{tissue} ({study})", tissue=Tissue, study=Study),
         model_sex_pretty=glue::glue("{model} ({sex})", model=Model, sex=Sex)) 

tissueStudySelectionsDF <- dForFilter %>% 
  select(tissue_study, tissue_study_pretty) %>% 
  distinct()

tissueStudySelections <- purrr::set_names(tissueStudySelectionsDF$tissue_study, 
                                          tissueStudySelectionsDF$tissue_study_pretty)


modelSexSelectionsDF <- dForFilter %>% 
  select(model_sex, model_sex_pretty) %>% 
  distinct()

modelSexSelections <- purrr::set_names(modelSexSelectionsDF$model_sex, 
                                       modelSexSelectionsDF$model_sex_pretty)

IMSR <- synGet(IMSRId) %>% 
  getFileLocation() %>% 
  read_feather()

geneFPKMLong <- synGet(fGeneFPKMLongId) %>% 
  getFileLocation() %>% 
  read_feather()

network <- readr::read_csv(synGet("syn11685347") %>% getFileLocation())

network2 <- network %>% 
  group_by(geneA_ensembl_gene_id, geneB_ensembl_gene_id,
                                 geneA_external_gene_name, geneB_external_gene_name) %>% 
  summarize(value=n_distinct(brainRegion)) %>% 
  #  regions=paste(unique(brainRegion), collapse=","))
  ungroup() %>% 
  distinct() %>% 
  filter(value > 1)

network2 <- network2 %>% left_join(network) %>% 
  group_by(geneA_ensembl_gene_id, geneB_ensembl_gene_id,
           geneA_external_gene_name, geneB_external_gene_name, value) %>% 
  summarize(title=paste(unique(brainRegion), collapse=",")) %>% 
  ungroup()
  
  
nodesForNetwork <- dplyr::bind_rows(network2 %>% 
                                      select(gene=geneA_ensembl_gene_id, 
                                             symbol=geneA_external_gene_name),
                                    network2 %>% 
                                      select(gene=geneB_ensembl_gene_id, 
                                             symbol=geneB_external_gene_name)) %>% 
  distinct() %>% 
  mutate(label=ifelse(is.na(symbol), gene, symbol)) %>% 
  select(id=gene, label) %>% 
  left_join(scoreData %>% select(id=ensembl.gene, Score)) %>% 
  mutate(title=sprintf("score: %0.1f", Score)) %>% 
  mutate(color=cut(Score, breaks=c(-Inf, 0, 2, 4, Inf), 
                   labels=c("red", "yellow", "orange", "green")))

edgesForNetwork <- network2 %>% 
  select(from=geneA_ensembl_gene_id, to=geneB_ensembl_gene_id, value, title)


gg <- graph_from_data_frame(network2)
