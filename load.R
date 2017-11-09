# Unreleased data
# druggabilityDataId <- "syn11420935"
# targetListOrigId <- "syn11421406"
# targetListDistinctId <- "syn11421426"
# targetManifestId <- "syn11318664"

# Public data
druggabilityDataId <- "syn11420935"
targetListOrigId <- "syn8656625"
targetListDistinctId <- "syn11318663"
targetManifestId <- "syn11421445"

geneExprDataId <- "syn11326321"
fGeneFPKMLongId <- 'syn11327106'
IMSRId <- "syn11420915"

druggabilityData <- synGet(druggabilityDataId) %>% 
  getFileLocation() %>% read_csv()

targetListOrig <- synGet(targetListOrigId) %>% 
  getFileLocation %>% 
  read_csv()

targetListDistinct <- synGet(targetListDistinctId) %>% 
  getFileLocation %>% 
  read_csv()

targetManifest <- synGet(targetManifestId) %>% 
  getFileLocation %>% 
  read_csv()

targetManifestTable <- targetManifest %>%
  left_join(druggabilityData %>% dplyr::select(Gene=GENE_SYMBOL,
                                               Assays=status_assays,
                                               `In vivo`=status_in_vivo_work,
                                               `Known ligands`=status_known_ligands)) %>%
  DT::datatable(options=list(lengthChange=FALSE,
                             autoWidth=TRUE, scrollX=TRUE,
                             pageLength=50, dom="ftp",
                             columnDefs = list(list(width = '75px', targets =c(0,3,4,5)))),
                rownames = FALSE,
                selection = list(mode='single', target='row')) %>%
  DT::formatStyle(c('Assays', 'Known ligands', 'In vivo'),
                  backgroundColor=DT::styleEqual(c("good", "medium", "bad", "unknown"), #levels(tmp$status),
                                                 c("green", "orange", "red", "grey")),
                  color=DT::styleEqual(c("good", "medium", "bad", "unknown"), #levels(tmp$status),
                                       c("green", "orange", "red", "grey"))) %>% 
  DT::formatStyle(columns = c(1, 2, 3), fontSize = '125%')

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
