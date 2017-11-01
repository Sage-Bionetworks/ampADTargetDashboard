druggabilityDataId <- "syn11318563"
druggabilityData <- synGet(druggabilityDataId) %>% 
  getFileLocation() %>% read_csv()

targetListOrigId <- "syn8656625"
targetListOrig <- synGet(targetListOrigId) %>% 
  getFileLocation %>% 
  read_csv()

targetListDistinctId <- "syn11318663"
targetListDistinct <- synGet(targetListDistinctId) %>% 
  getFileLocation %>% 
  read_csv()

targetManifestId <- "syn11318664"
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


geneExprDataId <- "syn11326321"
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

IMSRId <- "syn11318727"
IMSR <- synGet(IMSRId) %>% 
  getFileLocation() %>% 
  read_csv()

fGeneFPKMLongId <- 'syn11327106'
geneFPKMLong <- synGet(fGeneFPKMLongId) %>% 
  getFileLocation() %>% 
  read_feather()
