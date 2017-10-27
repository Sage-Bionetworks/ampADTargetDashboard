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
                             pageLength=50, dom="ftp"),
                rownames = FALSE,
                selection = list(mode='single', target='row')) %>%
  DT::formatStyle(c('Assays', 'Known ligands', 'In vivo'),
                  backgroundColor=DT::styleEqual(c("good", "medium", "bad", "unknown"), #levels(tmp$status),
                                                 c("green", "orange", "red", "grey")))

geneExprDataId <- "syn11318688"
geneExprData <- synGet(geneExprDataId) %>% 
  getFileLocation %>% 
  read_csv()

geneDF <- geneExprData %>%
  dplyr::select(Gene=hgnc_symbol, `ensembl.gene`=ensembl_gene_id) %>%
  dplyr::distinct()

# create list of filters
dForFilter <- geneExprData %>% 
  dplyr::distinct(Study, Tissue, Model, Sex)

IMSRId <- "syn11318727"
IMSR <- synGet(IMSRId) %>% 
  getFileLocation() %>% 
  read_csv()

fGeneFPKMLongId <- 'syn7555798'
geneFPKMLong <- synGet(fGeneFPKMLongId) %>% 
  getFileLocation() %>% 
  read_csv()
