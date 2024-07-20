
library(data.table)
library(clusterProfiler)

rds_file <-
  list.files(path = 'enrichment/',
             pattern = '*.Rds',
             full.names = T)

go_inf <- readRDS('OrySat/release52/GO_inf.Rds')

parallel::mclapply(rds_file, function(x) {
  alias_id <- stringr::str_replace_all(basename(x), '[.]Rds', '')
  tmp <- readRDS(x)
  
  go_up_res <- data.frame(tmp$go_up_res)
  
  if (dim(go_up_res)[1] != 0) {
    go_up_res$cate <- 'UP'
    if (!('ONTOLOGY' %in% colnames(go_up_res))) {
      go_up_res <- cbind(data.frame(ONTOLOGY = NA), go_up_res)
    }
    
    if (any(is.na(go_up_res$ONTOLOGY))) {
      go_up_res$ONTOLOGY <- plyr::mapvalues(
        from = go_inf$GOID,
        to = go_inf$ONTOLOGY,
        x = go_up_res$ID,
        warn_missing = F
      )
      go_up_res <- go_up_res[!is.na(go_up_res$Description), ]
    }
    
  }
  go_down_res <- data.frame(tmp$go_down_res)
  if (dim(go_down_res)[1] != 0) {
    go_down_res$cate <- 'DOWN'
    if (!('ONTOLOGY' %in% colnames(go_down_res))) {
      go_down_res <- cbind(data.frame(ONTOLOGY = NA), go_down_res)
    }
    
    if (any(is.na(go_down_res$ONTOLOGY))) {
      go_down_res$ONTOLOGY <- plyr::mapvalues(
        from = go_inf$GOID,
        to = go_inf$ONTOLOGY,
        x = go_down_res$ID,
        warn_missing = F
      )
      go_down_res <- go_down_res[!is.na(go_down_res$Description), ]
      
    }
    
  }
  
  data.table::fwrite(
    rbind(go_up_res, go_down_res),
    file = glue::glue('enrichment_file/{alias_id}.tsv.gz'),
    sep = '\t'
  )
  
  
}, mc.cores = 48)