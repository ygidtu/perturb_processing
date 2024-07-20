
library(data.table)
library(dplyr)
library(clusterProfiler)
library(pbapply)
library(glue)
library(stringr)

annot_inf <- fread('03.generate_pair/annotation.csv')
for (i in annot_inf$loc) {
    if (!file.exists(i)) {
        if (!is.element(i, .packages(all.available = TRUE)) ) {
            BiocManager::install(i)
        }
        library(i, character.only = TRUE)
    }
}
for (i in annot_inf$loc) {
    if (!str_starts(i, "/")) {
        if (!is.element(i, .packages(all.available = TRUE)) ) {
            BiocManager::install(i)
        }
        library(i, character.only = TRUE)
    }
}


pair_meta <- as.data.frame(fread('meta.csv.gz'))
rownames(pair_meta) = pair_meta$Perturb_ID

output_dir = "enrichments"
dir.create(output_dir, showWarnings = F)

gene_map_ss <- readRDS('Ref/SusScr/release101/gene_map.Rds')
gene_map_ec <- readRDS('Ref/EscCol/release52/gene_map.Rds')
go_inf <- readRDS('Ref/OrySat/release52/GO_inf.Rds')


pbapply(pair_meta, 1, function(row) {
    perturb_id = row["perturb_id"]
    f = glue("degs/{perturb_id}.csv.gz")

    tryCatch({
        de_res <- fread(f)
        colnames(de_res) <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "gene_id", "perturb_id")
        org_use <- annot_inf$loc[annot_inf$org == str_trim(row["org"])]
        if (org_use == "org.Ss.eg.db") {
          de_res$gene_id <- plyr::mapvalues(
            from = gene_map_ss$gene_id,
            to = gene_map_ss$gene_name,
            x = de_res$gene_id,
            warn_missing = F
          )
        }
        if (org_use == "org.EcK12.eg.db") {
          de_res$gene_id <- plyr::mapvalues(
            from = gene_map_ec$gene_id,
            to = gene_map_ec$gene_name,
            x = de_res$gene_id,
            warn_missing = F
          )
        }

        up_gene <- subset(de_res, pvalue < 0.05 & log2FoldChange >= log2(1.5)) %>% pull(gene_id)
        down_gene <- subset(de_res, pvalue < 0.05 & log2FoldChange <= -log2(1.5)) %>% pull(gene_id)
        key_type <- annot_inf$key[annot_inf$org == str_trim(row["org"])]
        if (endsWith(org_use, 'db')) {
          if (length(up_gene) == 0) {
            go_up_res = NULL
          } else {
            tryCatch({
                go_up_res <- enrichGO(
                  gene = up_gene,
                  OrgDb = org_use,
                  qvalueCutoff = 1,
                  pvalueCutoff = 1,
                  keyType = key_type,
                  ont = 'ALL'
                )
            }, error = function(e){ go_up_res = NULL })
          }
          if (length(down_gene) == 0) {
            go_down_res = NULL
          } else {
            tryCatch({
                go_down_res <- enrichGO(
                  gene = down_gene,
                  OrgDb = org_use,
                  qvalueCutoff = 1,
                  pvalueCutoff = 1,
                  keyType = key_type,
                  ont = 'ALL'
                )
            }, error = function(e) {
                go_down_res = NULL
            })
          }
        } else {
          go_up_res = enricher(
            up_gene,
            TERM2GENE = readRDS(org_use),
            TERM2NAME = go_inf[, c(1, 2)],
            pvalueCutoff = 1,
            qvalueCutoff = 1
          )

          go_down_res = enricher(
            down_gene,
            TERM2GENE = readRDS(org_use),
            TERM2NAME = go_inf[, c(1, 2)],
            pvalueCutoff = 1,
            qvalueCutoff = 1
          )
        }

        saveRDS(list(go_up_res = go_up_res, go_down_res = go_down_res), glue("{output_dir}/{perturb_id}.rds"))


        go_up_res <- data.frame(go_up_res)

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
            go_up_res <- go_up_res[!is.na(go_up_res$Description),]
          }

        }
        go_down_res <- data.frame(go_down_res)
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
            go_down_res <- go_down_res[!is.na(go_down_res$Description),]

          }

        }

        fwrite(
          rbind(go_up_res, go_down_res),
          file = glue("{output_dir}/{perturb_id}.tsv.gz"),
          sep = '\t'
        )

    }, error = function(e) {
        print(e)
    })
}, cl = 10)
