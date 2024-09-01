library(DESeq2)
library(data.table)
library(stringr)
library(pbapply)
library(glue)


meta = as.data.frame(fread("meta.csv.gz"))

geneexpr = "geneexpr"
cl = 30
output = "degs"
dir.create(output, showWarnings = F)


read_expr <- function(path) {
    expr = NULL
    for(x in list.files(path, full.names = T)) {
        temp = as.data.frame(fread(x))
        temp = temp[, 2:3]
        colnames(temp) = c("gene", basename(x))
        
        if (is.null(expr)) {
            expr = temp
        } else {
            expr = merge(expr, temp, by = "gene")
        }
    }

    rownames(expr) <- expr$gene
    return(expr[, colnames(expr) != "gene"])
}


pbapply(meta, 1, function(row) {
    perturb_id = row[["perturb_id"]]
    perturbs = str_split(row[["perturbs"]], ",")[[1]]
    controls = str_split(row[["controls"]], ",")[[1]]
    
    tryCatch({
        o = glue("{output}/{perturb_id}.rds")
        if (file.exists(o)) {
            dds = readRDS(o)
        } else {
            meta = data.frame(
                samples=c(perturbs, controls),
                condition=c(rep("Perturb", length(perturbs)), rep("Control", length(controls))),
                row.names = c(perturbs, controls)
            )
            expr = read_expr(glue("{geneexpr}/{perturb_id}"))

            dds = DESeqDataSetFromMatrix(as.matrix(expr), meta[colnames(expr), ], design = ~condition)


            dds = DESeq(dds, test = "Wald")
            saveRDS(dds, glue("{output}/{perturb_id}.rds"))
        }
        res = as.data.frame(results(dds, contrast = c('condition', 'Perturb', 'Control')))
        res$gene = rownames(res)
        res$perturb_id = perturb_id
        fwrite(res, glue("{output}/{perturb_id}.csv.gz"), col.names = F)
        return(NULL)
    }, error = function(e) {
        return(NULL)
    })
}, cl = cl)


