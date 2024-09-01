library(DESeq2)
library(data.table)
library(stringr)
library(parallelDist)
library(dplyr)

meta = as.data.frame(fread("meta.csv.gz"))

fs = Sys.glob("fordb/degs/*.rds")
names(fs) = sapply(fs, function(f) { str_split(basename(f), "\\.")[[1]][1] })

meta = meta[meta$perturb_id %in% names(fs), ]

lapply(unique(meta$sci_name), function(sci_name) {
    res = pbapply::pblapply(meta$perturb_id[meta$sci_name == sci_name], function(perturb) {
        dds = readRDS(fs[[perturb]])
        temp = as.data.frame(results(dds))
        temp$gene = rownames(temp)
        temp = temp[, c("gene", "log2FoldChange")]
        temp$perturb = perturb
        temp
    }, cl = 10)
    res = do.call(rbind, res)
    res = reshape2::dcast(res, gene~perturb, value.var = "log2FoldChange", fill = 0, fun.aggregate = mean)
    rownames(res) = res$gene
    res = res[, colnames(res) != "gene"]
    trqwe::mcsaveRDS(res, glue::glue("{sci_name}.rds"))
    rm(res)
    gc()
})


res = lapply(Sys.glob("expr/*.rds"), readRDS)
names(res) = lapply(Sys.glob("expr/*.rds"), basename)
names(res) = str_replace_all(names(res), ".rds", "")

dists = pbapply::pblapply(res, function(x) { as.matrix(parDist(t(x), method = "euclidean", threads = 100)) })

format_dists = pbapply::pblapply(names(dists), function(x) {
    temp = reshape2::melt(dists[[x]])
    temp$sci_name = x
    temp[temp$Var1 != temp$Var2, ]
})
format_dists = do.call(rbind, format_dists)
colnames(format_dists)[1:2] = c("source", "target")

format_dists = format_dists[order(format_dists$source, format_dists$value), ]
write.csv(format_dists, "all_dists.csv", row.names = F, quote = F)


format_dists1 = format_dists %>%
    group_by(source) %>%
    top_n(50, wt = -value) %>%
    as.data.frame()

format_dists1 = format_dists1[order(format_dists1$source, format_dists1$value), ]
write.csv(format_dists1, "format_dists.csv", row.names = F, quote = F)

