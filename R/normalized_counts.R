library(DESeq2)
library(reshape2)
library(pbapply)
library(stringr)
library(data.table)
library(glue)

outdir = "degs/normalized"
dir.create(outdir, showWarnings = F)

pblapply(Sys.glob("degs/*.rds"), function(path) {
    dds = readRDS(path)
    key = str_replace_all(basename(path), "\\.rds", "")

    temp = melt(counts(dds, normalized=T))
    temp$perturb_id = key
    temp = temp[temp$value > 0, ]

    lapply(unique(temp$Var2), function(x){
        if(!file.exists(glue("{outdir}/{x}.csv.gz"))) {
            fwrite(temp[temp$Var2 == x, ], glue("{outdir}/{x}.csv.gz"), col.names = F, row.names=F, quote=F)
        }
    })
}, cl = 10)