library(data.table)
library(pbapply)
library(reshape2)
library(stringr)
library(glue)
library(rjson)

inpath = "salmon_de_test/output"
outpath = "salmon_de.csv.gz"

read_single_file <- function(path) {
    df = fread(path)
    df = as.data.frame(df)
    df = reshape2::melt(df, id.vars = c("feature_id", "gene_id", "pvalue", "log2fcp"))
    colnames(df)[(ncol(df) - 1): ncol(df)] = c("experiment_accession", "delta")

    df$perturb_id = str_split(basename(path), "\\.")[[1]][1]

    for(i in c("feature_id", "gene_id")) {
        df[str_detect(df[, i], "^ENS"), i] = sapply(df[str_detect(df[, i], "^ENS"), i], function(x) {
            str_split(x, "\\.")[[1]][1]
        })
    }

    for (i in c("pvalue", "log2fcp", "delta")) {
        df[,i ] = round(df[,i], digits=20)
    }
    na.omit(df)
}


fs = Sys.glob(glue("{inpath}/*.gz"))

res = pblapply(fs, read_single_file, cl = 10)
res = do.call(rbind, res)

fwrite(res, outpath)