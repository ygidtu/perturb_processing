argvs <- commandArgs(T)

meta <- readRDS(argvs[2])
pert_label <- unique(meta$perturb_id)
output_file <- file.path('output', glue::glue("{pert_label}.txt.gz"))

if (file.exists(output_file)) {
  cat("File exists. Exiting the program.\n")
  quit(save = "no", status = 0)
}

library(data.table)
library(tximport)
library(DRIMSeq)
library(stageR)
library(dplyr)
library(reshape2)
ref <- fread(argvs[1])



## Load annotation
txdb.filename = ref$loc[ref$ID %in% unique(meta$sci_name)]
txdf = data.table::fread(txdb.filename,header=F)
colnames(txdf) <- c("TXNAME", "GENEID")




colnames(meta) <- c('perturb_id',
                    'sci_name',
                    'sample_id',
                    'group')

fileloc <- readRDS('salmon_de_test/file_loc.Rds')
fileloc <-
  fileloc[fileloc$experiment_accession %in% meta$sample_id,]

if (nrow(fileloc) != nrow(meta)) {
  cat("Sample not equal, quit!\n")
  quit(save = "no", status = 0)
  
}

fileloc$experiment_accession <-
  factor(fileloc$experiment_accession, levels = meta$sample_id)
fileloc <- fileloc[order(fileloc$experiment_accession), ]


files <- file.path(fileloc$path, 'quant.sf')
names(files) <- meta$sample_id

index_use <- file.exists(files)
files <- files[index_use]
meta <- meta[index_use,]

if (isTRUE(any(table(meta$group) < 2))) {
  cat("Sample not equal, quit!\n")
  quit(save = "no", status = 0)
  
}

txi <- tximport(
  files,
  type = "salmon",
  txOut = TRUE,
  countsFromAbundance = "scaledTPM"
)

cts <- txi$counts

cts <- cts[rowSums(cts) > 0,]


txdf.sub = txdf[match(rownames(cts),txdf$TXNAME),]
counts = data.frame(gene_id = txdf.sub$GENEID, feature_id = txdf.sub$TXNAME, cts)
counts <- counts[counts$gene_id %in% counts$gene_id[duplicated(counts$gene_id)], ]


d = dmDSdata(counts = counts, samples = meta)

n = nrow(meta)
n.small = min(table(meta$group))

d = dmFilter(
  d,
  min_samps_feature_expr = n.small,
  min_feature_expr = 10,
  min_samps_feature_prop = n.small,
  min_feature_prop = 0.1,
  min_samps_gene_expr = n,
  min_gene_expr = 10
)

design = model.matrix(~group, data = meta)

set.seed(12345)
system.time({
  d <- dmPrecision(d, design = design)    # Estimate the precision (Higher dispersion is associated with lower precision)
  d <- dmFit(d, design = design)          # Fit regression coefficients
  d <- dmTest(d, coef = "grouppert")     # Perform null hypothesis testing on the coefficient of interest
})

# saveRDS(d, file = 'd.Rds')

# Single p-value per gene
res.g = DRIMSeq::results(d)

# Single p-value per transcript
res.t = DRIMSeq::results(d, level = "feature")

no.na <- function(x) ifelse(is.na(x), 1, x)
res.g$pvalue <- no.na(res.g$pvalue)
res.t$pvalue <- no.na(res.t$pvalue)


smallProportionSD <- function(d, filter = 0.1) {
  # Generate count table
  cts = as.matrix(subset(counts(d), select = -c(gene_id, feature_id)))
  # cts[is.na(cts)] <- 0
  # Summarise count total per gene
  gene.cts = rowsum(cts, counts(d)$gene_id)
  # Use total count per gene as count per transcript
  total.cts = gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
  # Calculate proportion of transcript in gene
  props = cts / total.cts
  # props[is.na(props)] <- 0
  rownames(props) = rownames(total.cts)
  
  # Calculate standard deviation
  propSD = sqrt(rowVars(props))
  # Check if standard deviation of per-sample proportions is < 0.1
  propSD < filter
}


filt = smallProportionSD(d)

res.t.filt = DRIMSeq::results(d, level = "feature")
res.t.filt$pvalue[filt] = 1
res.t.filt$adj_pvalue[filt] = 1

strp <- function(x) substr(x,1,15)

pScreen = res.g$pvalue
names(pScreen) = strp(res.g$gene_id)

pScreen <- pScreen[pScreen < 0.05]
pScreen <- pScreen[names(pScreen) %in% names(pScreen)[duplicated(names(pScreen))]]

tmp_p <- res.t.filt$pvalue
tmp_p[is.na(tmp_p)] <- 1
pConfirmation = matrix(tmp_p, ncol = 1)
dimnames(pConfirmation) = list(strp(res.t.filt$feature_id), "transcript")


tx2gene = data.frame(res.t[, c("feature_id", "gene_id")],
                     res.t[, c("feature_id", "gene_id")])

for (i in 1:2) tx2gene[, i] = strp(tx2gene[, i])
# tx2gene_rm <-
#   tx2gene[tx2gene$gene_id %in% tx2gene$gene_id[duplicated(tx2gene$gene_id)],]
# 
# 
# 
# pConfirmation_new <- matrix(pConfirmation[tx2gene_rm$feature_id,], ncol=1)
# dimnames(pConfirmation_new) = list(strp(tx2gene_rm$feature_id), "transcript")
# 
# (names(pScreen) %in% tx2gene_rm$gene_id)


tx2gene_use <- tx2gene[tx2gene$gene_id %in% names(pScreen), ]
tmp_p <- pConfirmation[rownames(pConfirmation) %in% tx2gene_use$feature_id,]
pConfirmation = matrix(tmp_p, ncol = 1)
dimnames(pConfirmation) = list(names(tmp_p), "transcript")

stageRObj = stageRTx(
  pScreen = pScreen,
  pConfirmation = pConfirmation,
  pScreenAdjusted = FALSE,
  tx2gene = tx2gene_use[,c(1, 2)]
)

stageRObj = stageWiseAdjustment(stageRObj,
                                method = "dtu",
                                alpha = 0.05,
                                allowNA = TRUE)

drim.padj = unique(getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = TRUE)[, 1:3])
drim.padj$feature_id.1 <- plyr::mapvalues(
  from = tx2gene_use$feature_id,
  to = tx2gene_use$feature_id.1,
  x = drim.padj$txID,
  warn_missing = F
)
drim.padj$gene_id.1 <- plyr::mapvalues(
  from = tx2gene_use$gene_id,
  to = tx2gene_use$gene_id.1,
  x = drim.padj$geneID,
  warn_missing = F
)


drim.prop = reshape2::melt(counts[counts$feature_id %in% proportions(d)$feature_id, ], id = c("gene_id", "feature_id"))
drim.prop = drim.prop[order(drim.prop$gene_id, drim.prop$variable, drim.prop$feature_id), ]

system.time({
  drim.prop = drim.prop %>%
    group_by(gene_id, variable) %>%
    mutate(total = sum(value)) %>%
    group_by(variable, add = TRUE) %>%
    mutate(prop = value / total)
})

# Convert the data.frame to wide-format data with reshape2::dcast
drim.prop = reshape2::dcast(drim.prop[, c(1, 2, 3, 6)], gene_id + feature_id ~ variable)

drim.mean = as.data.frame(sapply(c('ctrl', 'pert'),
                                 function(lvl)
                                   rowMeans(proportions(d)[, 3:ncol(proportions(d))][, meta$group == lvl, drop = FALSE])))
# drim.mean$gene_id <- proportions(d)$
# log2 fold change in proportions
drim.log2fcp = log2(drim.mean[2]/drim.mean[1])
colnames(drim.log2fcp) = "log2fcp"
rownames(drim.log2fcp) = proportions(d)$feature_id

# Merge to create result data
drimData = cbind(drim.prop[,1:2], drim.prop[, 3:ncol(drim.prop)])

# drimDTU = merge(drim.padj[, c("feature_id.1", "gene_id.1", "gene", "transcript")], drim.log2fcp, by.x = "feature_id.1", by.y = "row.names")
drimDTU <- drim.padj[, c("feature_id.1", "gene_id.1", "gene")]
drimDTU$log2fcp <- drim.log2fcp[drimDTU$feature_id.1,]

colnames(drimDTU)[1:3] <- c("feature_id", "gene_id", "pvalue")

# saveRDS(list(exp = drimData,
#              p_val = drimDTU[, c(1:3, 5)]), file = 'final_res.Rds')

final_res <-
  merge(drimData, drimDTU, by = c("feature_id", "gene_id"))

data.table::fwrite(
  final_res,
  file = output_file,
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

