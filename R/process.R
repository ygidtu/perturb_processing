
dat <-
  lapply(c(
    'SE_protein_isoforms_neg_IncDiff.deduped.fa',
    'SE_protein_isoforms_pos_IncDiff.deduped.fa'
  ), function(x) {
    data.table::fread(file = x, sep  = '\t',header=T)
  })


dat <- do.call(rbind, dat)
colnames(dat) <- 'V1'
dat <- dat[startsWith(dat$V1, '>c'),]


lapply(dat$V1[-1], function(x) {
  nmd_inf <- sapply(strsplit(x, split = '\t'), '[', 2)
  info <- strsplit(strsplit(gsub(pattern = '^>', replacement = '', x), split = '[,][+]')[[1]][1], split=',')[[1]]
  res <- c(info, nmd_inf)
  res
}) -> res


res <- do.call(rbind, res)
colnames(res) <- c('chr', 'upstreamEE', 'downsteamES', 'Trans_symbol', 'score', 'Triger_NMD')

