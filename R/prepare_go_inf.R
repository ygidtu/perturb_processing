setwd("OrySat/release52/")
library(tidyr)
library(stringr)
library(GO.db)
GO <- as.list(GOTERM)

do.call(rbind, lapply(GO, function(x) {
  return(c(x@GOID, x@Term, x@Ontology))
})) %>% data.frame -> GO
colnames(GO) <- c('GOID', 'Des', 'Cate')

Func_Anno <-
  read.table(
    "IRGSP-1.0_representative_annotation_2019-06-26.tsv.gz",
    sep = "\t",
    quote = "",
    header = TRUE,
    fill = TRUE
  )
GO_Gene <- Func_Anno[, c(2, 10)]
#去掉空值行(基因对应GO id为空)
GO_Gene <- GO_Gene[-which(GO_Gene$GO == ""), ]
#拆分第二列数据(Note:一个逗号一个GO)
GO_Sep <-
  separate(
    GO_Gene,
    col = GO,
    sep = ',',
    remove = TRUE,
    into = as.character(c(1:50))
  )
#检查GO是否完全分开
which(!is.na(GO_Sep[, 51]))         #判断最后一列是否有非NA值
sum(is.na(GO_Sep[, 51]))            #判断最后一列为NA值行数是否与矩阵行一样
#删除全为NA的列
GO_Sep <- GO_Sep[, -which(apply(GO_Sep, 2, function(x)
  all(is.na(x))))]
#按行合并
GO_Gene_Matr <- data.frame(matrix(NA, 300000, 2))
for (i in 2:ncol(GO_Sep)) {
  if (i == 2) {
    tmp <- as.matrix(GO_Sep[, c(1, i)])
    GO_Gene_Matr <- tmp
  } else{
    tmp <- as.matrix(GO_Sep[, c(1, i)])
    GO_Gene_Matr <- rbind(GO_Gene_Matr, tmp)
  }
}
#取出含GO的行
GO_Gene_Matr <- GO_Gene_Matr[grep("GO:", GO_Gene_Matr[, 2]), ]
GO <- data.frame(str_sub(GO_Gene_Matr[, 2], start = -11L, end = -2L))
GO_Gene_Matr_Final <- cbind(GO_Gene_Matr, GO)
GO_Gene_Matr_Final <- GO_Gene_Matr_Final[, c(1, 3)]
GO_Gene_Matr_Final <- unique(GO_Gene_Matr_Final)
colnames(GO_Gene_Matr_Final) <- c("GeneID", "GO")

saveRDS(GO_Gene_Matr_Final, file = 'Gene_go_inf.Rds')
saveRDS(GO, file = 'GO_inf.Rds')


