args <- commandArgs(trailingOnly=TRUE)
print(args)


getIndels <- function(df){
        df[df$Variant_Type %in% c("INS","DEL"),]
}

getSnps <- function(df){
        df[!(df$Variant_Type %in% c("INS","DEL")),]
}

file <- args[1]
output <- args[2]

maf <- read.delim(file)

snps = getSnps(maf)
indels = getIndels(maf)


summary_table = NULL
summary_table$tp_snp = 0
summary_table$fp_snp = nrow(snps)
summary_table$novel_snp = 0
summary_table$fn_snp = 0

summary_table$tp_indel = 0
summary_table$fp_indel = nrow(indels)
summary_table$novel_indel = 0
summary_table$fn_indel = 0

write.table(summary_table, output, sep="\t", row.names= F, quote=F)

