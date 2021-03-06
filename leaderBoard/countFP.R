args <- commandArgs(trailingOnly=TRUE)
print("Args:")
print(args)


getIndels <- function(df){
        df[df$Variant_Type %in% c("INS","DEL"),]
}

getSnps <- function(df){
        df[!(df$Variant_Type %in% c("INS","DEL")),]
}

file <- args[1]
output <- args[2]

maf <- read.delim(file, comment="#")

snps = getSnps(maf)
indels = getIndels(maf)

coding <-   c("Missense_Mutation", "Nonstop_Mutation","Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del","In_Frame_Ins", "Nonsense_Mutation", "Splice_Site", "Silent")
snps <- snps[snps$Variant_Classification %in% coding,]
indels <- indels[indels$Variant_Classification %in% coding,]

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

