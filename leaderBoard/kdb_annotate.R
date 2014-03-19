#===========================================
#
# KDB ANNOTATE
# DESCRIPTION: Will add annotation column to maf or call stats column based on true status in kdb
# INPUT: maf or call stats, kdb as text file
# OUTPUT: annotated maf or call stats and false negative file
#
#=============================================
library(methods)

args <- commandArgs(trailingOnly=TRUE)
print(args)

libdir = args[1]
mut.file = args[2] #mutation file - can be either call stats or maf format
kdb.file = args[3] #file from the kdb
pair = args[4] #pair id for naming of output folder
type = args[5]
interval.file = args[6]

source(paste0(libdir,"functions_mara.R"))

print("Reading files")
mut = read.table(as.character(mut.file), sep="\t", header=TRUE, quote='', stringsAsFactors=FALSE)
kdb = read.table(as.character(kdb.file), sep="\t", header=TRUE, quote='', stringsAsFactors=FALSE)
kdb$key = paste(kdb$Chromosome, kdb$Start_position, kdb$End_position, kdb$Variant_Type, kdb$Tumor_Seq_Allele2, sep=":")


if (file.exists(interval.file)) {
	print("subset to interval. interval file must be in bed format")
	colnames(interval) = c("Chromosome", "Start_position", "End_position")
	interval = read.table(interval.file, sep="\t", header=F, quote='', )
	kdb = interval_subset(af, interval)
} else {
	print("no interval file to subset to")
}

if (toupper(type) %in% c("WEX", "EXOME", "CODING")) {
	print("subseting to coding regions")
	cols = c("Missense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "Nonsense_Mutation", "Splice_Site", "Silent")
	print(cols)
	mut <- mut[mut$Variant_Classification %in% cols,]
	kdb = kdb[kdb$Variant_Classification %in% cols,]
}

kdb_tp = subset(kdb, kdb$Status %in% "TP")
kdb_fp = subset(kdb, kdb$Status %in% "FP")

annotating_maf <- function(maf, kdb_tp, kdb_fp) {
#	main implementation function to annotate maf by kdb
	if(nrow(maf) > 0){ 
		maf$kdb_annotation <- "Novel"
		maf[maf$key %in% kdb_tp$key,"kdb_annotation"] <- rep("TP", sum(maf$key %in% kdb_tp$key))
		maf[maf$key %in% kdb_fp$key, "kdb_annotation"] <- rep("FP", sum(maf$key %in% kdb_fp$key))
	}
	return(maf)
}

fn_sub <- function(maf, kdb_fp) {
# extract out false negatives
	false_neg = subset(kdb_fp, !(kdb_fp$key %in% maf$key))
	return(false_neg)
}

a = colnames(mut)
if ("Chromosome" %in% a) {
	print("This file is a maf.")
	if (sum((a %in% c("Chromosome", "Start_position", "End_position", "Variant_Type", "Tumor_Seq_Allele2")))!=5) {
		print("Mutation file requires Chromosome, Start_position, End_position, Variant_Type, and Tumor_Seq_Allele2")
	}
	mut$key = paste(mut$Chromosome, mut$Start_position, mut$End_position, mut$Variant_Type, mut$Tumor_Seq_Allele2, sep=":")
	out = annotating_maf(mut, kdb_tp, kdb_fp)
	fn = fn_sub(mut, kdb)
}

if ("contig" %in% a) {
	print("This file is a callstats")
	mut$Variant_Type <- "SNP"
	mut$key = paste(mut$contig, mut$position, mut$position, mut$Variant_Type, mut$alt_allele, sep=":")
	out = annotating_maf(mut, kdb_tp)
	fn = fn_sub(mut, kdb)
}



getIndels <- function(df){
	df[df$Variant_Type %in% c("INS","DEL"),]
}

getSnps <- function(df){
	df[df$Variant_Type == "SNP",]
}

indels <- getIndels(out)
snps <- getSnps(out)

fn_indels <- getIndels(fn)
fn_snps <- getSnps(fn)

summary_table = NULL
summary_table$tp_snp = sum(snps$kdb_annotation %in% "TP")
summary_table$fp_snp = sum(snps$kdb_annotation %in% "FP")
summary_table$novel_snp = sum(snps$kdb_annotation %in% "Novel") 
summary_table$fn_snp = dim(fn_snps)[1]

summary_table$tp_indel = sum(indels$kdb_annotation %in% "TP")
summary_table$fp_indel = sum(indels$kdb_annotation %in% "FP")
summary_table$novel_indel = sum(indels$kdb_annotation %in% "Novel") 
summary_table$fn_indel = dim(fn_indels)[1]

write.delim(fn, file=paste0(pair, ".false_negatives.maf"))
write.delim(out, file=paste0(pair, ".kdb_annotated.maf"))
write.delim(summary_table, file=paste0(pair, ".summary_kdb.txt"))

##TEST
#Rscript kdb_annotate.R /xchip/cga/gdac-prod/cga/jobResults/Oncotator_v1/HCC1954-wex_agilent_2/4908204/HCC1954-wex_agilent_2.snp.capture.maf.annotated /xchip/cga_home/mara/projets/cga_kdb/mongo/HCC1954.anotated HCC1954-wex-agilent
