library(Nozzle.R1)
library(plyr)
dir.create( "reports", showWarnings=FALSE );

args <- commandArgs(trailingOnly=TRUE)
print(args)

stats_file = args[1] 
output_file = args[2]

stats <- read.delim(stats_file)

precision <- function(tp, fp){
	precision <- tp / (tp + fp)
	return(precision)
}

recall <- function(tp,fn){
	recall <- tp / (tp + fn)
	return(recall)
}

f1 <- function(precision, recall){
	f1 <- 2 * (precision * recall) / (precision + recall)
	f1
}

stats_summary <- ddply(stats, .(EvaluationGroup, Caller,Version,Time),
	 function(df) c(fp_snp = mean(df$fp_snp), tp_snp = mean(df$tp_snp), fn_snp = mean(df$fn_snp), novel_snp = mean(df$novel_snp),
		fp_indel = mean(df$fp_indel), tp_indel = mean(df$tp_indel), fn_indel = mean(df$fn_indel), novel_indel = mean(df$novel_indel) ))

stats_summary <- mutate(stats_summary, precision_snp = precision(tp_snp, fp_snp))
stats_summary <- mutate(stats_summary, precision_indel = precision(tp_indel, fp_indel))

stats_summary <- mutate(stats_summary, recall_snp = recall(tp_snp, fn_snp))
stats_summary <- mutate(stats_summary, recall_indel = recall(tp_indel, fn_indel))


stats_summary <- mutate(stats_summary, f1_snp = f1(precision_snp, recall_snp))
stats_summary <- mutate(stats_summary, f1_indel = f1(precision_indel, recall_indel))



# Phase 1: create report elements
r <- newCustomReport( "Caller LeaderBoard")
hcc1143 <- newSection( "HCC1143" );
hcc1954 <- newSection( "HCC1954" );
normalNormal <- newSection( "NormalNormal");
completeData <- newSection( "Complete");

common_cols <- c("Caller", "Version", "Time")
snp_cols = c(common_cols,c("tp_snp","fp_snp","fn_snp","novel_snp","precision_snp","recall_snp","f1_snp"))
indel_cols = c(common_cols, c("tp_indel","fp_indel","fn_indel","novel_indel","precision_indel","recall_indel","f1_indel"))
print(snp_cols)
print(indel_cols)

stat_1143_s <- stats_summary[stats_summary$EvaluationGroup == "HCC1143",snp_cols]
stat_1143_i <- stats_summary[stats_summary$EvaluationGroup == "HCC1143",indel_cols]

stat_1954_s <- stats_summary[stats_summary$EvaluationGroup == "HCC1954",snp_cols]
stat_1954_i <- stats_summary[stats_summary$EvaluationGroup == "HCC1954",indel_cols]

stat_nn <- stats_summary[stats_summary$EvaluationGroup =="NormalNormal", c(common_cols, c("fp_snp", "fp_indel"))]


t_1143_s <- newTable( stat_1143_s, "HCC1143 Snps" ); # w/ caption
t_1143_i <- newTable( stat_1143_i, "HCC1143 Indels" ); # w/ caption

t_1954_s <- newTable( stat_1954_s, "HCC1954 Snps" );
t_1954_i <- newTable( stat_1954_i, "HCC1954 Indels" );

t_nn <- newTable(stat_nn, "Normal-Normal Stats");

t_complete <- newTable(stats, "Complete unaggregated calls")

# Phase 2: assemble report structure bottom-up
hcc1143 <- addTo(hcc1143, t_1143_s,t_1143_i)
hcc1954 <- addTo(hcc1954, t_1954_s,t_1954_i)
normalNormal <- addTo(normalNormal, t_nn)
completeData <- addTo(completeData, t_complete)

r <- addTo(r, hcc1143, hcc1954, normalNormal, completeData)

# Phase 3: render report to file
writeReport( r, filename=output_file ); # w/o extension
