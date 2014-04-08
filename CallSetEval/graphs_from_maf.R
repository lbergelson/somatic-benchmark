###  Loading required libraries
print("Loading Libraries")
library(ggplot2)
library(plyr)
library(gridExtra)
library(gtools)
library(gdata)
library(scales)
library(intervals)
library(GenomicRanges)
library(reshape)
print("Done loading libraries")

options(error=traceback)


NAToFalse <- function(x){
    NAToUnknown(x, unknown=FALSE, force=TRUE)
}

sort_chromosomes <- function(df){
  return(factor(df$Chromosome, mixedsort(df$Chromosome))) 
}


cosmic_or_dbsnp <- function(is_cosmic, is_dbsnp){
  if (is_cosmic)
  {
    return("COSMIC")
  } else if(is_dbsnp){
    return("dbSNP")
  } else return("Unclassified")
}

calc_length <- function( ref, tumor, type){
  if(type == "INS"){
    return(nchar(tumor))
  } else if( type =="DEL"){
    return(nchar(ref))
  } else {
    return(0)
  }
}

coding_or_non_coding <- function(variant_classification){
    coding <-  c("Nonsense_Mutation","Nonstop_Mutation","Missense_Mutation", "Splice_Site", "Silent","Frame_Shift_Del",
         "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins")

    return( ifelse(variant_classification %in% coding, "Coding", "Non_Coding"))
}

df_to_granges <- function(df){
    gr <- GRanges( seqnames= Rle( df$Chromosome ),
                    ranges= IRanges( df$Start_position, end= df$End_position))
    return(gr)   
} 

##Read an interval_list file and return simplified GRanges object 
read_interval_file <- function(interval_file) {
    if( ! is.null(interval_file) & file.exists(interval_file)){
     print(paste("Loading interval file:", interval_file))
     intervals <- read.delim(interval_file)           
     return ( reduce(df_to_granges(intervals)) ) 

    } else {
      print(paste("Can't find interval file:", interval_file))
      quit()
    }
} 

in_interval <- function( chr, start, end, intervals) {
    variant <- GRanges( seqnames=Rle(c(chr)), ranges=(IRanges(c(start), end=c(end))))
    return( variant == intersect(variant, intervals))
}

### Preparing data
prepare_data <- function(inputfile, intervals) {
    print("Preparing Data")
    maf <- read.table(file=inputfile,header=TRUE, quote='', sep="\t", stringsAsFactors=FALSE)
    
    maf$Pair_ID <- paste0(gsub("-Tumor","",maf$Tumor_Sample_Barcode), "\n",  gsub("-Normal","",maf$Matched_Norm_Sample_Barcode))
    maf$Pair <- paste0(gsub("-Tumor","",maf$Tumor_Sample_Barcode),"-",gsub("-Normal","",maf$Matched_Norm_Sample_Barcode))
 
    maf$Chromosome <- sort_chromosomes(maf)
    maf$Chromosome <- droplevels(maf$Chromosome)

    maf <- maf[with(maf, order(Chromosome, Start_position)),]
    
   
    maf$allele_fraction <- maf$t_alt_count / (maf$t_alt_count+maf$t_ref_count)
    
    maf$Matches_COSMIC_Mutation <- NAToFalse( ! maf$COSMIC_overlapping_mutations == "")
    
    maf$Overlaps_DB_SNP_Site <- NAToFalse(! maf$dbSNP_RS =="" )
    
    maf$Classification <-  mapply(cosmic_or_dbsnp,maf$Matches_COSMIC_Mutation, maf$Overlaps_DB_SNP_Site)
    
    maf$Coding <- coding_or_non_coding( maf$Variant_Classification)
    
    maf$Indel_Length = mapply( calc_length, maf$Reference_Allele, maf$Tumor_Seq_Allele2, maf$Variant_Type)
    
    maf <- mutate(maf, Tumor_Depth = t_alt_count+t_ref_count)

    if(!is.null(intervals)){
        mafRanges <- df_to_granges(maf)
        maf$in_interval <- countOverlaps(mafRanges, intervals) > 0
    }     
    else{
        maf$in_interval <- FALSE
    }

    return(maf)
}

#### Making graphs
plot_percentage_and_count <- function(df, variable, name, outputdir){
  #Determine font size for name
  #name_length <- max(nchar(df$Pair_ID)) 
   
  samples <- length( unique(maf$Pair_ID)) 

  #sort factor by length
  df$Pair_ID <- with(df, reorder(Pair_ID, Pair_ID, function(x) length(x)))

  percent <- ggplot(df, aes(x = Pair_ID)) + geom_bar(aes_string(fill = variable ), position = 'fill') +
    coord_flip() +theme_bw(base_family='Helvetica')+ scale_fill_brewer(palette="Paired") +
    theme(legend.position = "none", strip.text.x = element_text(size=6), axis.text.y= element_text(size=4)) +labs(y="Percent")

  counts <- ggplot(df, aes(x = Pair_ID)) + geom_bar(aes_string(fill = variable)) + coord_flip() +
    theme_bw(base_family='Helvetica')+ scale_fill_brewer(palette="Paired") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank())

  g <- arrangeGrob(percent, counts, nrow=1, sub=textGrob(name, vjust=0.1)) 
  name_pieces <- c(outputdir, "/", name, ".pdf")
  filename <- paste(name_pieces, collapse='')
  print(paste("Saving",filename))
  ggsave(file=filename, g, height=max(samples/8,4), width=10, units="in", limitsize=FALSE)
}

add_title_info <- function(plot, title, footnote) {
    plot <- plot + ggtitle(gsub("_"," ", title))
    plot <- arrangeGrob(plot, sub=textGrob(footnote, x = 0, hjust = -0.1, vjust=0.1))
    return(plot)
}

bind_save_function <- function(outputdir, prefix){
    save_function <- function(name,height=10, width=10, plot=last_plot()) {
      name_pieces <- c(outputdir, "/",prefix,"_",name, ".pdf")
      filename <- paste(name_pieces, collapse='')
        
      plot <- add_title_info(plot, name, prefix)
      print(paste("Saving",filename))
      ggsave(file=filename, height=height, width=width, units="in", plot=plot, limitsize=FALSE)  
    }     
}

##############
#scale_x_tumor_depth(maf)
#takes a maf data frame as input, returns a log10 scaled x axis with ticks appropriate to the maf's tumor depth
##############
scale_x_tumor_depth <- function(maf) {
      breaks_to_use <- log_breaks(max(maf$Tumor_Depth))
      return(scale_x_log10(breaks=breaks_to_use, minor_breaks=c()))
}

log_breaks <- function(max){
    potential_breaks <- unique(c(seq(0,10),seq(10,100,10), seq(100,1000,100), seq(1000,10000,10000)))
    return( potential_breaks[potential_breaks <= max])
}



#graphs that make sense for snps and indels together or apart
shared_graphs <- function(maf, outputdir, prefix){
      save_with_name <- bind_save_function(outputdir, prefix) 
         
      samples <- length( unique(maf$Pair_ID))

      
      ggplot(maf, aes(x=Tumor_Depth, y= allele_fraction)) +stat_binhex(bins=100)+ scale_fill_gradientn( colours=c("white","red")) + scale_x_tumor_depth(maf)
      save_with_name("allele_fraction_vs_tumor_depth")

      ggplot(maf, aes(x=Tumor_Depth, y= t_alt_count)) +scale_y_log10(breaks = log_breaks(max(maf$t_alt_count)), minor_breaks = c()) +stat_binhex(bins=100)+ scale_fill_gradientn( colours=c("white","red")) + scale_x_tumor_depth(maf)
      save_with_name("altreads_vs_tumor_depth")

      qplot(data=maf,x=t_alt_count, fill=Classification, binwidth=1) + theme_bw()
      save_with_name("alt_reads_all_samples") 
 
      qplot(data=maf,x=allele_fraction, fill=Classification) + theme_bw()
      save_with_name("allele_fraction_all_samples") 
      
      qplot(data=maf, x=allele_fraction, fill=Classification) + facet_wrap(facets=~Pair_ID, ncol=4)+theme_bw() + theme(strip.text.x = element_text(size=4))
      save_with_name("allele_fraction_by_sample", height=max(samples/4,4))
      
      qplot(data=maf, x=allele_fraction, fill=Classification) + facet_wrap(facets=~Pair_ID, ncol=4, scales="free_y") +theme_bw() + theme(strip.text.x = element_text(size=4))
      save_with_name("allele_fraction_by_sample_normalized", height=max(samples/4,4), width=10)
      
      plot_percentage_and_count(maf, "Classification", paste0(prefix,"_COSMIC_overlap_by_sample"), outputdir)

      plot_percentage_and_count(maf, "in_interval", paste0(prefix,"_in_interval_by_sample"), outputdir)
      
      qplot(data=maf, x=Tumor_Depth, fill=Classification) + theme_bw() + scale_x_tumor_depth(maf)
      save_with_name("tumor_depth_all_samples")
      
      qplot(data=maf, x=Tumor_Depth, fill=Classification, facets= ~in_interval) + theme_bw() + scale_x_tumor_depth(maf)
      save_with_name("tumor_depth_all_samples_by_territory")


}



draw_graphs <- function(basedir, subdir, maf){
    
    print(paste("drawing graphs for", subdir))
    outputdir <- paste(basedir, subdir, sep="/")
    if( ! file.exists(outputdir) ){
        dir.create(outputdir)
    }
    
    save_with_name <- bind_save_function(outputdir, subdir)
    #first draw graphs for all variants
    shared_graphs(maf, outputdir, paste0(subdir,"_all_variants"))

    #then snps (and dnp, tnps for now...)
    snps_only <- maf[!(maf$Variant_Type %in% c("DEL","INS")), ]

    if( dim(snps_only)[1]!=0){
        shared_graphs(snps_only, outputdir,paste0(subdir,"_snps")) 
        
        if("i_t_lod_fstar" %in% colnames(snps_only)){
            qplot(data=snps_only, x=i_t_lod_fstar, fill=Classification)+scale_x_log10(breaks=log_breaks(max(snps_only$i_t_lod_fstar)), minor_breaks=c())+theme_bw()
            save_with_name("snp_lod_score")

            ggplot(snps_only, aes(x=Tumor_Depth, y= i_t_lod_fstar)) +stat_binhex(bins=100)+scale_y_log10(breaks=log_breaks(max(snps_only$i_t_lod_fstar)), minor_breaks=c())+ scale_fill_gradientn( colours=c("white","red")) + scale_x_tumor_depth(maf)
            save_with_name("snp_lod_score_vs_tumor_depth")    
        } else {
            print("no lod score in maf")
        }
    } else {
        print("no snps in maf")
    }

    #then indels only
    indels_only <- maf[maf$Variant_Type %in% c("DEL","INS"), ]
    
    if (dim(indels_only)[1]!=0) {
        shared_graphs(indels_only,outputdir, paste0(subdir,"_indels") )
        ggplot(data=maf)+ geom_bar(subset=.(Variant_Type == "DEL"), binwidth=1, aes(x=Indel_Length,y=-..count..,fill=Variant_Classification,stat="identity")) +
             geom_bar(subset=.(Variant_Type=="INS"), binwidth=1, aes(x=Indel_Length,y=..count.., fill=Variant_Classification, stat="identity")) +
             ylab("Deletions - Insertions") + theme_bw()
        save_with_name("stacked_indel_lengths", height=5, width=7)
        
        ggplot(data=indels_only)+ geom_bar( binwidth=1, aes(x=Indel_Length,y=..count..,fill=Variant_Classification,stat="identity"))+
            facet_wrap(facets=~Variant_Type, drop=TRUE) + theme_bw()
        save_with_name("indel_lengths_by_type", height=5, width=7)
        
        plot_percentage_and_count(indels_only, "Variant_Classification", "indels_by_type", outputdir)

    } else {
        print("no indels in maf")
    }


     
       
}

create_mutation_stats_report <- function(input_file_name, interval_file_name, maf, interval_size){
    require( Nozzle.R1 )
   
    maf[maf$Variant_Type %in% c("INS","DEL"),]$Variant_Type <- "INDEL"
    
    interval_mb <- interval_size / 1000000

    per_pair_counts <- ddply(maf, .( Coding, in_interval, Variant_Type, Pair), summarize, count = length(Variant_Type))    
    all_pair <- ddply(per_pair_counts, .(Coding, Variant_Type, Pair), summarize,  count=sum(count))

    sum_counts <- ddply(per_pair_counts, .(Coding, in_interval, Variant_Type), here(summarize), total = sum(count), mean = mean(count), rate = ifelse(all(in_interval), mean / interval_mb, NA))
    #all_counts <- ddply(sum_counts, .(Coding, Variant_Type), summarize, count = sum( count ), mean=sum(mean), rate = NA)
    ggplot(data=per_pair_counts, aes(x= in_interval, y=count,  fill=Variant_Type)) + geom_boxplot()+ facet_wrap(facets=~Coding) + geom_boxplot(data=all_pair, aes(x="All", y=count, fill=Variant_Type))
    ggsave("summary.pdf")  
    ggsave("summary.png", width=7, height=7) 
    # Phase 1: create report elements
    r <- newCustomReport( "Mutations Stats" )
    s <- newSection( "By Sample" )
    ss1 <- newSection( "Per Pair Counts" )
    
    boxplot <- newFigure("summary.png",fileHighRes="summary.pdf","Per pair counts of indels and snps.")
       
    
    ss2 <- newSection( "Totals" )
    pairs_table <- newTable(per_pair_counts , "Mutation Counts per Pair")
    totals_table <- newTable(sum_counts, "Overall Mutations") 
    p <- newParagraph( paste("Counts of coding and non-coding mutations from", input_file_name, ".\n Interval file:", interval_file_name, "contains", interval_mb, "Mb"))

   # Phase 2: assemble report structure bottom-up
    ss1 <- addTo( ss1, pairs_table ); # parent, child_1, ..., child_n 
    ss2 <- addTo( ss2, totals_table, boxplot );
    s <- addTo( s, ss1, ss2, p );
    r <- addTo( r, s );
    
    # Phase 3: render report to file
    writeReport( r); # w/o extension
}

#### Getting input file names  

### parse options
suppressPackageStartupMessages(require(optparse))
option_list <- list(                                                                                                                        
    make_option(c("-o", "--outputdir"), action="store", type="character", default=".",  help="Output directory"),
    make_option(c("-l", "--interval_file"), action="store", type="character", default=NULL, help="Territory file (must be in bed format)")
)

parsed <- parse_args(OptionParser(option_list=option_list, usage = "Rscript %prog [options] <maf_file>"), print_help_and_exit=TRUE, positional_arguments=1)
opt <- parsed$options
positional <- parsed$args

print(paste("positional:",positional))
print(paste("opts:",opt))

if(length(positional) != 1){
    stop("maf file is a required parameter")
}

inputfile <- positional[1]
outputdir <- opt$outputdir
interval_file <- opt$interval_file

print(paste("input maf =", inputfile))
print(paste("output directory =", outputdir))
print(paste("interval file =", interval_file))

if( ! file.exists(inputfile)){
    print("Input maf does not exist.  Exiting")
    stop()
} 
if( ! file.exists(outputdir) ){
    print("Output directory doesn't exist.  Creating it.")
    dir.create(outputdir)
}

 

intervals <- read_interval_file(interval_file)
maf <- prepare_data(inputfile, intervals)

create_mutation_stats_report(inputfile, interval_file, maf, sum(width(intervals)) )

draw_graphs(outputdir, "all", maf) 
draw_graphs(outputdir, "coding", maf[maf$Coding == "Coding",])
draw_graphs(outputdir, "non_coding",maf[maf$Coding == "Non_Coding",])

draw_graphs(outputdir, "in_interval", maf[maf$in_interval == TRUE,])
draw_graphs(outputdir, "not_in_interval", maf[maf$in_interval == FALSE,])
