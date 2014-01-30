# ### Getting input file names  
if( interactive() ){
  #for developement and testing purposes
  print("I see you're running interactively, setting default values")
  outputdir <- "."
  inputfile <- "~/Downloads/An_Histo_Only.final_analysis_set.maf"
} else {
  args <- commandArgs(trailingOnly=TRUE)

  if(length(args)!=2){
    print("Usage: Rscript graphs_from_mafs.R <input.maf> <outputdirectory> ")
    quit()
  }
  inputfile <- args[1]
  print(paste("input maf =", inputfile))
  outputdir <- args[2]
  print(paste("output directory =", outputdir))
 }

if( ! file.exists(inputfile)){
    print("Input maf does not exist.  Exiting")
    quit()
} 
if( ! file.exists(outputdir) ){
    print("Output directory doesn't exist.  Creating it.")
    dir.create(outputdir)
}

###  Loading required libraries
library(ggplot2)
library(plyr)
library(gridExtra)
library(gtools)
library(gdata)


### Defining functions
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

    if (variant_classification %in% coding){
        return("Coding")
    }
    else {
        return("Non_Coding")
    }
}

### Preparing data

maf <- read.table(file=inputfile,header=TRUE, quote='', sep="\t", stringsAsFactors=FALSE)

maf$Pair_ID <- paste(gsub("-Tumor","",maf$Tumor_Sample_Barcode), "-",  gsub("-Normal","",maf$Matched_Norm_Sample_Barcode))
samples <- length( unique(maf$Pair_ID)) 

maf$Chromosome <- sort_chromosomes(maf)

maf$allele_fraction <- maf$t_alt_count / (maf$t_alt_count+maf$t_ref_count)

maf$Matches_COSMIC_Mutation <- NAToFalse( ! maf$COSMIC_overlapping_mutations == "")

maf$Overlaps_DB_SNP_Site <- NAToFalse(! maf$dbSNP_RS =="" )

maf$Classification <-  mapply(cosmic_or_dbsnp,maf$Matches_COSMIC_Mutation, maf$Overlaps_DB_SNP_Site)

maf$Coding <- mapply(coding_or_non_coding, maf$Variant_Classification)

maf$Indel_Length = mapply( calc_length, maf$Reference_Allele, maf$Tumor_Seq_Allele2, maf$Variant_Type)

maf <- mutate(maf, Tumor_Depth = t_alt_count+t_ref_count)

#### Making graphs
plot_percentage_and_count <- function(df, variable, name, outputdir){
  #Determine font size for name
  #name_length <- max(nchar(df$Pair_ID)) 
   

  #sort factor by length
  df$Pair_ID <- with(df, reorder(Pair_ID, Pair_ID, function(x) length(x)))

  percent <- ggplot(df, aes(x = Pair_ID)) + geom_bar(aes_string(fill = variable ), position = 'fill') +
    coord_flip() +theme_bw(base_family='Helvetica')+ scale_fill_brewer(palette="Paired") +
    theme(legend.position = "none", strip.text.x = element_text(size=6)) +labs(y="Percent")

  counts <- ggplot(df, aes(x = Pair_ID)) + geom_bar(aes_string(fill = variable)) + coord_flip() +
    theme_bw(base_family='Helvetica')+ scale_fill_brewer(palette="Paired") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank())

  g <- arrangeGrob(percent, counts, nrow=1)
  name_pieces <- c(outputdir, "/", name, ".pdf")
  filename <- paste(name_pieces, collapse='')
  print(paste("Saving ",filename, sep=''))
  ggsave(file=filename, g, height=max(samples/8,4), width=10, units="in", limitsize=FALSE)
}


draw_graphs <- function(basedir, subdir, maf){
    print(paste("drawing graphs for ", subdir))
    outputdir <- paste(basedir, subdir, sep="/")
    dir.create(outputdir)

    save_with_name <- function(name,height=10, width=10) {
      name_pieces <- c(outputdir, "/",subdir,"_",name, ".pdf")
      filename <- paste(name_pieces, collapse='')
      print(paste("Saving ",filename, sep=''))
      ggsave(file=filename, height=height, width=width, units="in", limitsize=FALSE)  
    }
    ggplot(maf, aes(x=Tumor_Depth, y= allele_fraction)) +stat_binhex(bins=100)+ scale_fill_gradientn( colours=c("white","red"))
    save_with_name("allele_fraction_vs_tumor_depth")
    
    ggplot(maf[maf$Variant_Type == "SNP",], aes(x=Tumor_Depth, y= allele_fraction)) +stat_binhex(bins=50)+ scale_x_log10() + scale_fill_gradientn( colours=c("white","red"))
    save_with_name("allele_fraction_vs_tumor_depth_snps")

    ggplot(maf[maf$Variant_Type %in% c("INS","DEL"),], aes(x=Tumor_Depth, y= allele_fraction)) +stat_binhex(bins=50)+ scale_x_log10() + scale_fill_gradientn( colours=c("white","red"))
    save_with_name("allele_fraction_vs_tumor_depth_indels")
    
    qplot(data=maf,x=allele_fraction, fill=Classification) + theme_bw()
    save_with_name("allele_fraction_all_samples") 
    
    qplot(data=maf, x=allele_fraction, fill=Classification) + facet_wrap(facets=~Pair_ID, ncol=4)+theme_bw() + theme(strip.text.x = element_text(size=6))
    save_with_name("allele_fraction_by_sample", height=max(samples/4,4))
    
    qplot(data=maf, x=allele_fraction, fill=Classification) + facet_wrap(facets=~Pair_ID, ncol=4, scales="free_y") +theme_bw() + theme(strip.text.x = element_text(size=6))
    save_with_name("allele_fraction_by_sample_normalized", height=max(samples/4,4), width=10)
    
    plot_percentage_and_count(maf, "Classification", "COSMIC_overlap_by_sample", outputdir)
    
    qplot(data=maf, x=Tumor_Depth) + theme_bw() + scale_x_log10()
    save_with_name("tumor_depth")

    qplot(data=maf, x=allele_fraction) +theme_bw()
    save_with_name("allele_fraction")

    #subset_to_indels
    indels_only <- maf[maf$Variant_Type %in% c("DEL","INS"), ]
    
    if (dim(indels_only)[1]!=0) {
    ggplot(data=maf)+ geom_bar(subset=.(Variant_Type == "DEL"), aes(x=Indel_Length,y=-..count..,fill=Variant_Classification,stat="identity")) + geom_bar(subset=.(Variant_Type=="INS"), aes(x=Indel_Length,y=..count.., fill=Variant_Classification, stat="identity")) + ylab("Deletions - Insertions") + theme_bw()
    save_with_name("stacked_indel_lengths", height=5, width=7)
    
    ggplot(data=indels_only)+ geom_bar( aes(x=Indel_Length,y=..count..,fill=Variant_Classification,stat="identity")) + facet_wrap(facets=~Variant_Type, drop=TRUE) + theme_bw()
    save_with_name("indel_lengths_by_type", height=5, width=7)
    
    plot_percentage_and_count(indels_only, "Variant_Classification", "indels_by_type", outputdir)
    } else {
    print("no indels in maf")
    }


    #subset to snps
    snps_only <- maf[!(maf$Variant_Type %in% c("DEL","INS")), ]
    qplot(data=snps_only, x=i_t_lod_fstar)+theme_bw()
    save_with_name("lod_score")

}

draw_graphs(outputdir, "all", maf)

draw_graphs(outputdir, "coding", maf[maf$Coding == "Coding",])
draw_graphs(outputdir, "non_coding",maf[maf$Coding == "Non_Coding",])
