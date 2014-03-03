#!/bin/sh -e

if [ $# != 4 ]
then
    echo "Usage runMutect.sh <normal bam> <tumor bam> <reference> <output dir>"
    echo "Requires a working installation of oncotator."
    echo "Please edit this file to set the gatk path."
    exit 1 
fi

#must match <normal><tumor><reference><outputDir>
NORMALBAM=$1
TUMORBAM=$2
REFERENCE=$3
OUTPUTDIR=$4

#Create output directory
mkdir -p $OUTPUTDIR

#Fill in the commands to run your caller here

MutectJar="/xchip/cga_home/louisb/gatk-protected/dist/GenomeAnalysisTK.jar"
IntervalFile='/xchip/cga/reference/hg19/gaf_20111020+broad_wex_1.1_hg19.bed'

java -jar $MutectJar \
'--analysis_type' 'MuTect' \
'--normal_sample_name' "NORMAL" \
'-I:normal' $NORMALBAM \
'--tumor_sample_name' "TUMOR" \
'-I:tumor' $TUMORBAM \
'--reference_sequence' $REFERENCE \
'--dbsnp' '/xchip/cga/reference/hg19/dbsnp_134_b37.leftAligned.vcf' \
'--cosmic' '/xchip/cga/reference/hg19/hg19_cosmic_v54_120711.vcf' \
'--normal_panel' '/xchip/cga/reference/hg19/refseq_exome_10bp_hg19_300_1kg_normal_panel.vcf' \
'--out' "${OUTPUTDIR}/call_stats.txt" \
'--vcf' "${OUTPUTDIR}/final.snps.vcf" \
'--only_passing_calls' \
'--coverage_file' "${OUTPUTDIR}/coverage.wig.txt" \
'--power_file' "${OUTPUTDIR}/power.wig.txt" \
'--downsample_to_coverage' '1000' \
'--enable_extended_output' \
'--intervals' $IntervalFile \
'--intervals' 20:2000000-4000000 \
'--interval_set_rule' INTERSECTION

#convert to maflite
echo "converting to maflite"
"/usr/bin/perl" "/xchip/tcga/gdac_prod/applications/process_mgmt/firehose_task_registry/cga/CallstatsToMaflite/broadinstitute.org/cancer.genome.analysis/00162/14/call_stats_to_maflite.pl" \
$OUTPUTDIR/call_stats.txt \
"19" \
"FSTAR" \
"$OUTPUTDIR/maflite.maf" \
"tumor_f,init_t_lod,t_lod_fstar,t_alt_count,t_ref_count,judgement" 

#apply maf validation
echo "apply maf validation"
"/broad/software/free/Linux/redhat_5_x86_64/pkgs/sun-java-jdk_1.6.0-21_x86_64/bin/java" "-Xmx1g" \
"-jar" "/xchip/tcga/gdac_prod/applications/process_mgmt/firehose_task_registry/cga/ApplyMAFValidation/broadinstitute.org/cancer.genome.analysis/00163/23/ApplyMAFValidation.jar" \
"M=$OUTPUTDIR/maflite.maf" \
"OUTPUT_MAF=$OUTPUTDIR/maflite.maf.annotated" \
"MATCH_MODE=Sample" \
"V=/cga/tcga-gsc/svnreference/validation"

#oncotate
"sh" "/xchip/tcga/gdac_prod/applications/process_mgmt/firehose_task_registry/cga/Oncotator_v1/broadinstitute.org/cancer.genome.analysis/10202/47/oncotator.sh" \
"MAFLITE" "TCGAMAF" \
"$OUTPUTDIR/maflite.maf.annotated" \
"$OUTPUTDIR/final_calls.maf" \
"hg19" \
"/xchip/cga/reference/annotation/db/oncotator_v1_ds" \
"/xchip/cga/reference/annotation/db/tcgaMAFManualOverrides2.4.config" \
"/xchip/tcga/gdac_prod/applications/process_mgmt/firehose_task_registry/cga/Oncotator_v1/broadinstitute.org/cancer.genome.analysis/10202/47/" \
"CANONICAL" \
"/xchip/tcga/Tools/python_fh/python_fh_env/"

