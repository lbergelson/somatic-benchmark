#!/bin/sh -e

if [ $# != 5 ]
then
    echo "Usage runMutect.sh <normal bam> <tumor bam> <reference> <output dir>"
    exit 1 
fi

#must match <normal><tumor><reference><outputDir>
TOOLDIR=$1
NORMALBAM=$2
TUMORBAM=$4
REFERENCE=$4
OUTPUTDIR=$5

#Create output directory
mkdir -p $OUTPUTDIR

#Fill in the commands to run your caller here

MutectJar="/xchip/cga_home/louisb/gatk-protected/dist/Queue.jar"
IntervalFile='/xchip/cga/reference/hg19/gaf_20111020+broad_wex_1.1_hg19.bed'

java -jar -Xmx4g $MutectJar -S $TOOLDIR/MuTectPipeline.scala \
-tb $TUMORBAM \
-nb $NORMALBAM \
-o $OUTPUTDIR/out \
-config hg19-wex \
-c .01 \
--sg 12 \
-run -bsub


#convert to maflite
echo "converting to maflite"
"/usr/bin/perl" "/xchip/tcga/gdac_prod/applications/process_mgmt/firehose_task_registry/cga/CallstatsToMaflite/broadinstitute.org/cancer.genome.analysis/00162/14/call_stats_to_maflite.pl" \
$OUTPUTDIR/out.call_stats.txt \
"19" \
"FSTAR" \
"$OUTPUTDIR/out.maflite.maf" \
"tumor_f,init_t_lod,t_lod_fstar,t_alt_count,t_ref_count,judgement" 

#apply maf validation
echo "apply maf validation"
"/broad/software/free/Linux/redhat_5_x86_64/pkgs/sun-java-jdk_1.6.0-21_x86_64/bin/java" "-Xmx1g" \
"-jar" "/xchip/tcga/gdac_prod/applications/process_mgmt/firehose_task_registry/cga/ApplyMAFValidation/broadinstitute.org/cancer.genome.analysis/00163/23/ApplyMAFValidation.jar" \
"M=$OUTPUTDIR/out.maflite.maf" \
"OUTPUT_MAF=$OUTPUTDIR/out.maflite.maf.annotated" \
"MATCH_MODE=Sample" \
"V=/cga/tcga-gsc/svnreference/validation"

#oncotate
"sh" "/xchip/tcga/gdac_prod/applications/process_mgmt/firehose_task_registry/cga/Oncotator_v1/broadinstitute.org/cancer.genome.analysis/10202/47/oncotator.sh" \
"MAFLITE" "TCGAMAF" \
"$OUTPUTDIR/out.maflite.maf.annotated" \
"$OUTPUTDIR/final_calls.maf" \
"hg19" \
"/xchip/cga/reference/annotation/db/oncotator_v1_ds" \
"/xchip/cga/reference/annotation/db/tcgaMAFManualOverrides2.4.config" \
"/xchip/tcga/gdac_prod/applications/process_mgmt/firehose_task_registry/cga/Oncotator_v1/broadinstitute.org/cancer.genome.analysis/10202/47/" \
"CANONICAL" \
"/xchip/tcga/Tools/python_fh/python_fh_env/"

