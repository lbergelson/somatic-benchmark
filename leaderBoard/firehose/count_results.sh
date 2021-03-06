if [ "$#" -ne 12 ]; then
    echo "usage: $0 <libdir> <run_type=NN|NormalNormal|KDB|HCC1143|HCC1954> <evaluation_maf> <comparison maf> <output prefix> <individual> <caller> <caller version> <task workspace name> <tumor bam> <normal bam> <interval file>"
	exit 1
fi

libdir=$1
run_type=$2
evaluation_maf=$3
comparison_maf=$4
output_prefix=$5
individual=$6
caller=$7
version=$8
workspace=$9
tumor=${10}
normal=${11}
interval=${12}

echo "libdir=$libdir"
echo "run_type=$run_type"
echo "evaluation_maf=$evaluation_maf"
echo "comparison_maf=$comparison_maf"
echo "output_prefix=$output_prefix"
echo "individual=$individual"
echo "caller=$caller"
echo "version=$version"
echo "workspace=$workspace"
echo "tumor=$tumor"
echo "normal=$normal"


. /broad/tools/scripts/useuse
reuse R-3.0

case "$run_type" in
"NormalNormal"|"NN") 
	echo "Counting false positives"
	Rscript $libdir/countFP.R $evaluation_maf ${output_prefix}.summary_kdb.txt
	;;
"KDB"|"HCC1143"|"HCC1954")
	echo "Comparing to kdb."
	Rscript $libdir/kdb_annotate.R $libdir $evaluation_maf $comparison_maf $output_prefix WEX $interval
	;;
*)
	echo "$run_type is not recognized."
	echo "Expected NN or KDB"
	exit 1
	;; 
esac 

###### Extra important information
timestamp=$( date +'%Y-%d-%d:%T' )
bonusheader="EvaluationGroup\tIndividual\tCaller\tVersion\tTime\tMaf\tWorkspace\tTumor\tNormal\n"
bonusInfo="$run_type\t$individual\t$caller\t$version\t$timestamp\t$evaluation_maf\t$workspace\t$tumor\t$normal"


output=${output_prefix}.summary_kdb.txt

paste $output <( echo -e $bonusheader $bonusInfo ) > ${output}.annotated

echo "Updating database with new results"
echo "uploading ${output}.annotated"
bash $libdir/update_website.sh $libdir $output.annotated
