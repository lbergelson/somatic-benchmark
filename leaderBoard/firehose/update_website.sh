libdir=$1
new_summary=$2

$tmp

python $libdir/update_database.py $new_summary
python $libdir/download_database.py  
Rscript $libdir/makeReport.R "tmp.tsv" leaderboard 
cp *.html /cga/tcga-gsc/benchmark/data 
