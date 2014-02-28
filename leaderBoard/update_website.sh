new_summary=$1

$tmp

python update_database.py $new_summary
python download_database.py  
Rscript makeReport "tmp.tsv" leaderboard 
