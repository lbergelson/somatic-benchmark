new_summary=$1

$tmp

python update_database.py $new_summary
python download_database.py  
Rscript makeReport.R "tmp.tsv" leaderboard 
cp *.html /home/unix/louisb/private_html/
