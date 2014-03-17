space=$1

echo "Checking if the workspace:$space exists already."

if fiss space_exists -n $space
then
	echo "$space already exists!"
	echo "If you REALLY want to overwrite this space run this first:\n"
	echo "fiss space_delete $space"
	exit 2
else
	echo "creating $space"
	fiss space_new $space
fi

echo "uploading pairs"
fiss sample_import $space samples.tsv
fiss pair_import $space pairs.tsv
fiss pset_import $space pair_sets.tsv


