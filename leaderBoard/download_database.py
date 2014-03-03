import sys
import csv
import itertools
import subprocess

from pymongo import MongoClient

def main(args):

    client = MongoClient("cgamongo.broadinstitute.org", 27017 )
    db = client.KnowledgeDB_production
    collection = db.performance
    
    tmpFile = "tmp.tsv"

    with open(tmpFile, 'w') as tmp:
   
        items = collection.find()
        #We're going to iterate twice, so we just copy everything into memory
        storedItems = list(items)
        
        #Now we iterate through and find all the field names
        keylists = [i.keys() for i in storedItems]
        keys = itertools.chain.from_iterable(keylists)
        names = set(keys)
        
        #Create the writer
        dr = csv.DictWriter(tmp, delimiter="\t", fieldnames=names)

        #Write the values
        dr.writeheader()
        for row in storedItems:
            dr.writerow(row)
        
        




if __name__ == "__main__":
        main(sys.argv[1:])



