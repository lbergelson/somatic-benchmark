import sys
import csv

from pymongo import MongoClient

def main(args):
    if len(args) != 1:
        print("Usage: update_database.py input_file")
        exit(1)

    input_file = args[0]

    client = MongoClient("cgamongo.broadinstitute.org", 27017 )
    db = client.KnowledgeDB_production
    collection = db.performance
    
    with open(input_file, 'r') as input:
        dr = csv.DictReader(input, delimiter="\t")
        
        collection.insert(dr)    
        


if __name__ == "__main__":
        main(sys.argv[1:])



