#!/usr/bin/python


import json

filename = "genelist.txt"
reference = "str_reference.json"

with open(filename, 'r') as fin:
    genelist = fin.read().strip().split('\n')

# Get str_reference

with open(reference, 'r') as fin:
    variant_catalog = json.loads(fin.read())

finallist = list()
for i in variant_catalog:
    if i['locus'] in genelist:
        finallist.append(i['locus'])

with open("outputlist.txt", "w+") as fout:
    fout.write("\n".join(finallist))
