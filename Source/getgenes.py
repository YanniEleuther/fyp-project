#!/usr/bin/python

import json

with open("Data/Variant.json", "r") as a_file:
    variant_catalog = json.loads(a_file.read())

for i in variant_catalog:
    print(i["LocusId"])
