#!/usr/bin/python

import json

with open("Data/Variant.json", "r") as a_file:
    variantdata = json.loads(a_file.read())

normalrepeats = list()
rarerepeats = list()
for i in variantdata:
    if i["VariantType"] == "RareRepeat":
        rarerepeats.append(i)
    elif i["VariantType"] == "Repeat":
        normalrepeats.append(i)

with open(".Temp/Non-Rare.json", "w+") as fout:
    json.dump(normalrepeats, fout)
with open(".Temp/Rare.json", "w+") as fout:
    json.dump(rarerepeats, fout)



