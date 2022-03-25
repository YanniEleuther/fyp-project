
import json

with open("variant_catalog.json", "r") as a_file:
    variantdata = json.loads(a_file.read())

normalrepeats = list()
rarerepeats = list()
for i in variantdata:
    if i["VariantType"] == "RareRepeat":
        rarerepeats.append(i)
    elif i["VariantType"] == "Repeat":
        normalrepeats.append(i)

with open("nonrare_catalog.json", "w+") as fout:
    json.dump(normalrepeats, fout)
with open("rare_catalog.json", "w+") as fout:
    json.dump(rarerepeats, fout)



