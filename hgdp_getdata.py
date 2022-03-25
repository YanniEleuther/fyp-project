#!/usr/bin/python

from os import listdir
import json

alsworkdir = 'c9orf72_exphunter/'
fragilexdir = 'fmr1_exphunter/'
nonraresdir = 'nonrarevariants_exphunter/'
referencefile = "hgdp_references.tsv"
catalogfile = "nonrare_catalog.json"
headerlist = [
    "sample_id",
    "locus",
    "population",
    "call_1",
    "call_2",
    "\n",
]
with open(catalogfile, 'r') as fin:
    variants = json.loads(fin.read())
    locuslist = [i["LocusId"] for i in variants]

print(locuslist)

superpop_nameref = {
        "Central South Asia (HGDP)": "central_south_asia",
        "Europe (HGDP)": "europe",
        "America (HGDP)": "america",
        "Africa (HGDP)": "africa",
        "Middle East (HGDP)": "middle_east",
        "Oceania (SGDP),Oceania (HGDP)": "oceania",
        "Oceania (HGDP)": "oceania",
        "East Asia (HGDP)": "east_asia",
}

samplename_superpopref = dict()

with open(referencefile, 'r') as fin:
    for line in fin:
        templist = line.strip().split('\t')
        
        samplename_superpopref[templist[0]] = superpop_nameref[templist[3]]

# Get nonrares data
nonrares_dict = dict()
for filename in listdir(nonraresdir):
    samplename = filename.split('.')[0]
    with open(nonraresdir + filename, 'r') as fin:
        tempdict = json.loads(fin.read())
        
        nonrares_dict[samplename] = tempdict.pop("LocusResults")

# Get fmr1 data
fmr1_dict = dict()
for filename in listdir(fragilexdir):
    samplename = filename.split('.')[0]
    with open(fragilexdir + filename, 'r') as fin:
        tempdict = json.loads(fin.read())
        
        fmr1_dict[samplename] = tempdict.pop("LocusResults")
    
with open("fmr1output.tsv", "w+") as fout:
    fout.write('\t'.join(headerlist))
    
    for sample_id in fmr1_dict:
        try:
            inferred_calls = fmr1_dict[sample_id]["FMR1"]["Variants"]["FMR1"]["Genotype"].split("/")
        except KeyError:
            pass
        
        formatlist = [
            sample_id,
            samplename_superpopref[sample_id],
            inferred_calls[0],
            inferred_calls[1],
            '\n'
        ]
        fout.write('\t'.join(formatlist))
        

# Get c9orf72 data
c9orf72_dict = dict()
for filename in listdir(alsworkdir):
    samplename = filename.split('.')[0]
    with open(alsworkdir + filename, 'r') as fin:
        tempdict = json.loads(fin.read())

        c9orf72_dict[samplename] = tempdict.pop("LocusResults")
    
with open("c9orf72output.tsv", "w+") as fout:
   fout.write('\t'.join(headerlist))

   for sample_id in c9orf72_dict:
       try:
           inferred_calls = c9orf72_dict[sample_id]["C9ORF72"]["Variants"]["C9ORF72"]["Genotype"].split("/")
       except KeyError:
           pass

       formatlist = [
               sample_id,
               samplename_superpopref[sample_id],
               inferred_calls[0],
               inferred_calls[1],
               '\n']
       fout.write('\t'.join(formatlist))

with open("nonrares_output.tsv", "w+") as fout:
   fout.write('\t'.join(headerlist))

   for sample_id in nonrares_dict.keys():
       for locus in locuslist:
           try:
               inferred_calls = nonrares_dict[sample_id][locus]["Variants"][locus]["Genotype"].split("/")
           except KeyError:
               pass
            
           formatlist = [
               sample_id,
               locus,
               samplename_superpopref[sample_id],
               inferred_calls[0],
               inferred_calls[1],
               '\n'
           ]
           fout.write('\t'.join(formatlist))
