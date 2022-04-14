#!/usr/bin/python

import json
import os

# Get the truth data

truthfile = "benchmark_data/truthgenes.tsv"
truthgenes = "benchmark_data/genelist.txt"
referencefile = "benchmark_data/str_reference.json"
hipstrfile = "benchmark_data/hipstr_regions.vcf"
gangstrfile = "benchmark_data/gangstr_regions.vcf"
exphunter_dir = "benchmark_data/exphunterdict/"
tredparsefile = "benchmark_data/tredparsereport/repeatcalls.tsv"
outputfilename = "benchmark_data/benchmarkfile.json"

# Get name of first sample for indexing later on from .env (not in the repo)
# ALS patient data used for benchmarking cannot be shared.
FIRSTSAMPLENAME = os.getenv('FIRSTSAMPLENAME')

def getGenesforBenchmark(truth, catalog):
    with open(truth, 'r') as fin:
        startlist = fin.read().strip().split('\n')
    
    with open(catalog, 'r') as fin:
        compare_catalog = json.loads(fin.read())
    
    finallist = list()
    
    for i in compare_catalog:
        if i['locus'] in startlist:
            finallist.append(i['locus'])
        
    return finallist
    
def getSampleNames(truthfilename):
    with open(truthfilename, 'r') as fin:
        # Extract first column (names of samples)
        templist = [line.strip().split()[0] for line in fin]
    
    # Remove duplicates
    samplenames = [d for i, d in enumerate(templist) if d not in templist[:i]]
    
    return samplenames

def usecoordinatesforname(startcoord, catalog):
    for i in catalog:
        if i["start"] == startcoord:
            return i["locus"]
        
def getGangstrData():
    # Used later in the function to remove headers
    slicer = "GT:DP:Q:REPCN:REPCI:RC:ENCLREADS:FLNKREADS:ML:INS:STDERR:QEXP"

    # Get data from gangstr
    with open(gangstrfile, 'r') as fin:
        infoformats = list()
        for line in fin:
            if line.startswith("#CHROM"):
                vcfsamps = line.strip().split('\t')
            elif line.startswith("chr"):
                infoformats.append(line.strip())
            
    # Extract just the filenames
    firstsample = vcfsamps.index(FIRSTSAMPLENAME)
    vcfsamps = vcfsamps[firstsample::]
    
    # Get only the chromosome and the start coordinates of every STR
    print(len(infoformats))

    # Get dictionary for reference strs
    with open(referencefile, 'r') as fin:
        variants = json.loads(fin.read())
        
    gangstr_dict = dict()

    for i in infoformats:
        # Get locus name
        locus_name = usecoordinatesforname(i.split('\t')[1], variants)
        
        # Assign sample names to each locus
        gangstr_dict[locus_name] = dict()
        
        # Assign duplet to each sample name in the locus
        locusinfo = i.split('\t')
        locusinfo = locusinfo[locusinfo.index(slicer)+1::]
        
        for n, k in enumerate(locusinfo):
            sample = vcfsamps[n]
            targetlist = k.split(":")
            first, second = targetlist[3].split(",")
            
            gangstr_dict[locus_name][sample] = [first, second]
    
    return gangstr_dict

def getmotifnum(locusid, catalog):
    for i in catalog:
        if i["locus"] == locusid:
            return len(i["motif"])

def getHipstrData():
    # Used later in the function to remove headers
    slicer = "GT:GB:Q:PQ:DP:DSNP:DSTUTTER:DFLANKINDEL:PDP:PSNP:GLDIFF:AB:DAB:FS:ALLREADS:MALLREADS:FILTER"
    
    # Gather data from hipstr vcf
    with open(hipstrfile, 'r') as fin:
        infoformats = list()
        for line in fin:
            if line.startswith("#CHROM"):
                vcfsamps = line.strip().split("\t")
            elif line.startswith("chr"):
                infoformats.append(line.strip())
                
    
    # Extract just the filenames 
    firstsample = vcfsamps.index(FIRSTSAMPLENAME)
    vcfsamps = vcfsamps[firstsample::]
    
    # Get reference str json
    with open(referencefile, 'r') as fin:
        variants = json.loads(fin.read())
    
    hipstr_dict = dict()
    
    for i in infoformats:
        locusinfo = i.split(slicer, 1)[0]
        
        locusmetadata = locusinfo.split('.', 1)[0].split('\t')
        # Remove empty elements caused by split()
        del locusmetadata[-1]

        # Get locus name
        locus_name = locusmetadata[2]
        # Get motif length specific for 
        motif_num = getmotifnum(locus_name, variants)
        
        reference_allele = i.split('\t.\t.')[0].strip().split('\t')[3]
        templist = i.split('\t.\t.')[0].strip().split('\t')[4].split(',')
        alleles = [reference_allele] + templist
        
        for index, allele in enumerate(alleles):
            alleles[index] = int(len(allele) / motif_num)
            
        hipstr_dict[locus_name] = dict()
        
        # Subset for repeat count indices
        locusdata = i.split('\t')[i.split('\t').index(slicer)+1::]

        for n, k in enumerate(locusdata):
            sample = vcfsamps[n]
            targetlist = k.split(":")[0]
            
            # Use | as delimiter since there's no phasing
            subsetlist = targetlist.split("|")
            if subsetlist  == ['.']:
                first = None
                second = None
            else:
                first = alleles[int(subsetlist[0])]
                second = alleles[int(subsetlist[1])]
            
            hipstr_dict[locus_name][sample] = [first, second]

    return hipstr_dict        

def getExphunterData():
    exphunter_dict = dict()
    samplefiles = os.listdir(exphunter_dir)
    
    for filename in samplefiles:
        fulldir = "benchmark_data/exphunterdict/" + filename
        
        with open(fulldir, 'r') as fin:
            sampledict = json.loads(fin.read())
        
        # Get the name of the sample without the file extension
        samplename = os.path.splitext(filename)[0]
        exphunter_dict[samplename] = sampledict.pop("LocusResults")
    
    return exphunter_dict

def getTredparseData():
    tredparse_dict = dict()
    with open(tredparsefile, 'r') as fin:
        for line in fin:
            if line.startswith('SampleKey'):
                genelist = line.strip().split('\t')[2:]
            elif line.startswith('LP'):
                samplename = line.strip().split('\t')[0]
                callslist = line.strip().split('\t')[2:]
                tredparse_dict[samplename] = dict()
                for i,d in enumerate(genelist):
                    tredparse_dict[samplename][d] = callslist[i]
    
    with open('benchmark_data/tredparsereport/tredtest.json', 'w+') as fout:
        json.dump(tredparse_dict, fout, indent=4)
    
    return tredparse_dict
                

def mainbenchmark(hipstr, gangstr, exphunter, tredparse):            
    # Specific loci were manually removed because of a repeat motif that HipSTR considers too large
    # thus a custom loci list is used
    with open("benchmark_data/outputlist.txt", "r") as fin:
        genelist = [line.strip() for line in fin]
    
    gangstr_dict = gangstr
    hipstr_dict = hipstr
    exphunter_dict = exphunter
    tredparse_dict = tredparse
    
    with open('benchmark_data/exphuntertest.json', 'w+') as fout:
        json.dump(exphunter_dict, fout, indent=4) 

    # Initialise a nested dictionary into each sample name
    dataheaders = [
        "sample_id",
        "locus_id",
        "truth_1",
        "truth_2",
        "hipstr_1",
        "hipstr_2",
        "gangstr_1",
        "gangstr_2",
        "exphunter_1",
        "exphunter_2",
        "tredparse_1",
        "tredparse_2",
    ]

    with open(truthfile, 'r') as fin:
        with open('benchmark_data/outputdata.tsv', 'w+') as fout:
            
            # Write the headers
            fout.write('\t'.join(dataheaders) + "\n")

            for line in fin:

                linelist = line.strip().split('\t')
                
                missingkey = False
                if linelist[1] in genelist and int(linelist[2]) >= 0 and int(linelist[3]) >= 0:
                    targetsample = linelist[0]
                    targetlocus = linelist[1]
                    truth_counts = sorted(linelist[2:])
                    try:
                        hipstr_calls = hipstr_dict[targetlocus][targetsample]
                        gangstr_calls = gangstr_dict[targetlocus][targetsample]
                        exphunter_calls = exphunter_dict[targetsample][targetlocus]["Variants"][targetlocus]["Genotype"]
                        exphunter_calls = exphunter_calls.split('/')
                        tredparse_calls = tredparse[targetsample][targetlocus].replace('.', '0')
                        tredparse_calls = tredparse_calls.split('|')
                    except KeyError:
                        missingkey = True
                    
                    if missingkey:
                        continue
                    else:
                        formatlist = [
                            targetsample,
                            targetlocus,
                            str(truth_counts[0]),
                            str(truth_counts[1]),
                            str(hipstr_calls[0]),
                            str(hipstr_calls[1]),
                            str(gangstr_calls[0]),
                            str(gangstr_calls[1]),
                            str(exphunter_calls[0]),
                            str(exphunter_calls[1]),
                            str(tredparse_calls[0]),
                            str(tredparse_calls[1]),
                        ]
                        
                        if "None" not in formatlist:
                            fout.write("\t".join(formatlist) + "\n")

def orfbenchmark(gangstr, exphunter, tredparse):
    
    gangstr_dict = gangstr
    exphunter_dict = exphunter
    tredparse_dict = tredparse
    
    dataheaders = [
        "sample_id",
        "truth_1",
        "truth_2",
        "gangstr_1",
        "gangstr_2",
        "exphunter_1",
        "exphunter_2",
        "tredparse_1",
        "tredparse_2",
    ]
    with open(truthfile, 'r') as fin:
        with open('benchmark_data/c9orf72output_exphunter_gangstr.tsv', 'w+') as fout:
            fout.write("\t".join(dataheaders) + "\n")

            for line in fin:
                linelist = line.strip().split('\t')
                if linelist[1] == "C9ORF72":
                    targetsample = linelist[0]
                    truthcounts = linelist[2:]
                    
                    missingkey = False
                    try:
                        gangstr_genotype = gangstr_dict["C9ORF72"][targetsample]
                        exphunter_genotype = exphunter_dict[targetsample]["C9ORF72"]["Variants"]["C9ORF72"]["Genotype"]
                        exphunter_genotype = exphunter_genotype.split('/')
                        tredparse_genotype = tredparse_dict[targetsample]["C9ORF72"].replace('.', '0')
                        tredparse_genotype = tredparse_genotype.split('|')
                    except KeyError:
                        missingkey = True
                    
                    if missingkey:
                        continue
                    else:
                        formatlist = [
                            targetsample,
                            truthcounts[0],
                            truthcounts[1],
                            gangstr_genotype[0],
                            gangstr_genotype[1],
                            exphunter_genotype[0],
                            exphunter_genotype[1],
                            tredparse_genotype[0],
                            tredparse_genotype[1],
                        ]

                        if "None" not in formatlist:
                            fout.write("\t".join(formatlist) + "\n")

def main():
    gangstr_dict = getGangstrData()
    hipstr_dict = getHipstrData()
    exphunter_dict = getExphunterData()
    tredparse_dict = getTredparseData()

    mainbenchmark(hipstr_dict, gangstr_dict, exphunter_dict, tredparse_dict)
    orfbenchmark(gangstr_dict, exphunter_dict, tredparse_dict)
    
    
    
                
if __name__ == '__main__':
    main()
