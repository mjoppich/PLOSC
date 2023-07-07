import argparse
import gzip, sys
from collections import Counter
from Bio import SeqIO, bgzf

import itertools

def hammingDist(str1, str2, fillchar = '-'):
    return sum([ch1 != ch2 for (ch1,ch2) in itertools.zip_longest(str1, str2, fillvalue = fillchar)])


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument("-b1", "--barcode1", type=argparse.FileType('r'))
    parser.add_argument("-b2", "--barcode2", type=argparse.FileType('r'))
    
    parser.add_argument("-u", "--umi", type=argparse.FileType('r'))
    parser.add_argument("-r", "--reads", type=argparse.FileType('r'))

    parser.add_argument('-o', '--output', type=str)

    parser.add_argument('-gex', '--gex', type=str, nargs='+')
    parser.add_argument('-ab', '--ab', type=str, nargs='+')

    args = parser.parse_args()

    compressOutput = True

    outputABUmi = args.output + "_ab_S1_L001_R1_001.fastq" + ("" if not compressOutput else ".gz")
    outputABRead = args.output + "_ab_S1_L001_R2_001.fastq"+ ("" if not compressOutput else ".gz")

    outputGEXUmi = args.output + "_gex_S1_L001_R1_001.fastq"+ ("" if not compressOutput else ".gz")
    outputGEXRead = args.output + "_gex_S1_L001_R2_001.fastq"+ ("" if not compressOutput else ".gz")

    print(outputABUmi)
    print(outputABRead)
    print(outputGEXUmi)
    print(outputGEXRead)

    sys.stdout.flush()
    
    newGex = []
    for x in args.gex:
        assert x.count(";") == 1
        newGex.append( tuple( x.split(";") ) )
    args.gex = newGex
    
    print("GEX Indices", args.gex)
    
    newAb = []
    for x in args.ab:
        assert x.count(";") == 1
        newAb.append( tuple( x.split(";") ) )
    args.ab = newAb
    
    print("AB Indices", args.ab)

    #GEX = CACGCCTT	GTATATAG	TCTCGGGC	AGGATACA
    #AB = CTGCGGCT	GACTCAAA	AGAAACTC	TCTGTTGG

    unknownBCs = Counter()

    bcFQin1 = SeqIO.parse(gzip.open(args.barcode1.name, "rt"), "fastq")
    bcFQin2 = SeqIO.parse(gzip.open(args.barcode2.name, "rt"), "fastq")
    umiFQin = SeqIO.parse(gzip.open(args.umi.name, "rt"), "fastq")
    readFQin = SeqIO.parse(gzip.open(args.reads.name, "rt"), "fastq")

    #with bgzf.BgzfWriter(outputABUmi, "wb") as outABUmi, bgzf.BgzfWriter(outputABRead, "wb") as outABRead, bgzf.BgzfWriter(outputGEXUmi, "wb") as outGEXUmi, bgzf.BgzfWriter(outputGEXRead, "wb") as outGEXRead:
    if compressOutput:
        outABUmi = bgzf.BgzfWriter(outputABUmi, "wb")
        outABRead = bgzf.BgzfWriter(outputABRead, "wb")
        outGEXUmi = bgzf.BgzfWriter(outputGEXUmi, "wb")
        outGEXRead = bgzf.BgzfWriter(outputGEXRead, "wb")
    else:
        outABUmi = open(outputABUmi, "w")
        outABRead = open(outputABRead, "w")
        outGEXUmi = open(outputGEXUmi, "w")
        outGEXRead = open(outputGEXRead, "w")

    mismatchCounter = {"GEX": Counter(), "AB": Counter()}

    for elemID, (bcRecord1, bcRecord2, uRecords, rRecord) in enumerate(zip(bcFQin1, bcFQin2, umiFQin, readFQin)):
        
        bcSeq1 = str(bcRecord1.seq)
        bcSeq2 = str(bcRecord2.seq)

        if elemID % 100000 == 0:
            print(sum(mismatchCounter["GEX"].values()), "+", sum(mismatchCounter["AB"].values()), "=", sum(mismatchCounter["GEX"].values()) + sum(mismatchCounter["AB"].values()), "/", elemID)
            print("GEX", mismatchCounter["GEX"])
            print("AB", mismatchCounter["AB"])
            sys.stdout.flush()
            mismatchCounter = {"GEX": Counter(), "AB": Counter()}

        gexDistances1 = min([hammingDist(bcSeq1, x[0]) for x in args.gex])
        abDistances1 = min([hammingDist(bcSeq1, x[0]) for x in args.ab])

        gexDistances2 = min([hammingDist(bcSeq2, x[1]) for x in args.gex])
        abDistances2 = min([hammingDist(bcSeq2, x[1]) for x in args.ab])
        
        gexDistances = max([gexDistances1, gexDistances2])
        abDistances = max([abDistances1, abDistances2])

        if gexDistances > 1 and abDistances > 1:
            continue

        if gexDistances == abDistances:
            print("Same distances", gexDistances, abDistances, bcSeq1, bcSeq2)
            continue
        if gexDistances < abDistances:
            mismatchCounter["GEX"][gexDistances] += 1
            SeqIO.write(sequences=uRecords, handle=outGEXUmi, format="fastq")
            SeqIO.write(sequences=rRecord, handle=outGEXRead, format="fastq")
        elif abDistances < gexDistances:
            mismatchCounter["AB"][abDistances] += 1
            SeqIO.write(sequences=uRecords, handle=outABUmi, format="fastq")
            SeqIO.write(sequences=rRecord, handle=outABRead, format="fastq")
        else:
            pass
            #unknownBCs[bcSeq] += 1

    try:
        outGEXUmi.close()
    except:
        pass

    try:
        outGEXRead.close()
    except:
        pass

    try:
        outABUmi.close()
    except:
        pass

    try:
        outABRead.close()
    except:
        pass


    print("Unknown Barcodes")
    for x in unknownBCs:
        print(x, unknownBCs[x])

