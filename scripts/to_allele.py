#!/usr/bin/env python 

import sys 

from cphasing.alleles import PartigAllele, PartigRecords 

def main(args):
    fasta,  res, output = args
    pa = PartigAllele(fasta, output=output)
    pa.pr = PartigRecords(res)
    pa.pr.convert(fasta)
    
    pa.to_alleletable()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: to_allele.py <sample.fasta> <sample.res> <out.allele>")
        sys.exit()
    main(sys.argv[1:])