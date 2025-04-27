##This python script is intended to gather the trinucleotide context of all the
##mutations within a vcf file

##Usage:
##python extractMutSigs.py \
##        --vcf /path/to/vcf \
##        --ref /path/to/ref/fasta.fna \
##        --out /path/to/output.csv

#Import libraries---------------------------------------------------------------
import argparse
import gzip
import pandas as pd
from pyfaidx import Fasta
from collections import defaultdict
from tqdm import tqdm

#Reverse complement function----------------------------------------------------
def reverseComplement(seq):
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]

#DNA strand normalization function----------------------------------------------
def normalize(ref, alt, context):
    if ref in ['A', 'G']:
        refRc = reverseComplement(ref)
        altRc = reverseComplement(alt)
        contextRc = reverseComplement(context)
        return refRc, altRc, contextRc
    return ref, alt, context

#Parse VCF function-------------------------------------------------------------
def parseVcf(vcfPath, refGenome):
    trinucCounts = defaultdict(int)

    fasta = Fasta(refGenome, rebuild=False)

    opener = gzip.open if vcfPath.endswith('.gz') else open
    with opener(vcfPath, 'rt') as vcf:
        for line in tqdm(vcf):
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom, pos, _, ref, alt = fields[:5]

            #Focus only on SNPs and no INDELs
            if len(ref) == 1 and len(alt) == 1 and ref != alt:
                try:
                    pos = int(pos)
                    #Fetch the trinucleotide context (+- 1 nucleotide from SNP)
                    seq = fasta[chrom][pos - 2:pos + 1].seq.upper()
                    if len(seq) != 3 or 'N' in seq:
                        continue
                    
                    #Normalize strand orientation
                    normRef, normAlt, context = normalize(ref, alt, seq)

                    #Construct the mutation label as the full trinucleotide context before and after the change
                    beforeRef = context[0] + normRef + context[2]
                    afterRef = context[0] + normAlt + context[2]
                    mutLabel = f"{beforeRef}>{afterRef}"

                    #Count number of each unique trinuc mutation
                    trinucCounts[mutLabel] += 1
                
                except Exception as e:
                    continue

    return pd.DataFrame(list(trinucCounts.items()), columns=["Mutation", "Count"])

#Main---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, help="Path to VCF file")
    parser.add_argument("--ref", required=True, help="Path to reference genome (.fna/.fasta)")
    parser.add_argument("--out", required=True, help="Output CSV file path")
    args = parser.parse_args()

    print("Parsing VCF and extracting trinucleotide contexts...")
    df = parseVcf(args.vcf, args.ref)
    df.to_csv(args.out, index=False)
    print(f"Results saved to {args.out}")

#Start at main------------------------------------------------------------------
if __name__ == "__main__":
    main()
