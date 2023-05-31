from Bio import SeqIO
import sys

ref_fasta, cur_dir = sys.argv[1:3]

with open(f'{cur_dir}/chrom.sizes', 'w') as outfile:
    for rec in SeqIO.parse(ref_fasta, 'fasta'):
        print(f"{rec.id}\t{len(rec)}", file=outfile)


