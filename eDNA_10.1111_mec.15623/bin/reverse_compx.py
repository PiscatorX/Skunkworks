#!/usr/bin/env python
from Bio.Seq import Seq
from Bio import SeqIO
import sys
import argparse



def reverse_complimentx(seq_string, seq_file, informat, outfile):
    
    """
       
       get sequence and get reverse compliment

    """
    
    #custom script, light weight and easy customisatiom. Plus I had some free time
    #Numerous tools exist for this
    #for a better implementations see https://www.bioinformatics.nl/cgi-bin/emboss/revseq
    
    if seq_file:
        for rec in SeqIO.parse(seq_file, informat):
            rec.seq = rec.seq.reverse_complement()
            SeqIO.write(rec, outfile, informat)
            
    else:
        for seq in seq_string:
            my_seq = Seq(seq).reverse_complement()
            sys.stdout.write("{}\n".format(my_seq))


            
if __name__ == '__main__':
    parser = argparse.ArgumentParser("generates reverse complement")
    parser.add_argument('seq_string', nargs = "*", help ="sequence string DNA")
    parser.add_argument('-i','--informat', default="fasta")
    parser.add_argument('-s','--seq_file', type=argparse.FileType('r'))
    parser.add_argument('-o','--outfile', type=argparse.FileType('a'))
    
    args = parser.parse_args()
    err_msge = "provide one of two options. A sequence string e.g. 'ACGT' or sequence file e.g. '-s example.fasta'"
    if not args.seq_string:
        if not args.seq_file:
            parser.error(err_msge)
        elif not args.outfile:
            parser.error("provide an output filename e.g. '-o example_revcomp.fasta'")
    if (args.seq_string and args.seq_file):
        parser.error(err_msge)
        
    reverse_complimentx(args.seq_string, args.seq_file, args.informat, args.outfile)




