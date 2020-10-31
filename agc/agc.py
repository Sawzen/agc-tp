#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

# 1) Dé-duplication en séquence “complète”
def read_fasta(amplicon_file, minseqlen):

    with open(amplicon_file, "r") as f_fasta:

        for ligne in f_fasta:

            if not ligne.startswith(">"):

                if len(ligne) >= minseqlen:

                    yield ligne



def dereplication_fulllength(amplicon_file, minseqlen, mincount):

    list_seq = read_fasta(amplicon_file, minseqlen)
    dico_seq = {}

    for seq in list_seq:

        if seq not in dico_seq:
            dico_seq[seq] = 1

        else :
            dico_seq[seq] += 1

    for sequence, count in sorted(dico_seq.items(), key=lambda item: item[1], reverse = True):

        if seq_count >= mincount:
            yield [sequence, count]

# 2) Recherche de séquences chimériques par approche “de novo”
def get_chunks(sequence, chunk_size):
	list_seg =[]


	for i in (range(0, len(sequence) , chunk_size)):
		if i+chunk_size<=len(sequence):
			list_seg.append(sequence[i:i+chunk_size])
	if len(list_seg)>=4:
		return list_seg

def cut_kmer(sequence, kmer_size):
	for k in (range(len(sequence)-kmer_size+1)):
		yield sequence[k:kmer_size]

def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()


if __name__ == '__main__':
    main()