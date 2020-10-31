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

# 1) Dé-duplication en séquence “complète”
def read_fasta(amplicon_file, minseqlen):
	"""
	"""

    with open(amplicon_file, "r") as f_fasta:

        for ligne in f_fasta:

            if not ligne.startswith(">"):

                if len(ligne) >= minseqlen:

                    yield ligne



def dereplication_fulllength(amplicon_file, minseqlen, mincount):
	"""
	"""

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

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
	"""
	"""

	list_kmer = cut_kmer(sequence, kmer_size)

	for kmer in list_kmer:

		if kmer in kmer_dict:
			kmer_dict[kmer].append(id_seq)
		else :
			kmer_dict[kmer] = [id_seq]
	return kmer_dict

    

def search_mates(kmer_dict, sequence, kmer_size):
	"""
	"""

	return[i[0] for i in Counter([ids for kmer in cut_kmer(sequence, kmer_size) if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)]


def get_identity(alignment_list):
	"""
	"""
	nb_nuc_id = 0
	for i in range(len(alignment_list[0])):
		if alignment_list[0][i] == alignment_list[1][i]:
			nb_nuc_id += 1
	return nb_nuc_id/len(alignment_list[0])

def detect_chimera(perc_identity_matrix):
	"""
	"""
	som_et = 0
	seq1 = False
	seq2 = False

	for i in range(len(perc_identity_matrix)):
		som_et += statistics.stdev(perc_identity_matrix[i])
		if perc_identity_matrix[i][0] > perc_identity_matrix[i][1]:
			seq1 = True
		if perc_identity_matrix[i][0] < perc_identity_matrix[i][1]:
			seq2 = True
	if som_et/len(perc_identity_matrix) > 5 and seq1 and seq2:
		return True
	else :
		return False

def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
	"""
	"""
	dfr_lst = dereplication_fulllength(amplicon_file, minseqlen, mincount)

	for l in dfr_lst:
		chunks = get_chunks(l[0], chunk_size)

        chunk_mates = []
        for seq in chunks:
            mates = search_mates(kmer_dict, seq, kmer_size)
            chunk_mates.append(mates)
        com = []

        for j in range(len(chunk_mates)):
            com = common(com, chunk_mates[j])

        if len(com) > 1:
            for f in com[0:2]:
                sequ = get_chunks(no_chimere[f], chunk_size)
                perc_identity_matrix = [[]]
                for k, chunk in enumerate(chunks):
                    align = nw.global_align(chunk, sequ[k])
                    identite =  get_identity(align)
                    perc_identity_matrix[k].append(identite)
            chimera = detect_chimera(perc_identity_matrix)

		if not detect_chimera(perc_identity_matrix):
			yield l


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
	"""
	"""
	pass


def write_OTU(OTU_list, output_file):
	"""
	"""
	pass

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