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

__author__ = "Romain Pholoppe"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Romain Pholoppe"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Romain Pholoppe"
__email__ = "pholoppero@eisti.eu"
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

def read_fasta(amplicon_file, minseqlen):
    if amplicon_file.endswith("gz"):
        with gzip.open(amplicon_file, "rt") as fichier:
            sequence = ""
            for ligne in fichier :
                if ligne.startswith(">"):
                    if len(sequence)>= minseqlen:
                        yield sequence
                    sequence = ""
                else:
                    sequence = sequence + ligne.strip()
            yield sequence
    else :
        with open(amplicon_file, "r") as fichier:
            sequence = ""
            for ligne in fichier :
                if ligne.startswith(">"):
                    if len(sequence)>= minseqlen:
                        yield sequence
                    sequence = ""
                else:
                    sequence = sequence + ligne.strip()
            yield sequence




def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    fichierFasta = read_fasta(amplicon_file, minseqlen)
    dictionnaire = {}
    for sequence in fichierFasta:
        if sequence in dictionnaire :
            dictionnaire[sequence]= dictionnaire[sequence]+1
        else :
            dictionnaire[sequence]= 1
    for seq, nb in sorted(dictionnaire.items(), key = lambda t: t[1], reverse = True):
        if nb >= mincount :
            yield [seq, nb]


def get_chunks(sequence, chunk_size):
    liste=[]
    for i in range(0, len(sequence), chunk_size):
        sous_seq = sequence[i:chunk_size+i]
        if len(sous_seq)==chunk_size:
            liste.append(sous_seq)
    return liste


def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    for i in range(len(sequence)-kmer_size+1):
        yield sequence[i:kmer_size+i]

def get_identity(alignment_list):
    compteur = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            compteur = compteur+1
    pourcentage = compteur / len(alignment_list[0]) *100
    return pourcentage

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    gen = dereplication_fulllength(amplicon_file, minseqlen, mincount)

    pass

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    with open(output_file, "w") as fichier_sortie:
        for i, OTU in enumerate(OTU_list):
            fichier_sortie.write(">OTU_{} occurrence:{}\n".format(i+1, OTU[1]))
            fichier_sortie.write("{}\n".format(fill(OTU[0])))

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    OTU = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)
    write_OTU(OTU, args.output_file)

if __name__ == '__main__':
    main()
