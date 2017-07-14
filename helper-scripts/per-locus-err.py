#!/usr/bin/env python

# per-locus-err.py
# Written by PD Blischak

from __future__ import print_function
import numpy as np
from sys import argv
import argparse

PHRED = {'!': 1.00, '"': 0.79433, '#': 0.63096, '$': 0.50119, '%': 0.39811,
         '&': 0.31623, '\'': 0.25119, '(': 0.19953, ')': 0.15849, '*': 0.12589,
         '+': 0.10000, ',': 0.07943, '-': 0.06310, '.': 0.05012, '/': 0.03981,
         '0': 0.03162, '1': 0.02512, '2': 0.01995, '3': 0.01585, '4': 0.01259,
         '5': 0.01000, '6': 0.00784, '7': 0.00631, '8': 0.00501, '9': 0.00398,
         ':': 0.00316, ';': 0.00251, '<': 0.00200, '=': 0.00158, '>': 0.00126,
         '?': 0.00100, '@': 0.00079, 'A': 0.00063, 'B': 0.00050, 'C': 0.00040,
         'D': 0.00032, 'E': 0.00025, 'F': 0.00020, 'G': 0.00016, 'H': 0.00013,
         'I': 0.00010, 'J': 0.00008, 'K': 0.00006}

def _main(infile, nind):
    """

    """
    with open(infile) as f:
        lines = f.readlines()
        for l in lines:
            per_ind_phreds = []
            for i in range(5, ((nind+1) * 3) + 2, 3):
                if l.split()[i] != "*":
                    per_ind_phreds.append(l.split()[i])
                else:
                    pass
            qual_scores = "".join(per_ind_phreds)
            if len(qual_scores) > 0:
                print(np.mean([PHRED[qual_scores[q]] for q in range(len(qual_scores))]))
            else:
                print("0.01")

if __name__ == "__main__":
    """
    Parse command line arguments using argparse and run the _main() function.
    """
    parser = argparse.ArgumentParser(description="Extract mean per locus error rates from a SAMtools pileup file.",
                                     add_help=True)
    parser.add_argument('-i', '--infile', action="store", required=True,
                        metavar='\b', type=str, help="input pileup file")
    parser.add_argument('-n', '--nind', action="store", required=True,
                        metavar='\b', type=int, help="number of individuals")

    args = parser.parse_args()
    infile = args.infile
    nind   = args.nind
    _main(infile, nind)
