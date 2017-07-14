#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import argparse
import sys
import numpy as np

def _main():
    """

    """
    parser = argparse.ArgumentParser(description="Get read counts from a VCF file.",
                                     add_help=True)
    parser.add_argument('--vcf', action="store", type=str, help="Input VCF file.", default="-")
    parser.add_argument('-n', '--nind', action="store", type=int, help="Number of individuals.")
    parser.add_argument('-s', '--nsites', action="store", type=int, help="Number of sites.")
    parser.add_argument('--ad', action="store", type=int, help="Column with allele depth (AD) info.", default=2)
    parser.add_argument('--prefix', action="store", type=str, help="prefix for output files.", default="reads")

    args   = parser.parse_args()
    nind   = args.nind
    nsites = args.nsites
    vcf    = args.vcf
    ad     = args.ad
    prefix = args.prefix

    tmat = np.zeros((nind,nsites), dtype=np.int)
    amat = np.zeros((nind,nsites), dtype=np.int)

    totfile = open(prefix+"-tot.txt", 'wa')
    altfile = open(prefix+"-alt.txt", 'wa')

    if vcf == "-":
        lines = [l.strip() for l in sys.stdin.readlines() if not l.startswith("#")]
    else:
        with open(vcf) as f:
            lines = [l.strip() for l in f.readlines() if not l.startswith("#")]
    assert len(lines) == nsites
    line_idx = 0
    for line in lines:
        fields = line.split()[9:]
        assert len(fields) == nind
        field_idx = 0
        for fl in fields:
            if fl.startswith("."):
                tmat[field_idx, line_idx] = -9
                amat[field_idx, line_idx] = -9
                #print("-9\t", file=totfile, sep='', end='')
                #print("-9\t", file=altfile, sep='', end='')
            else:
                depths = fl.split(":")[ad-1].split(",")
                tmat[field_idx, line_idx] = int(depths[0]) + int(depths[1])
                amat[field_idx, line_idx] = int(depths[1])
                #print(tot, "\t", file=totfile, sep='', end='')
                #print(alt, "\t", file=altfile, sep='', end='')
            field_idx += 1
        #print("\n", file=totfile, sep='', end='')
        #print("\n", file=altfile, sep='', end='')
        line_idx += 1

    for i in range(nind):
        for j in range(nsites):
            print(tmat[i,j], "\t", file=totfile, sep='', end='')
            print(amat[i,j], "\t", file=altfile, sep='', end='')
        print("\n", file=totfile, sep='', end='')
        print("\n", file=altfile, sep='', end='')

if __name__ == "__main__":
    """
    Run the main script.
    """
    _main()
