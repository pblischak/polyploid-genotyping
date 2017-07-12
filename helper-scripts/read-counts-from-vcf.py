#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import argparse
import sys

def _main():
    """

    """
    parser = argparse.ArgumentParser(description="Get read counts from a VCF file.",
                                     add_help=True)
    parser.add_argument('--vcf', action="store", type=str, help="Input VCF file.", default="-")
    parser.add_argument('--ad', action="store", type=int, help="Column with allele depth (AD) info.")
    parser.add_argument('--minDP', action="store", type=int, help="filter by minimum read depth.")
    parser.add_argument('--prefix', action="store", type=str, help="prefix for output files.", default="reads")

    args   = parser.parse_args()
    vcf    = args.vcf
    ad     = args.ad
    minDP  = args.minDP
    prefix = args.prefix

    totfile = open(prefix+"-tot.txt", 'wa')
    altfile = open(prefix+"-alt.txt", 'wa')

    if vcf == "-":
        lines = [l.strip() for l in sys.stdin.readlines() if not l.startswith("#")]
    else:
        with open(vcf) as f:
            lines = [l.strip() for l in f.readlines() if not l.startswith("#")]
    for line in lines:
        fields = line.split()[9:]
        for fl in fields:
            if fl.startswith("."):
                print("-9\t", file=totfile, sep='', end='')
                print("-9\t", file=altfile, sep='', end='')
            else:
                depths = fl.split(":")[ad-1].split(",")
                tot = int(depths[0]) + int(depths[1])
                if tot < minDP or len(depths) > 2:
                    print("-9\t", file=totfile, sep='', end='')
                    print("-9\t", file=altfile, sep='', end='')
                else:
                    alt = int(depths[1])
                    print(tot, "\t", file=totfile, sep='', end='')
                    print(alt, "\t", file=altfile, sep='', end='')
        print("\n", file=totfile, sep='', end='')
        print("\n", file=altfile, sep='', end='')

if __name__ == "__main__":
    """
    Run the main script.
    """
    _main()
