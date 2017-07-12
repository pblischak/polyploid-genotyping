#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import argparse
import math

def _main():
    """

    """
    parser = argparse.ArgumentParser(description="Get genotypes from a VCF file.",
                                     add_help=True)
    parser.add_argument('--vcf', action="store", type=str, help="Input VCF file.")
    parser.add_argument('--prefix', action="store", type=str, help="prefix for output files.", default="gt")

    args   = parser.parse_args()
    vcf    = args.vcf
    prefix = args.prefix

    gtfile = open(prefix+"-out.txt", 'wa')

    with open(vcf) as f:
        lines = [l.strip() for l in f.readlines() if not l.startswith("#")]
        for line in lines:
            fields = line.split()[9:]
            for fl in fields:
                if fl.startswith("."):
                    print("NA\t", file=gtfile, sep='', end='')
                else:
                    gt = [int(g) for g in fl.split(":")[0].split("/")]
                    print(sum(gt), "\t", file=gtfile, sep='', end='')
            print("\n", file=gtfile, sep='', end='')

if __name__ == "__main__":
    """
    Run the main script.
    """
    _main()
