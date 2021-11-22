#!/usr/bin/env python
'''
DESCRIPTION

    After downstream unmixing, we get a probablity matrix for each cell across bins probability a read is assigned to mark1. Use this prob matrix with bam file to split a bam file into mark1 and mark2

FOR HELP

    python split_double_BAM_JY.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2020-03-20
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import collections
import os

import collections
import itertools
import numpy as np
import random
import networkx as nx
import pysam
import pysamiterators
import matplotlib.colors
from importlib import reload
import pandas as pd

from singlecellmultiomics.bamProcessing import sorted_bam_file
from singlecellmultiomics.bamProcessing.bamToCountTable import coordinate_to_bins


def main():
    parser = argparse.ArgumentParser(description='After downstream unmixing, we get a probablity matrix for each cell across bins probability a read is assigned to mark1. Use this prob matrix with bam file to split a bam file into mark1 and mark2')
    parser.add_argument('-inbam', metavar='INFILE',
                        help='Input bam file')
    # /hpc/hub_oudenaarden/jyeung/data/dblchic/double_staining_output_downstream/unfixed_louvain2/SplitReads/MF_BM_unfixed_louvain2_clstr_by_louvain_K4m1_K27m3.removeNA_FALSE-prob_mat.K4m1-K27m3_to_K4m1.txt
    parser.add_argument('-inprobmat', metavar='INFILE',
                        help='Tab sep matrix file. Columns are cell names (first fcolumn is ""). Rows are genomic bins. Values are probability of reads in bin assigned to mark1.')
    parser.add_argument('-outdir', metavar='OUTDIR',
                        help='Output directory for bams. Full name to be specified in script')
    parser.add_argument('-mapq', metavar='INTEGER 0 to 60', default=40, type=int,
                        help='Minimum quality of read to be considered')
    parser.add_argument('-binsize', metavar='Genomic binsize', default=50000, type=int,
                        help='Binsize of genomic bins to consider (assumes row names are defined by nearest 50kb bins)')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
    parser.add_argument('--logfile', '-l', metavar='LOGFILE', default = None,
                        help='Write arguments to logfile')
    args = parser.parse_args()

    # store command line arguments for reproducibility
    CMD_INPUTS = ' '.join(['python'] + sys.argv)    # easy printing later
    # store argparse inputs for reproducibility / debugging purposes
    args_dic = vars(args)
    # ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.iteritems()]  # for python2
    ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.items()]  # for python3
    ARG_INPUTS = ' '.join(ARG_INPUTS)

    # Print arguments supplied by user
    if not args.quiet:
        if args.logfile is not None:
            sys.stdout = open(args.logfile, "w+")
        print(datetime.datetime.now().strftime('Code output on %c'))
        print('Command line inputs:')
        print(CMD_INPUTS)
        print ('Argparse variables:')
        print(ARG_INPUTS)


        # pathToTables = '/hpc/hub_oudenaarden/mflorescu/data/mnase/mm/mergedBAMs/'

        # prob = pd.read_csv(pathToTables + 'MF_BM_unfixed_clstr_by_topics.K4m1-K27m3.removeNA_FALSE-prob_mat.K4m1-K27m3_to_K4m1.txt',
        #                   sep='\t')

        prob = pd.read_csv(args.inprobmat, sep="\t")

        new = prob["Unnamed: 0"].str.split(':|-', n=3, expand=True)

        prob['chr'] = new[0]
        prob['start'] = new[1]
        prob['end'] = new[2]

        prob.set_index(['chr','start','end'], inplace = True)

        prob.drop(["Unnamed: 0"], axis = 1, inplace = True)


        # bamFile = "/hpc/hub_oudenaarden/mflorescu/data/mnase/mm/mergedBAMs/all_BM_K4m1_K27m3_200119.bam"
        bamFile = args.inbam
        wrote = 0

        infboth = os.path.join(args.outdir, "both.bam")
        infA = os.path.join(args.outdir, "splitted_A.bam")
        infB = os.path.join(args.outdir, "splitted_B.bam")

        with pysam.AlignmentFile(bamFile) as f:
            with sorted_bam_file(infboth, f) as both, sorted_bam_file(infA,origin_bam=f) as a, sorted_bam_file(infB,origin_bam=f) as b:
                for readId,(R1,R2) in enumerate(pysamiterators.MatePairIterator(f)):
                    if R1.mapping_quality < args.mapq & R2.mapping_quality < args.mapq:
                        continue  # one of two reads should have sufficient MAPQ. Less stringent. Should be OK?

                    if R1.is_duplicate:
                        continue

                    bin_start, bin_end = coordinate_to_bins(R1.get_tag('DS'), args.binsize, args.binsize)[0]
                    # Obtain prob:
                    bin_name = (f'chr{R1.reference_name}',str(bin_start),str(bin_end))
                    if not bin_name in prob.index:
                        continue
                    if R1.get_tag('SM') not in prob.columns:
                        continue
                    p = prob.loc[bin_name, R1.get_tag('SM')]
                    wrote+=1
                    group = 'A' if np.random.random()<=p else 'B'
                    R1.set_tag('Gr',group)
                    R2.set_tag('Gr',group)
                    if group=='A':
                        a.write(R1)
                        a.write(R2)
                    else:
                        b.write(R1)
                        b.write(R2)
                    both.write(R1)
                    both.write(R2)
        print("Number of reads written:" + str(wrote))


if __name__ == '__main__':
    main()
