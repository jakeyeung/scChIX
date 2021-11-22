#!/usr/bin/env python
'''
DESCRIPTION

    Remove duplicated rows in a bam file

FOR HELP

    python remove_duplicated_rows_bams.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2021-07-31
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import os
import pysam
from collections import defaultdict
# import pysamiterators.iterators as pits
# copy instead below 

global qnames_seen1
global qnames_seen12

qnames_seen1 = set()
qnames_seen2 = set()

def read_pair_generator_remove_dupes(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name

        # handle duplicates, remove them
        if read.is_read1 and qname in qnames_seen1:
            continue
        if read.is_read2 and qname in qnames_seen2:
            continue

        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
                qnames_seen1.add(qname)
            else:
                read_dict[qname][1] = read
                qnames_seen2.add(qname)
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]
def main():
    parser = argparse.ArgumentParser(description='Remove duplicated rows in a bam file')
    parser.add_argument('infile', metavar='INFILE',
                        help='Input bam')
    parser.add_argument('outfileprefix', metavar='OUTFILEPREFIX',
                        help='Output bam file name without .bam suffix')
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

    outunsorted = args.outfileprefix + ".unsorted.bam"
    outsorted = args.outfileprefix + ".sorted.bam"

    assert not os.path.isfile(outunsorted)
    assert not os.path.isfile(outsorted)

    with pysam.AlignmentFile(args.infile, "rb") as inbam:
        with pysam.AlignmentFile(outunsorted, "wb", template = inbam) as outbam:
            for read1, read2 in read_pair_generator_remove_dupes( inbam ):
                outbam.write(read1)
                outbam.write(read2)

    pysam.sort(outunsorted, "-o", outsorted)
    pysam.index(outsorted)
    assert os.path.isfile(outsorted)
    os.remove(outunsorted)

if __name__ == '__main__':
    main()
