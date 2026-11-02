#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd

args = sys.argv[1:]

kmers = args[0]

mapping_bed = pd.read_csv(kmers, sep='\t', header=None)

mapping_bed.columns = ['chrom', 'chromStart', 'chromEnd', 'kmer', 'p_value']

bin_labels = [1,2,3,4,5,6,7,8,9,10][::-1]

mapping_bed['bins'] = pd.cut(mapping_bed['p_value'], bins=10, labels=bin_labels)

colour_dict = {
        1 : '94,79,162',
        2 : '50,136,189',
        3 : '102,194,165',
        4 : '171,221,164',
        5 : '230,245,152',
        6 : '254,224,139',
        7 : '253,174,97',
        8 : '244,109,67',
        9 : '213,62,79',
        10 : '158,1,66'
        }

mapping_bed['colours'] = mapping_bed['bins'].apply(lambda x: colour_dict[x])

mapping_bed['strand'] = '+'

mapping_bed.to_csv(sys.stdout, sep='\t', columns=['chrom', 'chromStart', 'chromEnd', 'kmer', 'p_value', 'strand', 'chromStart', 'chromEnd','colours'], index=False, header=False)
