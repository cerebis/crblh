#!/usr/bin/env python
import argparse
import csv
import glob
import os
from collections import OrderedDict

import numpy as np
from scipy import sparse

parser = argparse.ArgumentParser(description='Create occupancy matrix per cluster')
parser.add_argument('--input-sep', default='\t', help='Value separator in output table [tab]')
parser.add_argument('--output-sep', default='\t', help='Value separator in output table [tab]')
parser.add_argument('input_dir', help='Input directory containg *.memb files')
parser.add_argument('output_table', help='Output CSV table')
args = parser.parse_args()

genome_names = []
og_names = set()
og_per_org = OrderedDict()
for n, memb_file in enumerate(glob.glob(os.path.join(args.input_dir, '*.memb'))):
    base, ext = os.path.splitext(os.path.basename(memb_file))
    genome_names.append(base)
    with open(memb_file, 'r') as ih:
        # print memb_file
        og_per_org[n] = OrderedDict()
        reader = csv.reader(ih, delimiter=args.input_sep)
        reader.next()  # skip header
        for row in reader:
            if row[1].startswith('iso'):
                continue
            og = int(row[1])
            og_names.add(og)
            if og in og_per_org[n]:
                og_per_org[n][og] += 1
            else:
                og_per_org[n][og] = 1

max_og = 0
for k in og_per_org.keys():
    k_max = max(map(int, og_per_org[k].keys()))
    if k_max > max_og:
        max_og = k_max

spmat = sparse.lil_matrix((len(og_per_org), max_og + 1))
for k1, row in og_per_org.iteritems():
    for k2, value in row.iteritems():
        spmat[(int(k1), int(k2))] = int(value)

with open(args.output_table, 'w') as out_h:
    og_names = sorted(og_names)
    writer = csv.writer(out_h, delimiter=args.output_sep)
    writer.writerow(['CLID'] + genome_names)
    dat = np.transpose(np.delete(spmat.todense(), 0, 1)).astype(np.int)
    nr, nc = np.shape(dat)
    for i in xrange(nr):
        writer.writerow([i + 1] + dat[i, ].tolist()[0])
