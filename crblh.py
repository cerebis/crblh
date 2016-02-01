#!/usr/bin/env python
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.Alphabet import generic_protein
import re
import networkx as nx
import community as com
import numpy as np
import subprocess
import sys
import os
import os.path
import csv
import argparse

SEQFILE_MAP = "input_seqs.txt"

class Gene:
    def __init__(self, organism, name):
        self.organism = organism
        self.name = name
        self.id = organism + '-' + name

    def __repr__(self):
        return self.id

    def __str__(self):
        return self.id

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        if isinstance(other, Gene):
            return other.id == self.id
        return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        return self.id < other.id


# Calculate the modularity score for a given graph.
# - valued between 0 and 1, with 1 being highly modular.
def calc_modularity(g):
    return com.modularity(com.best_partition(g), g)


# Decompose a graph into constituent sub-graphs until
# all graphs are below a threshold modularity.
# returns a list of sub-graphs.
def decompose_graph(g, max_mod, decomposed = None):
    if decomposed is None:
        decomposed = []

    p = com.best_partition(g)
    if com.modularity(p, g) < max_mod:
        decomposed.append(g)
    else:
        # split communities
        part_ids = np.unique(p.values())
        for pi in part_ids:
            gi = g.copy()
            for n in g.nodes_iter():
                if p[n] != pi:  # remove all nodes not in partition
                    gi.remove_node(n)
            decompose_graph(gi, max_mod, decomposed)
    return decomposed


# Implementation of title2ids function for use with Bio.SeqIO.FastaIO.FastaIterator
#
# Splits FASTA headers at the first whitespace. Word[0] becomes the id, with the remainder
# forming the description.
def fasta_title2ids_simple(header):
    tok = re.split('\s', header, maxsplit=1)
    if len(tok) == 1:
        return tok[0], '', ''
    else:
        return tok[0], '', tok[1]

parser = argparse.ArgumentParser(description='Clustered Reciprocal Best Last Hit (RBLH)')
parser.add_argument('work_dir', metavar='WORK_DIR', nargs=1, help='Working directory')
parser.add_argument('output', metavar='OUTPUT', nargs=1, help='Output basename')
parser.add_argument('prot_files', metavar='PROTEIN_FILE', nargs='+', help='Proteins per organism')
parser.add_argument('--write-subgraphs', dest='writesubs', action='store_true', default=False, help='Write subgraphs')
parser.add_argument('--decomp', dest='decompose', action='store_true', default=False,
                    help='Decompose groups with modularity above max-mod')
parser.add_argument('--max-mod', type=float, metavar='THRESHOLD', nargs=1, default=0.2,
                    help='Maximum acceptable modularity score before decomposition (default=0.2)')
parser.add_argument('--bit-score', dest='bitscore', type=int, metavar='BITSCORE', nargs=1, default=[50],
                    help='Minimum acceptable bitscore for LAST homology search')
parser.add_argument('-C', dest='create_work_dir', action='store_true', default=False, help='Create working directory')
args = parser.parse_args()

num_genomes = len(args.prot_files)
if num_genomes < 2:
    print 'Analysis requires conparison between at least 2 organisms.'
    sys.exit(1)

#
# Set up the work directory
#
if not os.path.exists(args.work_dir[0]):
    if args.create_work_dir:
        os.mkdir(args.work_dir[0])
        # To maintain a consistent order between reruns, we create a table of
        # input sequence filenames and assign ids to use throughout analysis.
        # Eg.
        # 0 seq1.faa
        # 1 seq2.faa
        with open(os.path.join(args.work_dir[0], SEQFILE_MAP), 'w') as seq_h:
            for seq_idx, seq_fn in enumerate(args.prot_files):
                seq_h.write('{0}\t{1}\n'.format(seq_idx, seq_fn))

seqfile_map = {}
try:
    with open(os.path.join(args.work_dir[0], SEQFILE_MAP), 'r') as seq_h:
        for row in csv.reader(seq_h, delimiter='\t'):
            if len(row) != 2:
                raise IOError('Invalid number of columns in file index')
            seqfile_map[row[0]] = row[1]
except IOError as ex:
    print ex.message
    print 'Failed to open file index {0}'.format(os.path.join(args.work_dir[0], SEQFILE_MAP))
    sys.exit(1)

# Extract the basename of protein files to act as organism names.
# We're being a little tricky here and removing any file extension and internal whitespace.
orgs = [os.path.splitext(os.path.basename(pf))[0].replace(' ', '') for pf in args.prot_files]


# Tabulate the proteins across each genome and their descriptions.
print 'Building a list of sequences across the given genomes...'
seq_info_table = {}
isolate_genes = {}
for i in range(num_genomes):
    seq_info_table[orgs[i]] = {}
    # start with all genes as isolated
    isolate_genes[orgs[i]] = set()
    with open(args.prot_files[i], 'r') as fasta_h:
        for rec in FastaIterator(fasta_h, generic_protein, fasta_title2ids_simple):
            seq_info_table[orgs[i]][rec.id] = rec.description
            isolate_genes[orgs[i]].add(rec.id)
    print "Found {0} genes in {1}".format(len(isolate_genes[orgs[i]]), orgs[i])

# Make databases
print 'Database creation ...'
for i, fname in enumerate(args.prot_files):
    db_name = os.path.join(args.work_dir[0], os.path.basename(fname))
    if os.path.exists(db_name + ".prj"):
        print '{0}/{1} database already exists for {2}, skipping.'.format(i+1, num_genomes, fname)
        continue
    print '{0}/{1} creating database for {2}'.format(i+1, num_genomes, fname)
    subprocess.call(['lastdb', '-p', db_name, fname])

num_compare = num_genomes**2 - num_genomes
print 'About to perform {0} pair-wise comparisons.'.format(num_compare)

print 'Pair-wise Searching ...'
n = 0
for i in range(num_genomes):
    for j in range(num_genomes):
        if i == j:
            continue
        n += 1
        query = args.prot_files[i]
        subject = args.prot_files[j]
        out_name = os.path.join(args.work_dir[0], '{0}-{1}'.format(i, j))
        if os.path.exists(out_name):
            print '{0}/{1} already done {2} against {3}, skipping.'.format(n, num_compare, query, subject)
            continue
        print "{0}/{1} searching {2} against {3}".format(n, num_compare, query, subject)
        subject_db = os.path.join(args.work_dir[0], os.path.basename(subject))
        with open(out_name, 'w') as out_h:
            subprocess.call(['lastal', '-T', '1', '-f', '0', '-e', str(args.bitscore[0]), subject_db, query], stdout=out_h)

print "Finding best hits ..."
g = nx.Graph()
for i in range(num_genomes):
    for j in range(num_genomes):
        if i == j:
            continue

        result_name = os.path.join(args.work_dir[0], '{0}-{1}'.format(i, j))

        # read in Last search result for query/subject pair.
        # keep the highest scoring hit for each pairwise result
        with open(result_name, 'r') as result_h:
            best_hit = {}
            best_score = {}
            for line in result_h:
                line = line.strip()
                if line.startswith('#'):
                    continue
                tok = line.split()
                score = int(tok[0])
                subject_id = tok[1]
                query_id = tok[6]
                if not query_id in best_hit or best_score[query_id] < score:
                    best_hit[query_id] = subject_id
                    best_score[query_id] = score

        # insert edges into graph from best-hit result.
        # graph is undirected and duplicate edges increment their weight.
        # weights must exceed 1 to be considered reciprocal.
        for query_id in best_hit:
            v1 = Gene(orgs[i], query_id)
            v2 = Gene(orgs[j], best_hit[query_id])

            # string representation of nodes -- this is way faster but has limitations!
            #v1 = '{0}-{1}'.format(orgs[i], query_id)
            #v2 = '{0}-{1}'.format(orgs[j], best_hit[query_id])

            if g.has_edge(v1, v2):
                g[v1][v2]['weight'] += 1
            else:
                g.add_edge(v1, v2, weight=1)

print 'Raw graph contains {0} nodes and {1} edges.'.format(g.order(), g.size())

print 'Pruning edges for non-reciprocal (one-directional) hits ...'
for v1, v2, d in g.edges_iter(data=True):
    if d['weight'] > 2:
        raise RuntimeError('Edge ({0},{1} weight:{2}) error, weighting greater than 2.'.format(v1, v2, d['weight']))
    if d['weight'] != 2:
        g.remove_edge(v1, v2)

print 'Graph contained {0} nodes {1} edges after pruning.'.format(g.order(), g.size())

print 'Removing self-loops ...'
g.remove_edges_from(g.selfloop_edges())

print 'Removing isolated nodes ...'
g = nx.k_core(g, 1)

print 'Final graph contains {0} nodes and {1} edges'.format(g.order(), g.size())

#print 'Writing edge list'
#nx.write_edgelist(g, args.output[0] + '.edges')

subgraphs = list(nx.connected_component_subgraphs(g))
print 'There are {0} connected components'.format(len(subgraphs))
if args.decompose:
    print 'Decomposing connected components above {0} modularity threshold'.format(args.max_mod)
    decomp_subs = []
    for sg in subgraphs:
        decomp_subs.extend(decompose_graph(sg, args.max_mod))
    subgraphs = decomp_subs
    print 'After decomposition, there are {0} connected components'.format(len(subgraphs))

membership = {}
print 'Writing connected components as ortholog groups.'
with open(args.output[0] + '.txt', 'w') as output_h:
    ortho_wr = csv.writer(output_h, delimiter='\t')
    sg_order_hist = {}

    # sort subgraphs by decreasing number of nodes (order)
    subgraphs = sorted(subgraphs, key=lambda sg: sg.order(), reverse=True)

    for i, sg in enumerate(subgraphs):
        nodes = sorted(nx.nodes(sg))
        num_nodes = sg.order()
        if num_nodes in sg_order_hist:
            sg_order_hist[num_nodes] += 1
        else:
            sg_order_hist[num_nodes] = 1
        ortho_wr.writerow([i+1, num_nodes] + nodes)

        for n in nodes:
            isolate_genes[n.organism].remove(n.name)
            if n.organism not in membership:
                membership[n.organism] = []
            membership[n.organism].append([n.name, i+1])

    print '\nHistogram of component order'
    print 'order\tcount'
    for n in sorted(sg_order_hist):
        print '{0}\t{1}'.format(n, sg_order_hist[n])

print 'Writing subgraph info table.'
with open(args.output[0] + '.subgraph', 'w') as subinfo_h:
    subinfo_wr = csv.writer(subinfo_h, delimiter='\t')
    subinfo_wr.writerow(['id', 'size', 'order', 'density', 'eularian', 'binconn', 'modularity'])
    for i, sg in enumerate(subgraphs):
        part = com.best_partition(sg)
        subinfo_wr.writerow([i+1, sg.size(), sg.order(), nx.density(sg),
                             nx.is_eulerian(sg), nx.is_biconnected(sg), com.modularity(part, sg)])

if args.writesubs:
    print 'Writing all subgraphs ...'
    for i, sg in enumerate(subgraphs):
        out = os.path.join(args.work_dir[0], 'og{0}.graphml'.format(i+1))
        nx.write_graphml(sg, out)

iso_count = 0
for org in membership:
    print 'There remained {0} isolate genes in {1}'.format(len(isolate_genes[org]), org)
    with open(org + '.memb', 'w') as memb_h:
        memb_wr = csv.writer(memb_h, delimiter='\t')
        memb_wr.writerow(['gene','og'])
        for gn in membership[org]:
            memb_wr.writerow(gn)
        for ig in isolate_genes[org]:
            iso_count += 1
            memb_wr.writerow([ig, 'iso{0}'.format(iso_count)])
