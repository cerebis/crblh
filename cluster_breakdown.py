#!/usr/bin/env python
from Bio import SeqIO
from Bio.Alphabet import generic_protein
from Bio.Alphabet import generic_dna
from collections import Counter
import tqdm
import csv
import re
import sys
import argparse
import os


def fasta_id(header):
    """
    Implementation of title2ids function for use with Bio.SeqIO.FastaIO.FastaIterator
    Splits FASTA headers at the first whitespace. Word[0] becomes the id, with the remainder
    forming the description.

    :param header: fasta header line
    :return: sequence id
    """
    tok = re.split('\s', header, maxsplit=1)
    return tok[0]


def get_description(seq):
    """
    Return the description text from a Seq object
    :param seq: seq object
    :return: description string or empty string
    """
    m = re.match('^\S+\s(.*)$', seq.description)
    if m is None:
        return ''
    else:
        return m.group(1)


def take_after(s, delim='-'):
    """
    Simple method to split genome-seq string at a character boundary. This happens to get
    broken when users make use of cumbersome names.
    :param s: string to split
    :param delim: delimiter marking the end of the prefix
    :return: the suffix
    """
    return s[s.find(delim) + 1:]


def regex_after(pattern):
    """
    A regular expression parser for genome identifiers. Users can define an expression with a
    single capturing group, which will be used to extract the sequence identifier (only) from
    the genome-seq string.
    :param pattern:
    :return:
    """
    def _call(s):
        matcher = re.compile(pattern)
        if matcher.groups != 1:
            raise RuntimeError('Patterns must contain a single capturing group')
        return matcher.match(s).group(1)
    return _call

def count_lines(hndl):
    """
    Count the lines in a file. For the cluster table, this should equal the number of clusters.
    The file handle is returned to the start of the file.
    :param hndl: open readable file handle
    :return: number of lines in file.
    """
    n = 0
    while True:
        try:
            hndl.next()
            n += 1
        except StopIteration:
            break
    hndl.seek(0)
    return n



parser = argparse.ArgumentParser(description='Fetch sequences and annotation breakdown for protein cluster(s).')
parser.add_argument('cluster_table', metavar='ORTHO_TABLE', help='Cluster table.')
parser.add_argument('output', metavar='OUTPUT DIR', help='Output directory.')
parser.add_argument('fasta_files', metavar='FASTA_FILE', nargs='+', help='Multi-fasta sequences per organism.')
parser.add_argument('--id', dest='cluster_id', metavar='CLID', help='Find for a single cluster identifier.')
parser.add_argument('-p', '--protein', dest='isprotein', action='store_true', default=False,
                    help='Supplied sequences are translated proteins.')
parser.add_argument('--tmpdb', dest='tmpdb', default=':memory:',
                    help='Temporary sequence database [default in memory].')
parser.add_argument('--pattern', metavar='REGEX', help='Regex pattern for selecting sequence IDs')

try:
    args = parser.parse_args()

    if os.path.exists(args.output):
        print 'Error: {0} already exists'.format(args.output)
        sys.exit(1)
    else:
        os.mkdir(args.output)

    if args.tmpdb != ':memory:' and os.path.exists(args.tmpdb):
        print 'Error: {0} already exists'.format(args.tmpdb)

    print 'Preparing index ...'
    if args.isprotein:
        print 'Indexing aa seqs {0}'.format(args.fasta_files)
        seq_db = SeqIO.index_db(args.tmpdb, args.fasta_files, 'fasta', generic_protein, fasta_id)
    else:
        print 'Indexing nt seqs {0}'.format(args.fasta_files)
        seq_db = SeqIO.index_db(args.tmpdb, args.fasta_files, 'fasta', generic_dna, fasta_id)

    print 'Indexed {0} sequences'.format(len(seq_db))

    if args.pattern:
        get_id = regex_after(args.pattern)
    else:
        get_id = take_after

    with open(args.cluster_table, 'r') as clust_h:

        total_clusters = count_lines(clust_h)

        cl_table = csv.reader(clust_h, delimiter='\t')

        found_single = False

        for row in tqdm.tqdm(cl_table, "Processing", total=total_clusters):

            desc_tally = []

            clid = row[0]

            if found_single:
                sys.exit(0)

            if args.cluster_id:
                if clid == args.cluster_id:
                    found_single = True
                else:
                    continue

            seq_fname = os.path.join(args.output, 'cl{0}.faa'.format(clid))
            tally_fname = os.path.join(args.output, 'cl{0}.csv'.format(clid))
            with open(seq_fname, 'w') as seq_out, open(tally_fname, 'w') as tally_out:

                for genome_i in row[2:]:
                    # just the sequence identifier, no genome name
                    seq_id = get_id(genome_i)
                    if seq_id not in seq_db:
                        print 'Failed to find sequence {0}'.format(seq_id)
                    else:
                        d = get_description(seq_db[seq_id])
                        desc_tally.append(d)
                        SeqIO.write(seq_db[seq_id], seq_out, 'fasta')

                cd = Counter(desc_tally)
                tally_out.write('Breakdown of protein annotations:\n')
                for desc, count in cd.iteritems():
                    perc_desc = float(count)/len(desc_tally)*100
                    tally_out.write('CLID:{0} {1:.1f}\t{2}\t"{3}"\n'.format(args.cluster_id, perc_desc, count,
                                                                            desc if desc else 'blank'))

    if args.cluster_id and not found_single:
        print 'Did not find a record for CLID:{0} in cluster table'.format(args.cluster_id)

except Exception as e:
    print e

