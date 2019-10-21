# CRBLH

## Protein clustering by clustered reciprocal best LAST hit.

CRBLH is a python tool for clustering proteins (proxy for homologs/ortholog/etc) from an all-vs-all LAST search.

The initial set of associations generated by LAST are used to construct a graph, where nodes are protein sequences and edges are LAST hits. The graph is then pruned to only reciprocal associations.

Here, the disconnected components of the graph can then taken as the initial prediction of protein clusters. Optionally, these components can be further decomposed by application of the Louvain method, where a "best partitioning" of the graph is used to predict communities of tightly inter-associated sub-graphs. This decomposition is iteratively applied, until a user specified modularity threshold is reached. The set of decomposed sub-graphs are then taken as the predicted homology groups.

## Prerequisites:
Binaries
- LAST

Python Modules
- biopython
- networkx
- python-louvain
- scipy
- numpy
- tqdm 

## Installation

It is expected that users have previously installed LAST and that it is accessible on the path.

The simplest way to satisfy the python dependencies here is to use Pip. Assuming you have write access to the default installation of Python on your system, the following command would install on the dependent modules and any internal dependencies they may require.

```bash
pip install --upgrade -r requirements.txt
```

If you are using a shared system and the system-level Python installation, you will need to ask your administrator to add these modules or install a personal copy of Python. Unfortunately, at present we have not set-up crblh to maintain its own references to the needed modules.

## Clustering

At present, crblh has a very simple interface and can be a little tedious to start if a user is not fluent in the Unix commandline. It is expect that the sequences are translated into protein.

### A Basic Run

A basic run, where the initial disconnected components of the graph are assumed to be the clusters. The work directory ```work``` will be created (```-C```) and this will involve the sequences for genomes A, B and C.

From the crblh directory:
```bash
./crblh.py -C work myresults A.faa B.faa C.faa
```

### Decomposing Components
You can elect to decompose the components from the beginning, or run crblh a second time. The LAST search results will be re-used.

From the crblh directory:
```bash
./crblh.py --decomp work myresults_decomp A.faa B.faa C.faa
```

### Running Larger Jobs
In some cases, users may wish to analyse more than a handful of genomes. In this case, the basic UI we currently provide means a little work in list all the genome sequence files. To avoid this, we suggest placing all the targets in subdirectory and using Unix pipes and xargs to build the commandline for you.

Assuming I have ~20 genomes whose protein multi-fasta all end in .faa and which I have placed in the subdirectory ```input_prot```. The following will supply all protein files to crblh without the user having to explicitly type them all out. The list is sorted on input just for maintaining consistency.

From the crblh directory:
```bash
ls input_prot/*.faa | sort | xargs ./crblh.py -C work myresults
```

## Output files

  1. Cluster table: Each line of *myresults.txt* represents one cluster with the format "id, size, seq1, seq2, ... seqN". The sequence names are in the format "genome-seqid". The genome name is taken from the supplied file name. It is easiest if genome files does not use hyphens, as within the helper application ```cluster_breakdown.py``` the default strategy for splitting apart the genome name and sequence id is to find the first hyphen. If this is not possible, you can resort to specifying a regular expression.
  2. Graph statistics: *myresults.subgraph* contains information on each cluster from a graph perspective. The table lists cluster id, followed by graph size (number of edges), order (number of nodes/sequences) and density in the columns 1-4. Columns 5-6 simply test whether the graph for a cluster is Eulerian or biconnected, while column 7 lists the Louvain modularity for the cluster. When modularity is near 0, there is little structure suggesting a graph possesses more than a single community (no clumps). As a rule of thumb, a confident prediction of a protein cluster might be one where the graph is biconnected and with modularity=0.

## Primary Result Breakdown

Users can compile a set of sequences and a breakdown of sequence annotations using the tool ```cluster_breakdown.py```. By default, the tool with do this for the entire set of sequences. 

```bash
./cluster_breakdown.py myresults.txt outdir A.faa B.faa C.faa 
```

For large jobs, the same method can be applied as was explained above.

```bash
ls input_prot/*.faa | sort | xargs ./cluster_breakdown.py myresults.txt outdir
```

The result will be a directory containing 2 files for each cluster (clxx.faa, clxx.csv), containing the protein sequences belonging to the cluster and a table cataloging the various annotations given to each sequence. If the user has supplied sequences without annotation, the second file will not be of much value. If annotations were supplied, then this file is a handy means of infering function without resorting to additional searches or database lookups.

In any case, the columns are: 
  1. clid
  2. fraction (%) of genes with the given annotation
  3. number of genes observed with the given annotation
  4. the annotation in question

## Occupancy Table

An occupancy table can easily be created from the per-genome membership files using ```occupancy.py```. One use of this table is to consider which genomes share blocks of protein clusters.

## Hierarchical clustering of protein occupancy

The script ```heatmap.R``` is a simple tool to generate a heatmap, sorting both rows and columns using hierarchical clustering. In this case, the script expects the occupancy table generated both. The result is a heatmap, depicting strict "yes/no" occuapancy of each genome vs each protein cluster. This image is helpful in gaining a sense of the genes shared across the comparison set. This might also be regarded as a qualitative estimation of the "core" genes.

Required packages:
- pheatmap
- ape

Usage:
```bash
./heatmap.R occupancy.csv heatmap.pdf tree.pdf
```

Users are encourage to use this script as an example of what can be done. Rather than producing only a figure, the results of hierarchical clustering can be interogated to, for instance, select out sets of protein clusters or genomes with particular relatedness.


