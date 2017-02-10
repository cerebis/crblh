#!/usr/bin/env Rscript
library(pheatmap)
library(ape)

args = commandArgs(TRUE)
if (length(args) != 3) {
  stop("Usage: [occupancy] [heatmap.pdf] [tree.pdf]\n")
}

d = read.csv(args[1], sep=' ', header=TRUE, row.names = 1)
dbinary = apply(d, 2, function(x){as.numeric(x > 0)})
rownames(dbinary) = rownames(d)
h = pheatmap(t(dbinary),treeheight_row = 0, treeheight_col = 0,
            legend = FALSE, clustering_method = 'complete',
            show_colnames = FALSE, filename=args[2], width=15, height=10)

pdf(file=args[3])
plot(as.phylo(h$tree_row))
dev.off()
