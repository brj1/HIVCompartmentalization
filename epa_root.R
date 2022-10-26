library(ape)
library(phytools)
library(optparse)
library(dplyr)
library(tidyr)
library(tibble)

args <- OptionParser() %>%
  add_option("--alltree", type = 'character') %>%
  add_option("--plasmatree", type = 'character') %>%
  add_option("--rootedtree", type = 'character') %>%
  add_option("--outgroup", type = 'character') %>%
  parse_args()

all.tree.file <- args$alltree
outgroup <- args$outgroup
plasma.tree.file <- args$plasmatree
rooted.tree.file <- args$rootedtree


all.tree <- read.tree(all.tree.file)
plasma.tree <- read.tree(plasma.tree.file)

plasma.tree <- root(plasma.tree, outgroup)
plasma.tree <- drop.tip(plasma.tree, outgroup)

clade.edge <- which(plasma.tree$edge[, 1] == (Ntip(plasma.tree) + 1))
clade.node <- plasma.tree$edge[clade.edge, 2]

clade.len <- sum(plasma.tree$edge.length[clade.edge])

clades <- lapply(
  clade.node,
  function(x) {
    if (x > Ntip(plasma.tree)) {
      extract.clade(x, phy = plasma.tree) %>%
        `$`("tip.label")
    } else {
      plasma.tree$tip.label[x]
    }
  }
)

all.mrca <- sapply(
  clades,
  function(x) {
    if (length(x) == 1) {
      which(all.tree$tip.label == x)
    } else {
      getMRCA(x, phy = all.tree)
    }
  }
)

rep.clade <- which(all.mrca != (Ntip(all.tree) + 1))[1]

node <- all.mrca[rep.clade]
node.edge <- which(all.tree$edge[, 2] == node)

if (all.tree$edge[node.edge, 1] == Ntip(all.tree) + 1) {
  all.len <- sum(all.tree$edge.length[all.tree$edge[, 1] == Ntip(all.tree) + 1])
} else {
  all.len <- all.tree$edge.length[node.edge]
}

len <- plasma.tree$edge.length[clade.edge[rep.clade]] / clade.len * all.len

if (len > all.tree$edge.length[node.edge]) {
  len <- all.len - len
  node <- all.mtca[2]
}

rooted.tree <- reroot(
  all.tree,
  node,
  len
)

write.tree(rooted.tree, rooted.tree.file)
