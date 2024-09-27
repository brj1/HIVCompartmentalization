# HIVCompartmentalization

Scripts for:

Shahid, Aniqa, et al. (2022) "HIV proviral genetic diversity, compartmentalization and inferred dynamics in lung and blood during long-term suppressive antiretroviral therapy", *PLoS Pathogens*, in press.

Files:
* epa_root.R: EPA outgroup rooting
* FaithPD_not_normalized.R: Faith's phylogenetic diversity
* gcc_test.R: correlation coefficient (CC) test
* pat_dist.R: patristic distance
* sample_data: sample data for scripts

Documentation:
epa_root.R is used to find the root placement of a tree, plasmatree, by comparing to an EPA rooting of plasmatree, alltree, with outgroup, outgroup. The result is a rooted tree, rootedtree, with the same root position as alltree but the pranch lengths of plasmatree. EPA can be performed in RAxML.

`Rscript epa_root.R --alltree=alltree --plasmatree=plasmatree --rootedtree=rootedtree --outgroup=outgroup`

alltree: EPA rooted phlyogentic tree file in newick format with taxa of plasmatree plus outgroup
plasmatree: a phylogenetic tree file in newick format
rootedtree: outout phylgenetic tree of plasmatree rooted at the outgroup og alltree
outgroup: outgroup taxon name in alltree

---

FaithPD_not_normalized.R computes Faith's phylogenetic diversity (Faith (1992) Biological Conservation) of a tree over each group of taxa specified.

`Rscript Faith_PD_nor_normalized --trees=trees --info=info --output-output [--threads=threads]`

trees: phylogenetic tree file in newick format or multiline phylogenetic tree file in newick format
info: info CSV file (CSV with columns FULLSEQID (taxa) and TYPE (taxon group)
output: output CSV file with tree name, taxa type, and PD (phylogenetic diversity)
threads: number of threads for multithrading tree caluctations (default 2)

---

part_dist.R computes the mean patristic distance of pairs of sequences in a tree over each group of taxa taxa

`Rscript pat_dist --trees=trees --info=info --output=output [--threads=threads]`

trees: phylogenetic tree file in newick format or multiline phylogenetic tree file in newick format
info: info CSV file (CSV with columns FULLSEQID (taxa) and TYPE (taxon group)
output: output CSV file with tree name, taxa type, and TTD (mean patristic distance
threads: number of threads for multithrading tree caluctations (default 2)

---

gcc_test.R contains functions to perform corellation coefficient compartmentalization testing (Critchlow et al (2000) Math Comput Model). They are modelled after the Slatkin Maddisin tests in the repository hhtps://github.com/prmac/slatkin.maddison.

To run the compartmentalization test:

```R
source("gcc_test.R") # load in scripts

tree <- read.tree("sample_data/tree.nwk") # newick tree file

data <- read.csv("sample_data/info.csv") # CSV file (first column containing taxa names, another column containing groups)

gcc <- run_gcc(tree, data, trait = "TYPE", rep = 1000) # set trait to group column header, rep is the number of replicates for the permutation test

write.csv(t(gcc), "sample_data/gcc.csv")  # write results to file
```

run_gcc returns a table with
Levels: mnumber of groups
Coefficient: correlation coefficient (CC)
Null_model_median: median null CC,
Null_model_min = minimum null CC,
Null_model_max: maximum null CC,
p.value: permutation test p value
