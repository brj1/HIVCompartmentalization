library(ape)
library(phylobase)

# Computes the generalized correlation coefficient (GCC of Critchlow et al.
# (2000) Mathematical and Computer Modelling for dichotomous Y with distance 
# counting the number of internal nodes  and tests for significance using a 
# permutation test.
# Parameters:
# trait: vector containing trait membership in the order of the tip labels
# tree: phylo object
# rep: number of replicates for permutation test
# Output:
# Levels: number of distinct traits
# Coefficient: the generalized correlation coefficient
# Null_model_median: median GCC of the permutation test
# Null_model_min = minimum GCC of the permutation test
# Null_model_max = maximum GCC of the permutation test
# The p-value of the permutation test
gcc_test <- function(trait, tree, rep = 1000) {
  # function to compute gcc
  gcc <- function(trait, tree.dist) {
    # convert trait to distance
    trait.dist <- do.call(rbind, lapply(trait, function(x) ifelse(trait == x, 0, 1)))
    diag(trait.dist) <- NA
    trait.dist<- trait.dist[!is.na(trait.dist)]
    
    cor(trait.dist, tree.dist, use = 'everything', method = 'pearson')
  }
  
  # compute tree distance
  tree$edge.length <- rep(1, nrow(tree$edge))
  tree.dist <- cophenetic(tree) - 1
  diag(tree.dist) <- NA
  tree.dist<- tree.dist[!is.na(tree.dist)]
  
  # convert trait to integer
  trait <- as.integer(as.factor(trait))
  
  # compute gcc
  tree.stat <- gcc(trait, tree.dist)
  
  # do permutation test
  boot.gcc <- sapply(1:rep, function(i) gcc(sample(trait), tree.dist))
  
  # compute p-value
  p.value <- sum(boot.gcc >= tree.stat) / rep
  
  # output results
  t(data.frame(
    Levels = length(unique(trait)),
    Coefficient = tree.stat,
    Null_model_median = median(boot.gcc),
    Null_model_min = min(boot.gcc),
    Null_model_max = max(boot.gcc),
    p.value
  ))
}

# Wrapper function imitating run_sm copied from https://github.com/prmac/slatkin.maddison
run_gcc <- function (tree, data, trait, rep = 1000)
{
  tree <- drop.tip(tree, tree[["tip.label"]][!tree[["tip.label"]] %in%
                                               data[, 1]])
  data <- data[data[, 1] %in% tree[["tip.label"]], ]
  tr4 <- phylo4d(tree, data, label.type = "column")
  pruned <- as(extractTree(tr4), "phylo")
  ldata <- tdata(tr4, "tip")
  ldata <- ldata[pruned[["tip.label"]], trait]
  result <- gcc_test(ldata, pruned, rep = rep)
  result
}