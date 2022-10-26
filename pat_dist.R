## Install/load the following libraries ###

library(optparse)

library(ape)

library(tidyverse)

library(parallel)

library(adephylo)


op <- OptionParser()

op <- add_option(op, "--trees", type='character')

op <- add_option(op, "--info", type='character')

op <- add_option(op, "--output", type='character')

op <- add_option(op, "--threads", type='numeric', default=2)

args <- parse_args(op)


tree.file <- args$trees

info.file <- args$info

output.file <- args$output

threads <- args$threads


#### step 1: Upload trees file ####

trees <- read.tree("trees.nwk")


if (class(trees) == "phylo") {
  
  trees <- list(trees)
  
  names(trees) <- "tree"
  
}


#### step 2: Upload info file containing sequence categories and compute patristic distance ####

info <- read.csv("info.csv", stringsAsFactors=F)

info <- info[match(trees[[1]]$tip.label, info$FULLSEQID), ]

info.s <- split(info, info$TYPE)


ttd.data <- mclapply(
  
  1:length(trees),
  
  function(i) {
    
    tree <- trees[[i]]
    
    
    
    ttd <- sapply(
      
      info.s,
      
      function(x) {
        
        subtree <- keep.tip(tree, x$FULLSEQID)
        
        
        
        patristic.distance <- cophenetic(subtree)
        
        
        
        mean(patristic.distance)
        
      }
      
    )
    
    
    
    data <- data.frame(Treeid=names(trees)[i], Type=names(info.s), TTD=ttd)
    
    write.csv(data, paste(output.file, names(trees)[i], "csv", sep="."), row.names=F)
    
    
    
    data
    
  },
  
  mc.cores=threads
  
) %>% 
  
  do.call(rbind, .)


write.csv(ttd.data, "TTD.csv", row.names=F)