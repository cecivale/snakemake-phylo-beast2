# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  snakemake-phylo-beast2
#        V\ Y /V    Plot phylogenetic trees
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------


library(tidyverse)
library(lubridate)
library(treedataverse)
library(ggsci)
library(scales)

tree_file <- snakemake@input[["tree"]]
if (grepl("nexus", tree_file)) {
  tree_ref <- treeio::read.nexus(tree_file)
} else {
  tree_ref <- treeio::read.tree(tree_file)
}
tree <- drop.tip(tree_ref, snakemake@params[["outgroup"]])

ids <- read_tsv(snakemake@input[["ids"]])

d <- as_tibble(tree) %>% left_join(ids, by = c("label" = "sample_id"))

p <- ggtree(as.treedata(d), size = 0.4, 
            #aes(color = deme))  +
            color = "grey70")  +
  geom_tippoint(aes(color = deme), size = 1) +
  geom_nodepoint(aes(subset = as.numeric(label) >= 95), size = 1, shape = 15, color = "grey30") +
  scale_size_continuous(range = c(0, 2), name = "bootstrap") +
  scale_alpha(name = "bootstrap") +
  scale_color_manual(values = c("#ffd900", "#1464b5"), na.value = "grey70") +
  # theme_tree2()
  geom_treescale(x = 0, y = 220, label = "subs/site", fontsize = 3) 

ggsave(p, filename = snakemake@output[["fig"]], width = 297, heigh = 210, units = "mm")
