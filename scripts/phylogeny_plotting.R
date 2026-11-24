library("ape")
library("ggtree")
library("tidyverse")
library("tidytree")

# Code for generating phylogenies using in Figure 1 and Supplemental Figure 3.

# If you use ggtree in published research, please cite the most appropriate paper(s):
#   
# 1. Guangchuang Yu. Using ggtree to visualize data on tree-like structures. Current Protocols in Bioinformatics, 2020, 69:e96. doi:10.1002/cpbi.96
# 2. Guangchuang Yu, Tommy Tsan-Yuk Lam, Huachen Zhu, Yi Guan. Two methods for mapping and visualizing associated data on phylogeny using ggtree. Molecular Biology and Evolution 2018, 35(12):3041-3043. doi:10.1093/molbev/msy194
# 3. Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam. ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods in Ecology and Evolution 2017, 8(1):28-36. doi:10.1111/2041-210X.12628

setwd("~/data/uce_phylogeny/")

tiptable <- read.csv("species_table.txt", header=F)

#how many loci in phyparts analysis?
ngenes <- 2402

#Renaming function
rename_tips <- function(tree) {
  tree$tip.label <- tiptable[[1]][match(tree$tip.label, tiptable[[2]])]
  return(tree)
}

IQ_species.tree <- read.tree("species_tree/voles_iq_species.treefile")
IQ_species.tree <- rename_tips(IQ_species.tree)
IQ_bstrp <- tibble(node=1:Nnode(IQ_species.tree) + Ntip(IQ_species.tree), bootstrap=IQ_species.tree$node.label)
IQ_species.tree2 <- full_join(IQ_species.tree, IQ_bstrp, by="node")

IQ_species.tree3 <- root(IQ_species.tree, outgroup = "Arvicola amphibius", 
                         resolve.root=T, edgelabel=T)

s1 <- ggtree(IQ_species.tree3) +
  geom_tiplab(fontface=3) +
  geom_treescale(linesize=2) +
  #geom_text(aes(label=bootstrap), hjust=1.1, vjust=-0.7) +
  hexpand(0.3)
s1

#Version to drop/rename additional oregoni and longicaudus sequences
IQ_species.tree4 <- drop.tip(IQ_species.tree3, c("Microtus oregoni (m)", "Microtus longicaudus V1"))
IQ_species.tree4$tip.label <- sub("\\(f\\)", "", IQ_species.tree4$tip.label)
IQ_species.tree4$tip.label <- sub("V2", "", IQ_species.tree4$tip.label)


s2 <- ggtree(IQ_species.tree4, size=1.4) +
  geom_tiplab(fontface=3, size=10) +
  geom_treescale(linesize=2, fontsize=9, x=0, y=13) +
  #geom_text(aes(label=bootstrap), hjust=1.1, vjust=-0.7) +
  hexpand(0.3)
s2

ggsave(filename = "species_tree.png", device = "png", units = "in", width = 12, height = 12)

#Showing species for which we have genomes
IQ_species.tree5 <- drop.tip(IQ_species.tree4, c("Microtus richardsoni", "Microtus obscurus",
                                                 "Microtus agrestis", "Alexandromys fortis",
                                                 "Alexandromys oeconomus", "Neodon fuscus",
                                                 "Chionomys nivalis", "Arvicola amphibius"))

tipcols <- data.frame(sp = IQ_species.tree5$tip.label,
                      tipcol = c("B", "A", "A", "A",
                              "A", "A", "A"))

s3 <- ggtree(IQ_species.tree5, size=1.4) +
    geom_treescale(linesize=2, fontsize=9, x=0, y=6) +
  #geom_text(aes(label=bootstrap), hjust=1.1, vjust=-0.7) +
  hexpand(0.35)
s3 %<+% tipcols + 
  geom_tiplab(aes(colour = factor(tipcol)), fontface=3, size=12) +
  scale_color_manual(values= c(A="#006633", B="black")) +
  guides(colour="none")

ggsave(filename = "species_tree_reduced.png", device = "png", units = "in", width = 14, height = 12)


#Concordance tree
scf.tree <- read.tree("concordance_factors/voles_iq_sCF.cf.tree")
scf.tree <- rename_tips(scf.tree)
scf.tree$node.label <- gsub("100/100/","", scf.tree$node.label)

gcf.tree <- read.tree("concordance_factors/voles_iq_gCF.cf.tree")
gcf.tree <- rename_tips(gcf.tree)
gcf.tree$node.label <- gsub("100/100/","", gcf.tree$node.label)

concord_factors <- tibble(node=1:Nnode(scf.tree) + Ntip(scf.tree), bootstrap=paste(scf.tree$node.label, "/", gcf.tree$node.label, sep = ""))
concord_factors[1,2] <- ""
concord.tree <- full_join(scf.tree, concord_factors, by="node")
concord.tree <- drop.tip(concord.tree, c("Microtus oregoni (m)", "Microtus longicaudus V1"))
concord.tree@phylo$tip.label <- sub("\\(f\\)", "", concord.tree@phylo$tip.label)
concord.tree@phylo$tip.label <- sub("V2", "", concord.tree@phylo$tip.label)
concord.tree@phylo <- root(concord.tree@phylo, outgroup = "Arvicola amphibius")

cf1 <- ggtree(concord.tree) +
  geom_tiplab(fontface=3, size=9) +
  geom_treescale(x=0, y=15, linesize=2, fontsize=5) +
  geom_text(aes(label=bootstrap), hjust=1.1, vjust=-0.7, size=5) +
  hexpand(0.3)
cf1

ggsave(filename = "species_tree_concordance.png", device = "png", units = "in", width = 12, height = 12)

#Phyparts tree
phypart_tree <- read.tree("phyparts/murine_phyparts.concon.tre")
phypart_tree <- lapply(phypart_tree, rename_tips)
class(phypart_tree) <- "multiPhylo"

nconcord <- as.numeric(phypart_tree[[1]]$node.label)
nconflict <- as.numeric(phypart_tree[[2]]$node.label)

phypart_tree2 <- as.treedata(phypart_tree[[1]])

concon <- tibble(node=1:Nnode(phypart_tree2) + Ntip(phypart_tree2), concord = nconcord, conflict = nconflict)
phypart_tree2 <- full_join(phypart_tree2, concon, by="node")

phypart_hist <- read.table("phyparts/murine_phyparts.hist", header = FALSE, sep = ",", 
                           fill = TRUE, row.names = 1)
phypart_hist[1,][!is.na(phypart_hist[1,])]

phypart_data <- data.frame(matrix(ncol = 5, nrow = 16))
colnames(phypart_data) <- c("node", "nconcord", "nconflict_top", "nconflict_rest", "no_supp")
phypart_data$node <- Ntip(phypart_tree2)+2:Nnode(phypart_tree2)
phypart_data$nconcord <- nconcord[2:length(nconcord)]
phypart_data$no_supp <-  rep(ngenes, times = length(nconcord)-1) - nconcord[2:length(nconcord)] - nconflict[2:length(nconflict)]

for (i in 1:nrow(phypart_hist)) {
  tempdata <- phypart_hist[i,][!is.na(phypart_hist[i,])]
  if (length(tempdata) == 2) {
    phypart_data[i,3] <- 0
    phypart_data[i,4] <- 0
  }
  else if (length(tempdata) == 3) {
    phypart_data[i,3] <- tempdata[2]
    phypart_data[i,4] <- 0
  }
  else {
    phypart_data[i,3] <- tempdata[2]
    phypart_data[i,4] <- nconflict[i+1] - tempdata[2]
  }
}

pies <- nodepie(phypart_data, cols = 2:5, alpha = 0.9, 
                color = c("#1C6AA8", "red", "#65AE5F", "gray"))

s3 <- ggtree(phypart_tree2) + 
  geom_tiplab(fontface=3) +
  geom_text(aes(label=concord, x=x-0.9), hjust="left", vjust=-0.6) +
  geom_text(aes(label=conflict, x=x-0.9), hjust="left", vjust=1.2) +
  hexpand(0.2)
s3

inset(s3, pies, width = 0.08, height = 0.08)

ggsave(filename = "manuscript/figures/figureS3.png", device = "png", units = "in", width = 8, height = 9)

