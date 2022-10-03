library(phangorn)

a <- read.tree("IE2011_Cognates_rel_ANNOT.nwk")
a
write.tree(unroot.phylo(a), file="test3.nwk")
