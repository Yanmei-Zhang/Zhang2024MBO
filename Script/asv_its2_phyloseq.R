# This is code to replicate the analyses and pictures from my DNA extraction method paper 
# *Code developed by Yanmei zhang*

# ITS analysis for DADA2 output with phyloseq

# Load all packages
sapply(c("phyloseq", "Biostrings", "dplyr", "DECIPHER", "msa","phangorn", "ggplot2"), require, character.only = TRUE)

## Import data into phyloseq and build trees (MSI)

# Define meta file path and name
meta_file <- "data/metadata.tsv"
metadata <- read.delim(meta_file, header=T, row.names = 1)

# Define asv table file and name, convert to matrix
asv_table <- "../03_tabletax/dada2_its2_counts.txt"
asvmat <- as.matrix(read.delim(asv_table))

# Define taxa table file and name, convert to matrix
taxa_table <- "../03_tabletax/dada2_its2_taxa_newdatabase1.txt"
taxmat <- as.matrix(read.delim(taxa_table))

# Construct a phyloseq subject
ps.its <- phyloseq(otu_table(asvmat, taxa_are_rows=T),
               sample_data(metadata),
               tax_table(taxmat))

# Store the DNA sequences of our ASVs in the refseq slot of the phyloseq object, and then rename our taxa to a short string
dna <- Biostrings::DNAStringSet(taxa_names(ps.its))
names(dna) <- taxa_names(ps.its)
# Merge refseq into phyloseq object
ps.its <- merge_phyloseq(ps.its, dna)
taxa_names(ps.its) <- paste0("ASV", seq(ntaxa(ps.its)))

# Recover DNA sequence and save in fasta format
rep.seq <- refseq(ps.its)
Biostrings::writeXStringSet(rep.seq,  file = "results/rep_seqs.fasta", format = "fasta")

# Construct phylogenetic tree
# sequence alignment
mult <- msa(rep.seq, method="ClustalW", type="dna", order="input")
phang.align <- as.phyDat(mult, type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

# Merge phylo tree into phyloseq
ps.its <- merge_phyloseq(ps.its, phy_tree(fitGTR$tree))

save.image(file="asv_its2_phyloseq.RData")
saveRDS(ps.its, "results/ps.its.rds")

