# Taxonomy
library("Rcpp")
library("dada2")
library("ShortRead")
library("Biostrings")
### File Paths
path <- "/lustre/project/svanbael/bolivar/Mimulus_sequences/Aponte_8450_23052601"
out_dir <- "/lustre/project/svanbael/bolivar/Mimulus_sequences" #An "out" directory to avoid crowding main directory of new files. Also keeps intact raw files.

set.seed(755) #random number generator for reproducibility
seqtab.nochim <- readRDS(file.path(out_dir, "ASV_tables", "/ASV_seqtab_nochim.rds"))
unite.ref <- file.path(out_dir, "Taxonomy", "sh_general_release_dynamic_29.11.2022.fasta")  #DOI: 10.15156/BIO/2483911 # CHANGE ME to location on your machine
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, minBoot = 50, tryRC = TRUE) #Multithread = FALSE in Windows. TRUE in Mac/Linux.

# Inspecting the taxonomic assignments:
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.csv(taxa.print, file.path(out_dir, "Taxonomy", "/assign_tax_mim2.csv"))

#Done