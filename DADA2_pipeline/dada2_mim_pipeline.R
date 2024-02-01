#DADA2 pipeline
#Modified by Bolívar Aponte Rolón for bioinformatic analyses of ITS amplicon sequences
#14/june/2023

# Loading packages
# Activate commands as needed.
# DADA2 package and associated
#library("usethis") #Doesn't install on the HPC cluster for some reason.
#library("devtools") #Doesn't install on the HPC cluster for some reason.
library("Rcpp")
library("dada2")
library("ShortRead")
library("Biostrings")  #This will install other packages and dependencies.


### File Paths
path <- "/lustre/project/svanbael/bolivar/Mimulus_sequences/Aponte_8450_23052601"
out_dir <- "/lustre/project/svanbael/bolivar/Mimulus_sequences" #An "out" directory to avoid crowding main directory of new files. Also keeps intact raw files.

list.files(path)
fnFs <- sort(list.files(path, pattern = "R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2.fastq", full.names = TRUE))

# Identifying primers
#CHANGE primers accordingly
FWD<-"CACTCTTTCCCTACACGACGCTCTTCCGATCTCTTGGTCATTTAGAGGAAGTAA"# Forward ITS1f_adapt from IDT 
nchar(FWD) #Number of primer nucleotides.
REV<-"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGCTGCGTTCTTCATCGATGC"# Reverse primer ITS2r_adapt from IDT
nchar(REV)

### Verifying the orientation of the primers
allOrients <- function(primer) {# Create all orientations of the input sequence
                   require(Biostrings)
                   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
                   orients <- c(Forward = dna, 
                                Complement = Biostrings::complement(dna), 
                                Reverse = Biostrings::reverse(dna),
                                RevComp = Biostrings::reverseComplement(dna))
                   return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
#FWD2 <- FWD.orients[["Complement"]] #Use if you suspect an orientation mix-up.
#REV2 <- REV.orients[["Complement"]]

### Filter Reads for ambiguous bases (N)
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered forward read files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs)) #Reverse reads
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) #multithread = TRUE on Mac OS, FALSE in Windows

### Checking for primer hits
set.seed(123)
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

### Cutadapt: removal of primers
#Once this has been completed there is no need to run again when working on the script
cutadapt <-  "/home/baponterolon/.conda/envs/virtual_env/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(out_dir, "cutadapted") #Remember where this "out" directory path leads to.
#all.cut <- file.path(out_dir, "FastQC") # Path to concatenated files 
print(path.cut) #Checking if the path is correct.
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

### Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 removes FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) #input files
}

#Re-inspecting if all primers were removed.  
#Once this has been completed there is no need to run again when working on the script
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients,
      primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits,fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
                                                                                                        
#Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2.fastq", full.names = TRUE))
                                          
#Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_R")[[1]][1] #String in commas needs to be updated according to naming convention. If you have multiple underscores in the name then select the underscore next to the "R", like above, or any other unique identifier in the character string.
              
sample.namesF <- unname(sapply(cutFs, get.sample.name))
sample.namesR <- unname(sapply(cutRs, get.sample.name))
          
### Inspect the read quality
#plotQualityProfile(cutFs[20:21])
#plotQualityProfile(cutRs[20:21])

### Filter and trim
filtFs <- file.path(out_dir, "Filtered", basename(cutFs))
filtRs <- file.path(out_dir, "Filtered", basename(cutRs))
names(filtFs) <- sample.namesF
names(filtRs) <- sample.namesR

#Truncating Forward and Reverse reads to 230 and 180 respectively. The reverse maxEE is relaxed due to overall low quality of reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                             truncLen=c(230,180),
                             maxN=0,
                             maxEE=c(2,3),
                             truncQ=2,
                             minLen=50,
                             rm.phix=TRUE,
                             compress=TRUE,
                             multithread=TRUE) # minLen: Remove reads with length less than minLen. minLen is enforced after trimming and truncation. #enforce min length of 50 bp

#Once it is completed there is no need to run again unless you are changing parameters.
#Save file
saveRDS(out, file.path(out_dir, "Preprocess", "/out.rds"))

### Dereplication
# Dereplication combines all identical reads into one unique sequence with a corresponding abundance equal to the number of reads with that unique sequence. It is done because it reduces computation time by eliminating redundancy -- From DADA2 tutorial, v1.8
derepFs <- derepFastq(filtFs, n = 1000, verbose=TRUE) #n prevents it from reading more than 1000 reads at the same time. This controls the peak memory requirement so that large fastq files are supported. 
derepRs <- derepFastq(filtRs, n = 1000, verbose=TRUE)
saveRDS(derepFs, file.path(out_dir, "Preprocess", "/derepFs.rds"))
saveRDS(derepFs, file.path(out_dir, "Preprocess", "/derepFs.rds"))

#derepFS <- readRDS(file.path(out_dir, "Preprocess", "/derepFs.rds"))
#derepRS <- readRDS(file.path(out_dir, "Preprocess", "/derepRs.rds"))

# name the dereplicated reads by the sample names
names(derepFs) <- sample.namesF
names(derepRs) <- sample.namesR

### Learn error rates from dereplicated reads
set.seed(123)
errF <- learnErrors(derepFs, randomize = TRUE, multithread=TRUE) #multithread is set to FALSE in Windows. Unix OS is =TRUE.
errR <- learnErrors(derepRs, randomize = TRUE, multithread=TRUE)

# Save file
saveRDS(errF, file.path(out_dir, "Preprocess", "/errF.rds"))
saveRDS(errR, file.path(out_dir, "Preprocess", "/errR.rds"))

#errF <- readRDS(file.path(out_dir, "Preprocess", "/errF.rds"))
#errR <- readRDS(file.path(out_dir, "Preprocess", "/errR.rds"))

### Plot errors
#plotErrors(errF, nominalQ=TRUE)

### Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
saveRDS(dadaFs, file.path(out_dir, "Preprocess", "/dadaFs.rds"))
saveRDS(dadaFs, file.path(out_dir, "Preprocess", "/dadaFs.rds"))

#dadaFS <- readRDS(file.path(out_dir, "Preprocess", "/dadaFs.rds"))
#dadaRS <- readRDS(file.path(out_dir, "Preprocess", "/dadaRs.rds"))

### Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 20, maxMismatch = 0, verbose=TRUE)
saveRDS(mergers, file.path(out_dir, "Preprocess", "/mergers.rds"))

#mergers <- readRDS(file.path(out_dir, "Preprocess", "/mergers.rds"))
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

### Construct ASV Sequence Table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Save file R object, and .csv
saveRDS(seqtab, file.path(out_dir, "ASV_tables", "/ASV_seqtab.rds")) #Functions to write a single R object to a file, and to restore it.
write.csv(seqtab, file.path(out_dir, "ASV_tables", "/ASV_seqtab.rds"))

# Open from here in case R crashes
#seqtab<- readRDS(file.path(out_dir, "ASV_tables", "/ASV_seqtab.rds"))

### Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread = FALSE, verbose=TRUE) #Multithread = FALSE in Windows

# Frequency of chimeras
sum(seqtab.nochim)/sum(seqtab)

# Save file
saveRDS(seqtab.nochim, file.path(out_dir, "ASV_tables", "/ASV_seqtab_nochim.rds"))
write.csv(seqtab.nochim, file.path(out_dir, "ASV_tables", "/ASV_nochim_denoise_filt.csv")) # Long file name but it indicates this file has gone through all the steps in the pipeline.
#seqtab.nochim <- readRDS(file.path(out_dir, "ASV_tables", "/ASV_seqtab_nochim.rds"))

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab.nochim)))

### Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN),
rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.namesF
head(track)

# Save file
saveRDS(track, file.path(out_dir, "Preprocess", "/track.rds"))

### Taxonomy
set.seed(755) #random number generator for reproducibility
unite.ref <- file.path(out_dir, "Taxonomy", "sh_general_release_dynamic_29.11.2022.fasta")  # CHANGE ME to location on your machine
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, minBoot= 50, tryRC = TRUE) #Multithread = FALSE in Windows. TRUE in Mac/Linux.

#Loading from the files saved. In case it crashes, we start from here.
#seqtab.nochim2 <- readLines(file.path(out_dir, "output.txt")) 
#seqtab.nochim2 <- read.csv(file.path(out_dir, "ASV_tables", "/ASV_nochim_denoise_filt.csv")) 
#seqtab.matrix <- as.matrix(seqtab.nochim2) #assignTaxonomy needs a vector matrix

# Inspecting the taxonomic assignments:
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.csv(taxa.print, file.path(out_dir, "Taxonomy", "/assign_tax_mim2.csv"))

#Done