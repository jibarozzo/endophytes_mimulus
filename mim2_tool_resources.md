# Environment Setup

## Objective

The tools and resources presented here are necessary and useful to carry out the bioinformatic pipeline(s) used in the Van Bael Lab for ITS and 16S sequence data. The main repository can be found [here](https://github.com/VanBaelLab/VBL_DADA2). 
The pipeline below is mainly focused on DADA2 pipeline presented by Benjamin Callahan et al. in its [ITS](https://benjjneb.github.io/dada2/ITS_workflow.html) and [16S](https://benjjneb.github.io/dada2/tutorial.html) variants. The pipeline follows closely the DADA2 pipeline, the work of Mareli Sánchez Juliá for her MAMF project, and Farrer Lab at Tulane University. A hybrid approach with the USEARCH and VSEARCH tools can be employed and then assign taxonomy with the dada2 function `assignTaxonomy`. 

## Resources

The list below is to help guide your search and aid your bioinformatics journey. Each project and sequencing run is different and should be treated as such. Determine the parameters to filter, cut, trim and truncate accordingly.

### Van Bael Files
   + [VBL_Bioinformatics](https://drive.google.com/open?id=1Z4jHQDcS4dOpG6hlkXVWMw72tZeCYYRH&usp=drive_fs) folder in the Google Drive were you can find scripts and notes on bioinfomatic pipelines from past graduate students and post-docs.
   + [Farrer Lab](https://github.com/ecfarrer/LAmarshGradient2/blob/master/BioinformaticsITS.R) repository. They work on similar data sets.
   
### DADA2
   + [DADA2](https://benjjneb.github.io/dada2/index.html) Main page for the DADA2 workflow. Other links are scattered through the document. 
   + [Callahan, B., McMurdie, P., Rosen, M. et al. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581–583 (2016).](https://doi.org/10.1038/nmeth.3869)
   + [Callahan, B., McMurdie, P. & Holmes, S. Exact sequence variants should replace operational taxonomic units in marker-gene data analysis. ISME J 11, 2639–2643 (2017).](https://doi.org/10.1038/ismej.2017.119)
   + [Bioinformatics Cookbook](https://bioinformaticsworkbook.org/dataAnalysis/Metagenomics/Dada2.html#gsc.tab=0)

### Bionformatic software tools and environments 

#### R and Python3

A lot of the bioinformatic pipelines take place in a Unix/Linux environment. If you have Mac OS then half the troubles are gone. If you have Windows OS you will need to install a Virtual Machine to operate a Linux platform (e.g. Ubuntu). Although various Unix-like environments and command-line interfaces for Windows exist they have limitations and not all packages or modules are supported by them.

   + [Virtual Machine installation](https://www.virtualbox.org/)
   + [Ubuntu installation](https://ubuntu.com/tutorials/how-to-run-ubuntu-desktop-on-a-virtual-machine-using-virtualbox#1-overview)
   + [Virtual machine shared folders](https://averagelinuxuser.com/virtualbox-shared-folder/)
   + [Python installation](https://www.python.org/downloads/)
   + [Miniconda3](https://docs.conda.io/en/latest/miniconda.html#installing)
      + [Windows OS](https://www.codecademy.com/article/install-python3)
      + [Mac and Linux OS](https://engineeringfordatascience.com/posts/install_miniconda_from_the_command_line/)
   
   
Python comes installed in most OS and different versions can coexist and are usually employed simultaneously by apps and software. Verify the minimum version needed for your application. FastQC, MultiQC, Bioconductor and other bioinformatics tools use Python3+ (e.g. > 3.7).
If you are not familiar with the Python ecosystem, stop and take a step back. Think about where you want to install these tools. In your virtual machine, laptop, lab computer? What file paths (location in computer)?
Try to think about these things before installing to avoid installing in random places and then not knowing where it is and being unable to call a command. It is best to execute the bioinformatic pipelines in a HPC cluster or local computer. Avoid having the files in the lab's Google Drive and trying to access it this way. It is possible but the shortcuts provided by Google Drive make for long file paths and are not accessible through Virtual Machines.

#### Cutadapt

   + [Cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html) finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads. Sequencing cores usually de multiplex your data but do no remove the primers and adapters. This is the tool for the job.

   +  [cutadapt installation](https://cutadapt.readthedocs.io/en/stable/installation.html#installation-on-windows): You will need to install this application to complete this pipeline. 

#### FIGARO
   + [FIGARO](https://github.com/Zymo-Research/figaro#figaro) "FIGARO will quickly analyze error rates in a directory of FASTQ files to determine optimal trimming parameters for high-resolution targeted microbiome sequencing pipelines, such as those utilizing DADA2 and Deblur."

#### Bioconductor
The Bioconductor project purpose is to develop, support, and disseminate free open source software that facilitates rigorous and reproducible analysis of data from current and emerging biological assays.
     + [Installation](https://www.bioconductor.org/install/)
     
#### USEARCH and VSEARCH search and clustering algorithms
   + [VSEARCH](https://github.com/torognes/vsearch): From their website [...]"supports de novo and reference based chimera detection, clustering, full-length and prefix dereplication, rereplication, reverse complementation, masking, all-vs-all pairwise global alignment, exact and global alignment searching, shuffling, subsampling and sorting. It also supports FASTQ file analysis, filtering, conversion and merging of paired-end reads."
   + [USEARCH](https://www.drive5.com/usearch/): From their website "USEARCH offers search and clustering algorithms that are often orders of magnitude faster than BLAST."
   
#### FastQC and MultiQC
   + [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/: tool for reading the [fastq](https://en.wikipedia.org/wiki/FASTQ_format) file format and extracting the quality scores. 
      + [Installation]: [conda](https://anaconda.org/bioconda/fastqc [video](https://www.youtube.com/watch?v=Umo1pRuT0OI)
      + [How does FastQC work?](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf)
   + [MultiQC](https://multiqc.info/): MultiQC doesn’t do any analysis for you - it just finds results from other tools that you have already run and generates nice reports.

### Other resources
   + [Dr. Rachel Lappan's pipeline and workflow](https://rachaellappan.github.io/16S-analysis/index.html)