# Bioinformatics Toolkit Cheat Sheet

Welcome to the Bioinformatics Toolkit Cheat Sheet! This guide provides a quick overview of essential data types, formats, and tools used in bioinformatics. It's designed to help beginners and experts alike navigate the complex world of bioinformatics with ease.

## Data Types & Formats

- **FASTA**  
  Used for nucleotide and protein sequences. Simple format with a '>' symbol preceding the header line, followed by the sequence data.  
  [More Info](https://en.wikipedia.org/wiki/FASTA_format)

- **FASTQ**  
  Stores sequences and their quality scores. Used in high-throughput sequencing.  
  [More Info](https://en.wikipedia.org/wiki/FASTQ_format)

- **SAM/BAM**  
  Sequence Alignment/Map format; BAM is the binary version of SAM. Used for storing large nucleotide sequence alignments.  
  [More Info](https://samtools.github.io/hts-specs/)

- **VCF (Variant Call Format)**  
  Used to describe gene variants, such as SNPs and insertions/deletions, in genomic data.  
  [More Info](https://samtools.github.io/hts-specs/VCFv4.2.pdf)

- **GFF/GTF**  
  General Feature Format/ Gene Transfer Format used to describe gene structures.  
  [More Info](https://www.ensembl.org/info/website/upload/gff.html)

## Tools & Software

### Sequence Analysis

- **BLAST (Basic Local Alignment Search Tool)**  
  For aligning query sequences against a database.  
  [Visit BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)

- **Bowtie/Bowtie2**  
  Widely used for aligning sequencing reads to long reference sequences.  
  [Visit Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

### Multiple Sequence Alignment (MSA)

- **Clustal Omega/ClustalW**  
  Popular tools for multiple sequence alignment.  
  [Visit Clustal](http://www.clustal.org/)

- **MAFFT**  
  Another widely used tool for sequence alignments, known for its speed and efficiency.  
  [Visit MAFFT](https://mafft.cbrc.jp/alignment/software/)
  
### Sequence Clustering and Searching

- **MMseqs2 (Many-against-Many sequence searching)**  
  Fast and sensitive searching and clustering of large protein sequence datasets.  
  [Visit MMseqs2](https://github.com/soedinglab/MMseqs2)

### Structural Biology

- **PDB (Protein Data Bank)**  
  Stores 3D structures of molecules; widely used for proteins, nucleic acids, and complex assemblies.  
  [Visit PDB](https://www.rcsb.org/)

- **PyMOL**  
  An open-source molecular visualization tool.  
  [Visit PyMOL](https://pymol.org/2/)

- **AlphaFold**  
  A deep learning approach to protein structure prediction, developed by DeepMind.  
  [Visit AlphaFold](https://github.com/deepmind/alphafold)  
  Use ColabFold for an accessible implementation of AlphaFold online.  
  [Use ColabFold](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/beta/AlphaFold2_advanced.ipynb)

- **Chimera**  
  An extensible program for interactive visualization and analysis of molecular structures and related data.  
  [Visit Chimera](https://www.cgl.ucsf.edu/chimera/)

### Genome Analysis

- **GenBank**  
  The NIH genetic sequence database, an annotated collection of all publicly available DNA sequences.  
  [Visit GenBank](https://www.ncbi.nlm.nih.gov/genbank/)
  
- **BEDTools**  
  A powerful tool for genome arithmetic; it allows one to compare, manipulate, and annotate genomic features in BED format.  
  [Visit BEDTools](https://bedtools.readthedocs.io/en/latest/)

- **GATK (Genome Analysis Toolkit)**  
  Widely used for variant discovery in high-throughput sequencing data.  
  [Visit GATK](https://gatk.broadinstitute.org/hc/en-us)

- **SAMtools**  
  Provides various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing, and generating alignments in a per-position format.  
  [Visit SAMtools](http://www.htslib.org/)

### Phylogenetics

- **MEGA (Molecular Evolutionary Genetics Analysis)**  
  A comprehensive tool for comparative genomics, especially phylogenetics.  
  [Visit MEGA](https://www.megasoftware.net/)

- **PhyML**  
  A tool for maximum likelihood phylogeny estimation from DNA or amino acid sequence data.  
  [Visit PhyML](http://www.atgc-montpellier.fr/phyml/)

### Bioinformatics Pipelines and Workflows

- **Galaxy**  
  An open, web-based platform for accessible, reproducible, and transparent computational research.  
  [Visit Galaxy](https://usegalaxy.org/)

- **Bioconductor**  
  An open-source project that provides tools for the analysis and comprehension of high-throughput genomic data.  
  [Visit Bioconductor](https://www.bioconductor.org/)

