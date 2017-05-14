# BTyper
## A computational tool for virulence-based classification of *Bacillus cereus* group isolates using nucleotide sequencing data

## Overview

BTyper is a command-line tool that employs a combination of (i) virulence gene-based typing, (ii) multi-locus sequence typing (MLST), (iii) panC clade typing, and (iv) rpoB allelic typing to rapidly classify *B. cereus* group isolates using nucleotide sequencing data.

The program, as well as the associated databases, can be downloaded from https://github.com/lmc297/BTyper.

Post issues at https://github.com/lmc297/BTyper/issues


### Citation

If you found the BTyper tool, its source code, and/or any of its associated databases useful, please cite:
  
Carroll, Laura M., Jasna Kovac, Rachel A. Miller, Martin Wiedmann. 2017. Rapid, high-throughput identification of anthrax-causing and emetic *Bacillus cereus* group genome assemblies using BTyper, a computational tool for virulence-based classification of *Bacillus cereus* group isolates using nucleotide sequencing data. Submitted to *Applied and Environmental Microbiology*.


------------------------------------------------------------------------
  
  
  ## Quick Start
  
  #### Command structure:
  
```
btyper -t [input data type] -i [input file(s)] -o [output directory] [options...]
```

For help, type `btyper -h` or `btyper --help`

For your current version, type `btyper --version`

#### Sample commands and analyses:

**Using a fasta/multifasta containing 1 or more closed genomes:**
  
```
btyper -t seq -i my_genomes.fasta -o /path/to/output_directory
```

**Using a fasta file containing contigs or scaffolds (concatenates contigs/scaffolds into a pseudochromosome):**
  
```
btyper -t seq -i my_draftgenome.fasta -o /path/to/output_directory --draft_genome
```

**Using ILLUMINA paired-end reads in fastq.gz format (calls SPAdes to assemble):**
  
```
btyper -t pe -i forward_reads.fastq.gz reverse_reads.fastq.gz -o /path/to/output_directory
```

**Using ILLUMINA single-end reads in fastq.gz format (calls SPAdes to assemble):**
  
```
btyper -t se -i illumina_reads.fastq.gz -o /path/to/output_directory
```

**Using ILLUMINA single- or paired-end reads in SRA (sequence read archive) format (calls SPAdes to assemble):**
  
```
btyper -t sra -i illumina_reads.sra -o /path/to/output_directory
```

**Using an SRA accession number corresponding to ILLUMINA single- or paired-end reads (downloads reads from SRA, calls SPAdes to assemble):**
  
```
btyper -t sra-get -i SRAXXXXXXX -o /path/to/output_directory
```

------------------------------------------------------------------------
  
  
  ## Installation
  ### Install BTyper using Homebrew (macOS users)
  
  BTyper and its dependencies can be installed using <a href="https://brew.sh/">Homebrew</a>.

1. First, install Homebrew, if necessary, by running the following command from your terminal:
  
```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

2. Tap Homebrew Science, if necessary, by running the following command from your terminal:
  
```
brew tap homebrew/science
```

3. Tap BTyper by running the following command from your terminal:
  
```
brew tap lmc297/homebrew-btyper
```

4. Install BTyper and its dependencies by running the following command from your terminal:
  
```
brew install btyper
```

### Download and run BTyper using source file (macOS and Ubuntu)

1. To run BTyper, please download and install the following dependencies, if necessary:
  
  <a href="https://www.python.org/downloads/"> Python v. 2.7</a>
  
  <a href="http://biopython.org/DIST/docs/install/Installation.html"> Biopython v. 1.6.9</a>
  
  <a href="https://github.com/Homebrew/homebrew-science/blob/master/blast.rb">BLAST+ v. 2.4.0</a>
  
  <a href="https://github.com/Homebrew/homebrew-science/blob/master/spades.rb">SPAdes v. 3.9.0</a>
  
  <a href="https://github.com/ncbi/sra-tools/wiki/HowTo:-Binary-Installation">SRA toolkit v. 2.8.0</a>
  
  2. Download the BTyper's source file, and store it in your directory of choice:

https://github.com/lmc297/BTyper/raw/master/archive/btyper-1.0.0.tar.gz

3. Extract BTyper program/databases

```
tar -xzvf btyper-1.0.0.tar.gz
```

Note: to ensure that BTyper works correctly, make sure database directories (beginning with "seq_") remain in the same directory as the BTyper executable (stored as "btyper").

4. To run BTyper, call Python 2.7 and supply the full path to the btyper executable:

```
python /path/to/executable/btyper [options...]
```

Note: In the examples below, BTyper commands are shown as ```btyper [options...]```. If you are calling BTyper from the source file (i.e. you didn't install BTyper using Homebrew), keep in mind that you may have to call python and supply the path to btyper to execute the program: ```python btyper [options...]```.


------------------------------------------------------------------------
  
  
  ## Usage and Options
  ### Input File Formats
  
  The command line version of BTyper supports the following file formats as input:
  
1. **Nucleotide sequences in fasta or multifasta format** (1 fasta file per genome or multiple genomes per fasta file)

2. **Draft genome contigs and scaffolds in fasta format** (1 fasta file per genome)

3. **Paired-end ILLUMINA reads in gzipped fastq (.fastq.gz) format** (two separate fastq.gz files per genome, one for forward reads and one for reverse reads)

4. **Single-end ILLUMINA reads in gzipped fastq (.fastq.gz) format** (one fastq.gz file per genome)

5. **Paired-end or single-end ILLUMINA reads in sra format** (one sra file per genome)

6. **Sequence Read Archive (SRA) accession numbers corresponding to genomes sequenced with paired-end or single-end ILLUMINA reads** (one SRA accession number per genome)

### Required Arguments

BTyper can be run from your terminal with the following command line:
  
```
btyper -t [input data type] -i [input file(s)] -o [output directory] [options...]
```

Required arguments are:
  
  
**-t/-\-type [seq, pe, se, sra, sra-get]**
Input data type. Specify the format of the input data using one of the following strings: seq (a file in fasta or multifasta format), pe (ILLUMINA paired-end reads, forward and reverse, in two separate fastq.gz files), se (ILLUMINA single-end reads in one fastq.gz file), sra (Single- or paried-end ILLUMINA reads in one SRA file), sra-get (SRA accession number associated with a genome sequenced using single- or paried-end ILLUMINA reads)

**-i/-\-input [string]**
Path to input fasta file, fastq.gz file(s), sra file, or SRA accession number. For paired-end ILLUMINA reads in separate files, the path to the fastq.gz file containing forward reads should be specified first, followed by the path to the fastq.gz file containing reverse reads, with the two paths separated by a single space.

**-o/-\-output [string]**
Path to desired output directory. Specify the path to the output directory where a results directory (btyper_final_results) containing output files will be created.

### Optional Arguments

Options that can be specified in BTyper include the following:
  
**-\-draft_genome**
For use with a draft genome (contigs or scaffolds) in fasta format (-\-input my_contigs.fasta -\-type seq). If -\-draft_genome is included in the command, BTyper creates a pseudochromosome by concatenating contigs or scaffolds in a fasta file and inserting a spacer sequence ("NNnnNNnnNNnnNNnn") between them. This option is ommitted by default.

**-v/-\-virulence [True or False]**
Virulence gene typing. Performs typing based on presence of virulence genes using a database of genes specific to *Bacillus cereus* group species. Reports genes present at greater than specified percent identity/coverage thresholds. Default is set to True.

**-v_db/-\-virulence_database [aa or nuc]**
Virulence gene database to use for virulence gene detection. Optional argument for use with virulence typing (-\-virulence True). Specify -\-virulence_database nuc to use a nucleotide sequence database and nucleotide blast (blastn), or -\-virulence_database aa to use an amino acid sequence database and translated nucleotide blast (tblastn). Default is set to aa (amino acid sequence database).
Note: tblastn compares a protein query sequence against a nucleotide sequence database dynamically translated in all six reading frames. As a result, -\-virulence_database aa will take longer than -\-virulence_database nuc (virulence gene detection for 70 *B. cereus* group closed assemblies with -\-virulence_database aa and default thresholds takes ~3.5 minutes, while -\-virulence_database nuc and default thresholds takes ~1 minute).

**-aa_p/-\-amino_acid_p [integer between 0 and 100]**
Minimum percent identity for amino acid database. Optional argument for use with virulence typing using an amino acid database (-\-virulence True -\-virulence_database aa). Specify the minimum percent amino acid identity needed for a virulence gene to be considered present in a sequence. Default is set to 50.

**-aa_q/-\-amino_acid_q [integer between 0 and 100]**
Minimum query coverage for amino acid database. Optional argument for use with virulence typing using an amino acid database (-\-virulence True -\-virulence_database aa). Specify the minimum percent coverage needed for a virulence gene to be considered present in a sequence when using an amino acid database. Default is set to 70.

**-nuc_p/-\-nucleotide_p [integer between 0 and 100]**
Minimum percent identity for nucleotide database. Optional argument for use with virulence typing using a nucleotide database (-\-virulence True -\-virulence_database nuc). Specify the minimum percent nucleotide identity needed for a virulence gene to be considered present in a sequence. Default is set to 75.

**-nuc_q/-\-nucleotide_q [integer between 0 and 100]**
Minimum query coverage for nucleotide database. Optional argument for use with virulence typing using a nucleotide database (-\-virulence True -\-virulence_database nuc). Specify the minimum percent coverage needed for a virulence gene to be considered present in a sequence when using a nucleotide database. Default is set to 90.

**-m/-\-mlst [True or False]**
Multilocus sequence typing (MLST). Performs MLST using nucleotide blast (blastn) and the *Bacillus cereus* MLST scheme from PubMLST. Reports highest-scoring alleles using BLAST bit score and their associated sequence type, if available. Default is set to True.

**-r/-\-rpoB [True or False]**
*rpoB* allelic typing. Performs *rpoB* allelic typing using nucleotide blast (blastn) and the *rpoB* allelic typing scheme from Cornell's Food Safety Lab. Reports highest-scoring allele using BLAST bit score. Default is set to True.

**-p/-\-panC [True or False]**
*panC* clade typing. Performs *panC* clade typing using nucleotide blast (blastn) and BTyper's *panC* sequence database. Reports highest-scoring *panC* clade using BLAST bit score. Default is set to True.

**-s [True or False]**
16S rDNA typing. Performs 16S rDNA typing using nucleotide blast (blastn) and 16S genes from nine *B. cereus* group type strains. This method is NOT recommended for typing *B. cereus* group species, as 16S sequences cannot differentiate members of this group at the species level. However, it is included as an option, as some users may be interested in this locus. Reports highest-scoring 16S rDNA gene using BLAST bit score. Default is set to False.

**-e/-\-evalue [float greater than or equal to 0]**
Maximum blast e-value. Optional argument for use with any typing scheme(s). Specify the maximum e-value needed for a BTyper run. Default is set to 1e-5.
Note: this threshold is applied to an entire BTyper run specified by a single command; for example, if you are performing both virulence typing and MLST (-\-virulence True -\-mlst True), and you specify an e-value threshold of 1 (-\-evalue 1), BTyper will apply this maximum threshold to both virulence genes and MLST alleles.

**-\-spades_m [integer]**
Memory limit for SPAdes in Gb. Optional argument for use with ILLUMINA reads (-\-type pe, -\-type se, -\-type sra, or -\-type sra-get). BTyper passes this parameter to the -m/-\-memory option in SPAdes. Default is set to 250, the default for SPAdes.

**-\-spades_t [integer]**
Number of threads for SPAdes. Optional argument for use with ILLUMINA reads (-\-type pe, -\-type se, -\-type sra, or -\-type sra-get). BTyper passes this parameter to the -t/-\-threads option in SPAdes. Default is set to 16, the default for SPAdes.

**-\-spades_k [integer,integer,integer,...]**
Comma-separated list of k-mer sizes to be used for SPAdes (all values must be odd, less than 128 and listed in ascending order). Optional argument for use with ILLUMINA reads (-\-type pe, -\-type se, -\-type sra, or -\-type sra-get). BTyper passes this parameter to the -k option in SPAdes. Default is set to 77.
Note: We recommend selecting optimum k-mer size(s) for your specific data set by consulting the SPAdes documentation. Currently, SPAdes recommends using -k 21,33,55,77 for 150 bp ILLUMINA paired-end reads, and -k 21,33,55,77,99,127 for 250 bp ILLUMINA paired-end reads. 


------------------------------------------------------------------------
  
  
## Output Directories and Files
  
A single BTyper run will deposit the following in your specified output directory (-\-output):
  
**btyper_final_results**
*directory*
Final results directory in which BTyper deposits all of its output files. BTyper creates this directory in your specified output directory (-\-output) 

***your_genome_final_results.txt***
*file*
Final results text file, 1 per input genome. BTyper creates this final results text file, which contains the following, depending on which typing methods you have selected to perform:
  
* **If virulence typing is being performed (-\-virulence True):**
A tab-separated list of virulence genes detected in the genome with the respective e-value, percent identity, and percent coverage for each gene. If a gene is detected multiple times in a genome, BTyper reports only the highest-scoring hit based on its BLAST bit score.

* **If *panC* clade typing is being performed (-\-panC True):**
A tab-separated line, containing the closest-matching *panC* clade (clade1, clade2, ... clade7), the closest-matching *B. cereus* group genome, percent identity, and percent coverage. A *panC* gene that does not match any gene in the database at \geq 75\% identity gives a clade designation of "None" (your isolate may not be a member of the *B. cereus* group), while a *panC* gene that is present at \geq 75\% identity but \leq 90\% identity gives a clade designation of "?" (a *panC* clade could not be determined for your isolate).

* **If MLST is being performed (-\-mlst True):**
A tab-separated line, containing the isolate's (i) sequence type (ST), (ii) *glp* allelic type (AT), (iii) *gmk* AT, (iv) *ilv* AT, (v) *pta* AT, (vi) *pur* AT, (vii) *pyc* AT, and (viii) *tpi* AT. The best-matching allele is reported at each locus; an allele that does not match with 100\% identity or coverage is denoted by an asterisk (\*), while an allele that is not detected in the genome at the given e-value threshold is denoted by "?". If a sequence type cannot be determined using the 7 best-matching allelic types, a "?" is listed in its place. A ST that is detemined using any best-matching alleles that did not match with 100\% identity or coverage is denoted by \*, regardless of whether all 7 alleles could be associated with a ST or not.

* **If *rpoB* allelic typing is being performed (-\-rpoB True):**
A tab-separated line, containing information about the best-matching *rpoB* allele, percent identity, and percent coverage. The *rpoB* allelic information is presented in the following format:

rpoB|ATXXXX|FSL-ID|Genus|species

where *rpoB* refers to the gene,
"ATXXXX"" refers to the allelic type,
"FSL-ID"" refers to Cornell's Food Safety Lab identification number,
Currently, "Genus" refers to the genus of the best-matching organism using NCBI's BLAST server (https://blast.ncbi.nlm.nih.gov/Blast.cgi),
and "species" refers to the species of the best-matching organism using NCBI's BLAST server (https://blast.ncbi.nlm.nih.gov/Blast.cgi)

* **If 16S rDNA typing is being performed (-\-s True):**
A tab-separated line containing strain information of the best-matching of 9 *B. cereus* group type strains, percent identity, and percent coverage. Interpret this information at your own risk, as 16S rDNA sequencing is NOT recommended for typing *B. cereus* group isolates.

**genefiles**
*directory*
Directory in which BTyper deposits genefiles, (multi)fasta files which contain the sequences of all genes detected in a run. BTyper creates this directory within the btyper_final_results directory within your specified output directory (output_directory/btyper_final_results/genefiles).

***some_gene_genefile.fasta***
*file*
BTyper genefiles, (multi)fasta files which contain the sequences of all genes detected in a run. For virulence gene typing, a file will be created for each virulence gene detected in a genome that meet your specified thresholds. The sequence of the best-matching gene/allele using its BLAST bit score is printed to the genefile. If -\-virulence_database aa is selected (the default for BTyper), the amino acid sequence is printed. If -\-virulence_database nuc, is selected, the nucleotide sequence is reported. For MLST, 7 files are created (one for each of the 7 alleles). *rpoB*, *panC*, and 16S typing each produce one genefile. If BTyper is run using more than 1 genome as input (either in multifasta format, or if BTyper is run in a loop), genes from each genome are aggregated together in each genefile. These files are formatted so you can easily input them into your favorite aligner, phylogenetic tree construction program, the NCBI BLAST server, etc.

**isolatefiles**
*directory*
Directory in which BTyper deposits results directories for individual genomes. BTyper creates this directory within the btyper_final_results directory within your specified output directory (output_directory/btyper_final_results/isolatefiles).

***your_genome_results***
*directory*
Directory in which BTyper deposits additional results files for each input genome. BTyper creates this directory within the isolatefiles directory (output_directory/btyper_final_results/isolatefiles/*your_genome*_results). Within this directory, you'll find detailed tab-separated results files for each typing analysis performed, as well as fasta files containing genes extracted from the genome in question. If you're interested in virulence genes present in multiple copies in a genome, the location of each gene in a genome, sequences of all virulence genes detected in a particular isolate, alleles other than the best-matching one, etc., they will be deposited here.


------------------------------------------------------------------------
  
  
## Frequently Asked Questions
  
* **Can I use partial nucleotide sequences (plasmid sequences, MLST genes, *rpoB* alleles, etc.) as input for BTyper?**
  
Sure! You don't have to use whole-genome sequences as input; you can technically use any nucleotide sequencing data and treat it as an assembly, as long as it is in fasta or multifasta format (any of the options that require assembly with SPAdes are designed for bacterial genomes). Although it's not necessary, you may want to adjust the options to only perform typing methods you're interested in to make your output easier to read (i.e. if you're interested in detecting virulence genes in a plasmid sequence, don't waste your time performing other typing methods; just set -\-mlst, -\-rpoB, and -\-panC to False).

* **Can I use whole-genome sequencing data from organisms that don't belong to the *Bacillus cereus* group?**
  
Technically yes! You can use whole-genome sequencing data from any bacterial species as input; in fact, BTyper's *rpoB* allelic type database provided by Cornell's Food Safety Lab actually contains allelic types for non-*Bacillus* species. However, your results from certain typing analyses (i.e. MLST using the *B. cereus* group typing scheme) may be completely meaningless. Be extra cautious when interpreting them, and always take identity and coverage values for detected genes into consideration.


------------------------------------------------------------------------
  
  
## BTyper Tutorial #1: Characterizing a *B. cereus* isolate using its draft genome
  
1. First, let's download our isolate's draft genome from NCBI by clicking the follwoing link:
  
https://www.ncbi.nlm.nih.gov/Traces/wgs/JHQN01

If not already selected, click the "Download" tab. Click on the link for the FASTA file to download the contigs in fasta format: JHQN01.1.fsa_nt.gz. This should download the file into your "Downloads" directory.

2. If you haven't done so already, open your terminal. For Mac users, type **command-space** to open your search bar, type **terminal** in your search bar, and press **Enter**. For Ubuntu users, type **Ctrl-Alt-t** (assuming you haven't changed your default shortcuts).

3. From the command line in your terminal, move to your "Downloads directory" by typing the following, and then hitting **Enter**:
  
```
cd ~/Downloads
```

4. Now that we're in our Downloads directory, let's unzip our contigs file by typing the following command, and hitting **Enter**:
  
```
gunzip JHQN01.1.fsa_nt.gz
```

5. Let's create an output directory in which we can store our BTyper results. That way, once we're done with this tutorial, we can easily delete everything. To create a directory called "btyper_tutorial_1" in our Downloads directory, type the following command and hit **Enter**:
  
```
mkdir ~/Downloads/btyper_tutorial_1
```

6. Now, let's run BTyper on our contigs, directly from our Downloads directory. Because this genome is made up of multiple contigs in multifasta format, rather than a single chromosome, we want to make sure we include the -\-draft_genome option in our command. Because we don't know much about this genome, and we want as much information as possible, let's perform all of the default typing methods with their default settings (virulence gene detection using an amino acid database, MLST, *rpoB* allelic typing, and *panC* clade typing, all using default thresholds), as well as 16S gene detection, just for fun. To run BTyper, type the following command, and press **Enter**:

```
btyper -t seq -i ~/Downloads/JHQN01.1.fsa_nt -o ~/Downloads/btyper_tutorial_1 -s True --draft_genome
```

Here are the options we selected, explained:

* **-t seq**
Because our sequence is in fasta format, we're declaring our input type to be seq

* **-i ~/Downloads/JHQN01.1.fsa_nt**
We're directing BTyper to our input file (JHQN01.1.fsa_nt)

* **-o ~/Downloads/btyper_tutorial_1**
We're telling BTyper where to store the output files it produces (our directory, btyper_tutorial_1)

* **-s True**
Because 16S typing is not performed by default, we're setting this option to True

* **-\-draft_genome**
We're working with a single genome formed by multiple contigs in a single file. We're going to tell BTyper to concatenate these contigs into a pseudochromosome, rather than treat each as a separate genome.

7. Once the program is finished running (this should take about 11 seconds or so), we can take a look at our final results file to get detailed results about our isolate. Open the "JHQN01_final_results.txt" file in any text editor, either by searching for it or opening it in your "~/Downloads/btyper_tutorial_1/btyper_final_results" directory.

8. Now that we've opened the final results file for our isolate, we can see information about our isolate, divided into several parts:
  
* **Predicted Virulence Proteins**
  
This is a list of proteins detected in our genome at 70\% coverage and 50\% identity using an amino acid sequence database (BTyper's default settings for virulence typing). Looking through the list of virulence genes, it seems that our isolate possesses a couple of *B. anthracis*-associated genes with really high similarity, including anthrax toxin genes *cya, lef,* and *pagA*!

* **Predicted *panC* Clade Designation**
This corresponds to our isolate's *panC* clade assignment. It looks like the closest-matching *panC* sequence was that of *B. cereus* 03BB87, which belongs to *panC* Clade III, the same clade as *B. anthracis*.

* **Predicted MLST Profile**

This section contains our isolate's allelic types for 7-gene MLST, as well as the associated sequence type. It looks like our *B. cereus* isolate belongs to ST 78, matching with 100\% identity.

* **Predicted *rpoB* Allelic Type**

This section contains the closest-matching *rpoB* allelic type for our isolate. Our isolate's *rpoB* allele appears to match allelic type 365 exactly, which matched *B. cereus* most closely when the NCBI BLAST server was used.

* **Predicted 16S Type**

This section contains the closest-matching 16S gene of *B. cereus* group type strains. While out isolate looks like it matches *B. cereus* ATCC 14579 most closely, we should interpret all 16S typing results with caution.

Sure enough, this genome is actually that of *B. cereus* strain BcFL2013, which was isolated from a patient with an anthrax-like skin lesion in Florida (Gee, et al., 2014, *Genome Announcements*). It's been shown to possess anthrax toxin genes and belong to ST 78, the same ST as several other anthrax-causing *B. cereus* strains (Gee, et al., 2014, *Genome Announcements*).

9. If you want to delete the results from this tutorial, just go to your Downloads folder and delete the "btyper_tutorial_1" directory there.

------------------------------------------------------------------------

## BTyper Tutorial #2: Extracting *plcR* nucleotide sequences from 3 *B. cereus* group assemblies using their closed chromosomes

1. Open your terminal, as described in Tutorial #1

2. First, let's move to our Downloads directory by typing the following, and pressing **Enter**:

```
cd ~/Downloads
```

3. Let's create a project directory; type the following command to create a directory called "btyper_tutorial_2":

```
mkdir btyper_tutorial_2
```

4. Let's download 3 *B. cereus* group chromosomes from NCBI, starting with *B. anthracis* Ames Ancestor:

* Click on the following link to go to the genome page: https://www.ncbi.nlm.nih.gov/nuccore/NC_007530.2

* In the upper-right corner, click "Send"->"Complete Record"->"File"->"FASTA". The sequence should begin downloading as a file named "sequence.fasta".

* Move your sequence to the btyper_tutorial_2 directory and rename it to ames_ancestor.fasta by typing the following command and pressing **Enter:**

```
mv sequence.fasta btyper_tutorial_2/ames_ancestor.fasta
```
                                                                                                              
* Repeat steps (i) through (iii), replacing the URL and new file name (ames_ancestor.fasta) with the following:
                                                                                                                               https://www.ncbi.nlm.nih.gov/nuccore/NC_004722.1 (rename to atcc_14579.fasta), 
                                                                                                                             https://www.ncbi.nlm.nih.gov/nuccore/NC_005957.1 (rename to konkukian.fasta)
                                                                                                                              
                                                                                                                              
5. Move to your btyper_tutorial_2 by typing the following command and pressing **Enter:**

```
cd ~/Downloads/btyper_tutorial_2
```

                                                                                                                             
6. Next, let's concatenate all 3 chromosomes together to form a multifasta called "bacillus.fasta" by typing the following command and pressing **Enter:**

```
cat *.fasta > bacillus.fasta
```
                                                                                                                             
7. Run BTyper to extract the *plcR* nucleotide sequence from each of our 3 genomes by typing the following command: 

```
btyper -t seq -i bacillus.fasta -o . -m False -r False -p False -v_db nuc -nuc_p 0 -nuc_q 0
```

Here are the options we selected, explained:

* **-t seq**

Because our sequence is in multifasta format, we're declaring our input type to be seq

* **-i bacillus.fasta**

We're directing BTyper to our input file

* **-o .** 

We're telling BTyper where to store the output files in our current directory (btyper_tutorial_2)

* **-m False**

We only want to perform virulence typing, so we are setting MLST to False

* **-r False**

We only want to perform virulence typing, so we are setting *rpoB* allelic typing to False

* **-p False**

We only want to perform virulence typing, so we are setting *panC* clade typing to False

* **-v_db nuc**

We want to extract virulence gene nucleotide sequences, so we are telling BTyper to use the nucleotide sequence database

* **-nuc_p 0**

We're lowering the percent identity threshold for virulence gene detection using a nucleotide database

* **-nuc_q 0**

We're lowering the percent coverage threshold for virulence gene detection using a nucleotide database

8. We can find our virulence gene nucleotide sequences deposited in the genefiles folder in our btyper_final_results directory. Our *plcR* nucleotide sequences will be labeled as plcR_genefile.fasta. We can then align the sequences in this file, build a phylogenetic tree, call SNPs, etc.

9. If you want to delete the results from this tutorial, just go to your Downloads folder and delete the "btyper_tutorial_2" directory there.

------------------------------------------------------------------------

## BTyper Tutorial #3: Assembling and characterizing the genome of a clinical *B. cereus* isolate using its SRA accession number

1. First, open up your terminal, as described in Tutorial #1

2. Let's create a new directory named btyper_tutorial_3 in our Downloads directory by typing the following command into our terminal and pressing **Enter:**

```
mkdir ~/Downloads/btyper_tutorial_3
```

3. We're going to be assembling the genome and characterizing the following isolate, so let's make sure it meets BTyper's criteria for assembly: https://www.ncbi.nlm.nih.gov/sra/ERX1840887

We need to make sure this genome was sequenced using ILLUMINA reads (either single- or paired-end reads are fine; BTyper can infer this); according to SRA, it was sequenced using NextSeq 500, an ILLUMINA platform, so we should be good to go! Our isolate's SRA accession number appears to be ERR1775894.

4. It looks like this genome was sequenced with 150 bp paired-end reads. If we were in a hurry, the default k-mer size parameters that BTyper passes to SPAdes to assemble the genome (-\-spades_k 77) should suffice. However, SPAdes is a great assembler, and we want to take advantage of that; let's try to produce an optimal assembly for 150 bp reads by testing multiple k-mer sizes and having SPAdes pick the best one. To do this, type the following command and press **Enter.** This will download sequence data from SRA, assemble and correct mismatches in the genome using SPAdes, and perform all default typing methods using BTyper...but be prepared to wait a little while (about 30 minutes to run the example below, but this depends on your computer, the k-mer sizes you select, your memory/thread parameters, etc.)

```
btyper -t sra-get -i ERR1775894 -o ~/Downloads/btyper_tutorial_3 --spades_k 21,33,55,77 --spades_m 8 --spades_t 8
```

Here is our command, explained:

* **-t sra-get**

Tell BTyper that our sequencing data type is an SRA accession number associated with ILLUMINA reads

* **-i ERR1775894**

This is our SRA accession number that we'd like BTyper to search for.

* **-o ~/Downloads/btyper_tutorial_3**

This is the path to the output directory in which we want BTyper to store our results (including our assembled genomes)

* **-\-spades_k 21,33,55,77**

This is a list of k-mer sizes we want SPAdes to try, recommended by the developers for use with 150 bp paired-end ILLUMINA reads (http://spades.bioinf.spbau.ru/release3.5.0/manual.html#sec3.4)

* **-\-spades_m 8**

This is the memory limit we want to pass to SPAdes, in Gb. Here, we're setting ours to 8 Gb, but feel free to change this, depending on your machine!

* **-\-spades_t 8**

This is the number of threads we want SPAdes to use. Here, we're using 8 threads, but feel free to change this, depending on your machine!

5. Once our command is finished running, we can head to our output directory (btyper_tutorial_3). In addition to the usual btyper_final_results directory that BTyper produces, we should also see some additional files:

* **ERR1775894_1.fastq.gz**, our forward reads from SRA

* **ERR1775894_2.fastq.gz**, our reverse reads from SRA; if we were to input an SRA accession number associated with ILLUMINA single-end reads, we would only have 1 fastq.gz file of reads

* **ERR1775894_1_spades_assembly.fasta**, our file of contigs produced using SPAdes. If more than one k-mer size is tested (like we did here), BTyper chooses the optimal one, as determined by SPAdes. In this case, a k-mer size of 77 (BTyper's default option for -\-spades_k) produced the best assembly, so BTyper selected this assembly for its analysis.

* **ERR1775894_1_pseudochrom.fasta**, our pseudochromosome, formed by concatenating the contigs of our optimal assembly (in this case, ERR1775894_1_spades_assembly.fasta) using BTyper's -\-draft_genome option.

* **spades_assembly**, a directory produced by SPAdes (http://spades.bioinf.spbau.ru/release3.5.0/manual.html), containing its output files                                                                                                                           
6. We can then go into our btyper_final_results directory to look at our typing results. From a quick glance at the ERR1775894_1_final_results.txt file, it looks like our isolate belongs to MLST sequence type 26, *rpoB* allelic type 125, and *panC* clade 3. Not only that, it looks like our isolate may produce the *B. cereus* emetic toxin; *cesABCD* were all detected in the assembly at high identity/coverage!

7. To delete results from this tutorial, just delete the btyper_tutorial_3 folder in your Downloads folder.


------------------------------------------------------------------------


## References

#### Dependencies

Bankevich, Anton, et al. SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing. *Journal of Computational Biology* 2012 May; 19(5): 455-477.

Camacho, Christiam, et al. BLAST+: architecture and applications. *BMC Bioinformatics* 2009 10:421.

Cock, Peter J., et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 2009 June 1; 25(11): 1422-1423.

Leinonen, Rasko, et al. The Sequence Read Archive. *Nucleic Acids Research* 2011 Jan; 39(Database issue): D19–D21.

#### Typing Methods

Guinebretière, MH, et al. Ecological diversification in the *Bacillus cereus* Group. *Environmental Microbiology* 2008 April; 10(4): 851-65.

Guinebretière, MH, et al. Ability of *Bacillus cereus* Group Strains To Cause Food Poisoning Varies According to Phylogenetic Affiliation (Groups I to VII) Rather than Species Affiliation. *Journal of Clinical Microbiology* 2010 September; 8(9):3388-91.

Ivy, RA, et al. Identification and characterization of psychrotolerant sporeformers associated with fluid milk production and processing. *Applied and Environmental Microbiology* 2012 March; 78(6):1853-64.

Kovac, Jasna, et al. Production of hemolysin BL by *Bacillus cereus* group isolates of dairy origin is associated with whole-genome phylogenetic clade. *BMC Genomics* 2016 17:581.

PubMLST *Bacillus cereus* MLST database (https://pubmlst.org/bcereus/), based on Jolley, Keith A. and Martin CJ Maiden. BIGSdb: Scalable analysis of bacterial genome variation at the population level. *BMC Bioinformatics* 2010 11:595.

Rossi-Tamisier, M., et al. Cautionary tale of using 16S rRNA gene sequence similarity values in identification of human-associated bacterial species. *International Journal of Systematic and Evolutionary Microbiology* (2015), 65, 1929–1934.

#### Tutorial Genomes

Gee, JE, et al. Draft Genome Sequence of *Bacillus cereus* Strain BcFL2013, a Clinical Isolate Similar to G9241. *Genome Announcements* 2014 May 29;2(3).


------------------------------------------------------------------------


Disclaimer: BTyper and BMiner are pretty neat! However, no tool is perfect, and BTyper and BMiner cannot definitively prove whether a *B. cereus* group isolate is pathogenic or not. As always, interpret your results with caution. We are not responsible for taxonomic misclassifications, misclassifications of an isolate's pathogenic potential, and/or misinterpretations (biological, statistical, or otherwise) of BTyper and/or BMiner results.











