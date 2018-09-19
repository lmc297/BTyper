# BTyper

A computational tool for virulence-based classification of *Bacillus cereus* group isolates and/or antimicrobial resistance gene detection using nucleotide sequencing data

## Overview

BTyper is a command-line tool that employs a combination of *in silico* (i) virulence gene detection, (ii) multi-locus sequence typing (MLST), (iii) *panC* clade assignment, and (iv) *rpoB* allelic typing to rapidly classify *B. cereus* group isolates using nucleotide sequencing data.

Antimicrobial resistance (AMR) gene detection was added in BTyper version 2.0.0 (released 2017-06-29), a function that can be used with sequencing data from any bacterial species.

Most recently, an average nucleotide identity blast (ANIb) option was added in BTyper version 2.3.0 (released 2018-09-18)

The program, as well as the associated databases, can be downloaded from https://github.com/lmc297/BTyper.

Post issues at https://github.com/lmc297/BTyper/issues

BTyper v. 1.0.0 output files for the 662 genomes used in the manuscript can be downloaded from https://github.com/lmc297/BTyper/raw/master/archive/sample_data/ncbi_btyper_final_results.zip.


### Citation

#### If you found the BTyper tool, its source code, and/or any of its associated databases useful, please cite:
  
Carroll, Laura M., Jasna Kovac, Rachel A. Miller, Martin Wiedmann. 2017. Rapid, high-throughput identification of anthrax-causing and emetic *Bacillus cereus* group genome assemblies using BTyper, a computational tool for virulence-based classification of *Bacillus cereus* group isolates using nucleotide sequencing data. *Applied and Environmental Microbiology* 2017 Sep 1; 83(17): e01096-17.

#### Additionally, if you found BTyper's antimicrobial resistance (AMR) gene detection methods useful (particularly if you used them in non-*Bacillus cereus* group species), you may find the following paper to be helpful; in it, we implement and validate the AMR gene detection method employed by BTyper in *Salmonella enterica*:

Carroll, Laura M., Martin Wiedmann, Henk den Bakker, Julie Siler, Steven Warchocki, David Kent, Svetlana Lyalina, Margaret Davis, William Sischo, Thomas Besser, Lorin D. Warnick, Richard V. Pereira. Whole-Genome Sequencing of Drug-Resistant *Salmonella enterica* Isolated from Dairy Cattle and Humans in New York and Washington States Reveals Source and Geographic Associations. *Applied and Environmental Microbiology* 2017 May 31;83(12).


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

**Aggregate multiple BTyper final results files into one matrix/text file:

```
btyper2matrix.py -i </path/to/directory/btyper_final_results/> -o </path/to/desired/output/directory/>
```

------------------------------------------------------------------------
  
  
  ## Installation
  ### Install BTyper using Homebrew (macOS users)
  
  BTyper and its dependencies can be installed using <a href="https://brew.sh/">Homebrew</a>.

1. First, install Homebrew, if necessary, by running the following command from your terminal:
  
```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

2. Install pip, if necessary, by running the following command from your terminal:

```
sudo easy_install pip
```

3. Install Biopython, if necessary, by running the following command from your terminal:

```
pip install biopython
```
Note: if you don't have permissions, you may need to use sudo:

```
sudo pip install biopython
```

4. Tap brewsci/science, if necessary, by running the following command from your terminal:

```
brew tap brewsci/science
```

5. Tap BTyper by running the following command from your terminal:
  
```
brew tap lmc297/homebrew-btyper
```

6. Install BTyper and its dependencies by running the following command from your terminal:
  
```
brew install btyper
```

7. Optional: in BTyper version 2.3.0 (released 18-09-18), ANIb has been added as an optional typing method (see below for details). If you want to use BTyper for ANIb, download the "published" or "effective" *B. cereus* group species ANIb database(s) by running one of the following commands from your terminal:

For published database (recommended; needs about 118M disk space):

```
build_btyper_anib_db.py -db published
```

For effective database (not recommended; needs about 237M disk space):
```
build_btyper_anib_db.py -db effective
```

After running either command, follow the instructions in your terminal.


### Download and run BTyper using source file (macOS and Ubuntu)

1. To run BTyper, please download and install the following dependencies, if necessary:
  
  <a href="https://www.python.org/downloads/"> Python v. 2.7</a>
  
  <a href="http://biopython.org/DIST/docs/install/Installation.html"> Biopython v. 1.6.9</a>
  
  <a href="https://github.com/Homebrew/homebrew-science/blob/master/blast.rb">BLAST+ v. 2.4.0 or higher</a>
  
  <a href="https://github.com/Homebrew/homebrew-science/blob/master/spades.rb">SPAdes v. 3.9.0 or higher</a>
  
  <a href="https://github.com/ncbi/sra-tools/wiki/HowTo:-Binary-Installation">SRA toolkit v. 2.8.0 or higher</a>
  
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

Note: In the examples below, BTyper commands are shown as ```btyper [options...]```. If you are calling BTyper from the source file (i.e. you didn't install BTyper using Homebrew), keep in mind that you may have to call python and supply the path to btyper to execute the program or related scripts: ```python btyper [options...]```.

5. Optional: in BTyper version 2.3.0 (released 18-09-18), ANIb has been added as an optional typing method (see below for details). If you want to use BTyper for ANIb, download the "published" or "effective" *B. cereus* group species ANIb database(s) by running one of the following commands from your terminal:

For published database (recommended; needs about 118M disk space):

```
python /path/to/executable/build_btyper_anib_db.py -db published
```

For effective database (not recommended; needs about 237M disk space):
```
python /path/to/executable/build_btyper_anib_db.py -db effective
```

After running either command, follow the instructions in your terminal.

Note: In the examples below, BTyper commands are shown as ```btyper [options...]```. If you are calling BTyper from the source file (i.e. you didn't install BTyper using Homebrew), keep in mind that you may have to call python and supply the path to btyper to execute the program or related scripts: ```python btyper [options...]```.


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

**-b/-\-anib [True or False]**
Average nucleotide identity blast (ANIb). Performs ANIb using the algorithm described in the "Average Nucleotide Identity BLAST (ANIb) Q&A" section below, treating the input genome as the query and selected *B. cereus* group genomes as references (see -b_db/-\-anib_db option below). Reports reference genome with the highest ANIb value. Default is set to False.
Note: ANIb is a good method for predicting bacterial species. However, running ANIb (-\-anib True) will significantly increase your analysis time. For example: a single *B. cereus* group draft genome (-\-draft_genome) takes under 10 seconds to analyze using BTyper's default settings; using -\-anib True (with the default published genome database), it takes about 1 minute.

**-b_db/-\-anib_db [published or effective]**
Reference genome database to use for ANIb. Optional argument for use with ANIb (-\-anib True). Specify reference genomes for ANIb calculation: published for 18 published *Bacillus cereus* group species, or effective to use the published database plus 21 effective *Bacillus* cereus group species that have been proposed in the literature but not published as a novel species (39 species total). Default is set to published.
Note: Running ANIb with the effective species genome database (-\-anib_db effective) will significantly increase your analysis time. For example: a single *B. cereus* group draft genome (-\-draft_genome) takes about 1 minute to run using the published (default) genome database; using -\-anib_db effective, it takes about 2 minutes.

**-m/-\-mlst [True or False]**
Multilocus sequence typing (MLST). Performs MLST using nucleotide blast (blastn) and the *Bacillus cereus* MLST scheme from PubMLST. Reports highest-scoring alleles using BLAST bit score and their associated sequence type, if available. Default is set to True.

**-r/-\-rpoB [True or False]**
*rpoB* allelic typing. Performs *rpoB* allelic typing using nucleotide blast (blastn) and the *rpoB* allelic typing scheme from Cornell's Food Safety Lab. Reports highest-scoring allele using BLAST bit score. Default is set to True.

**-p/-\-panC [True or False]**
*panC* clade typing. Performs *panC* clade typing using nucleotide blast (blastn) and a BTyper *panC* sequence database (see -panC_db/-\-panC_database option below). Reports highest-scoring *panC* clade using BLAST bit score. Default is set to True.

**-panC_db/-\-panC_database [type, latest, or legacy]**
*panC* gene database to use for *panC* clade typing. Optional argument for use with *panC* clade typing (-\-panC True). In version 2.2.2, BTyper uses *panC* sequences from 9 type strains (*B. anthracis, B. cereus sensu stricto, B. cytotoxicus, B. mycoides, B. pseudomycoides, B. thuringiensis, B. toyonensis, B. weihenstephanensis, B. wiedmannii*) to perform 7-clade *panC* clade assignment by default. **It is highly recommended you use this method, as it shows high correlation with whole-genome phylogenetic clade (I-VII).** Default is set to type (7-clade *panC* database with 9 species type strains).

If you want to use *panC* sequences of 18 *B. cereus* group type strains, specify -\-panC_database latest to use the most recent *panC* database, which contains *panC* sequences from type strains for 18 *B. cereus* group species. This new database includes the 9 original species type strains (*B. anthracis, B. cereus sensu stricto, B. cytotoxicus, B. mycoides, B. pseudomycoides, B. thuringiensis, B. toyonensis, B. weihenstephanensis, B. wiedmannii*), as well as the 9 new species proposed by <a href="https://www.ncbi.nlm.nih.gov/pubmed/28792367"> Liu, et al., 2017</a>. Instead of outputting numerical clades (clade1, clade2, ... clade7), the new database assigns clades by species name to avoid arbitrary numbering of new clades (e.g. cladeAnthracis, cladeCereus, cladeTropicus, etc.). If you are interested in the *panC* gene itself, this may be a valuable option; however, if you are interested in inferring whole-genome phylogenetic clade from *panC*, preliminary results indicate that, for some of the 9 new species, whole-genome phylogenetic clade does not correlate with whole-genome clade. This option is not recommended.

In addition, the original, 7-clade typing scheme using BTyper's original *panC* allele database can be specified using -\-panC_database legacy. The legacy option is ***not*** recommended, as there are some *panC* alleles from non-type strains that, despite giving the same results as the original *panC* clade assignment tool described by <a href="https://www.tools.symprevius.org/bcereus/">Guinebretière, et al. 2010</a>, do not correlate with whole-genome clade. If you want to perform 7-clade *panC* assignment, please use the 9-species type strain (-panC_db type) method that is implemented by default in BTyper version 2.2.2 and up.

**-s [True or False]**
16S rDNA typing. Performs 16S rDNA typing using nucleotide blast (blastn) and 16S genes from 18 *B. cereus* group type strains. This method is NOT recommended for typing *B. cereus* group species, as 16S sequences cannot differentiate members of this group at the species level. However, it is included as an option, as some users may be interested in this locus. Reports highest-scoring 16S rDNA gene using BLAST bit score. Default is set to False.

**-s_db/-\-s_database [latest or legacy]**
16S rDNA gene database to use for 16S rDNA typing. Optional argument for use with 16S rDNA typing (-s True). Specify -\-s_database latest to use the most recent 16S rDNA database, which contains 16S rDNA sequences from type strains for 18 *B. cereus* group species. This new database includes the 9 original species (*B. anthracis, B. cereus sensu stricto, B. cytotoxicus, B. mycoides, B. pseudomycoides, B. thuringiensis, B. toyonensis, B. weihenstephanensis, B. wiedmannii*), as well as the 9 new species proposed by <a href="https://www.ncbi.nlm.nih.gov/pubmed/28792367"> Liu, et al., 2017</a>. The original, 9-species database can be specified using -\-s_database legacy. However, as always, **interpret any 16S rDNA typing results for the *B. cereus* group with extreme caution!** Default is set to latest (18-species 16S rDNA database). 

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

**-a/-\-amr [True or False]**
Antimicrobial resistance gene detection. Detects antimicrobial resistance genes, and reports genes present at greater than specified percent identity/coverage thresholds. Default is set to True. Note: in BTyper version 2.2.0, the -\-amr option can be used to detect plasmid replicons in nucleotide sequences using the PlasmidFinder database; see below for more information.

**-amr_db/-\-amr_database [argannot or megares]**
Antimicrobial resistance (AMR) gene database to use for antimicrobial resistance gene detection. Optional argument for use with antimicrobial resistance gene detection (-\-amr True). Specify -\-amr_database argannot to use the ARG-ANNOT AMR gene nucleotide database or -\-amr_database megares to use the MEGARes AMR gene nucleotide database. Default is set to argannot (ARG-ANNOT database). Note: in BTyper version 2.2.0, plasmid replicon detection can be performed using the PlasmidFinder database by specifying -\-amr_database plasmidfinder; all options that can be used with AMR gene detection (e.g. minimum percent identity and coverage thresholds, pruning method, overlap) described below can be applied to plasmid replicon detection.

**-\-amr_p  [integer between 0 and 100]**
Minimum percent nucleotide identity for antimicrobial resistance gene detection. Optional argument for use with antimicrobial resistance gene detection (-\-amr True). Specify the minimum percent nucleotide identity needed for an antimicrobial resistance gene to be considered present in a sequence. Default is set to 75.

**-\-amr_q  [integer between 0 and 100]**
Minimum nucleotide query coverage for antimicrobial resistance gene detection. Optional argument for use with antimicrobial resistance gene detection (-\-amr True). Specify the minimum percent query coverage needed for an antimicrobial resistance gene to be considered present in a sequence. Default is set to 50.

**-\-amr_blast [blastn or tblastx]**
BLAST algorithm to use for antimicrobial resistance gene detection. Optional argument for use with antimicrobial resistance gene detection (-\-amr True). Specify -\-amr_database blastn to use nucleotide BLAST against the selected antimicrobial resistance gene database, or -\-amr_database tblastx to use a translated nucleotide database and a translated nucleotide query. Default is set to blastn (nucleotide BLAST).
Note: tblastx compares a nucleotide query sequence against a nucleotide sequence database, both of which are dynamically translated in all six reading frames. As a result, -\-amr_blast tblastx will take significantly longer than -\-amr_blast blastn (antimicrobial resistance gene detection using the ARG-ANNOT database for 9 *B. cereus* group draft assemblies (contigs) with -\-amr_blast tblastx and default thresholds takes ~24 minutes, while -\-amr_blast blastn and default thresholds takes ~30 seconds).

**-\-prune [cluster or location]**
Pruning method to use for antimicrobial resistance gene detection. Optional argument for use with antimicrobial resistance gene detection (-\-amr True). Specify pruning method for final results files, which determines which alleles are reported in final results files. Specify -\-prune cluster to report top-scoring gene within a detected cluster of similar genes (obtained by clustering the appropriate database at 0.8 identity using cd-hit-est, with the top-scoring allele reported from each cluster), or -\-prune location to report the top-scoring gene for genes overlapping past some threshold (see -\-overlap below). Default is set to location. Note: BTyper versions 2.0.0 to 2.1.0 used the -\-prune cluster method; however, we have determined the -\-prune location method to be more accurate, hence its implementation as the default pruning method in BTyper version 2.2.0.

**-\-overlap [float between 0 and 1]**
Maximum proportion of overlap for overlapping antimicrobial resistance genes to be considered separate genes, rather than alleles of the same gene. Optional argument for use with antimicrobial resistance gene detection (-\-amr True) and location-based pruning method (-\-prune location). Hits overlapping above this threshold are considered to be alleles of the same gene, and the highest-scoring hit is reported in the final results file. Hits with overlap proportions below this threshold are considered separate genes, and both are reported in the final results file. Default is set to 0.7. 


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

* **If antimicrobial resistance gene detection is being performed (-\-amr True):**
A tab-separated list of antimicrobial resistance genes detected in the genome with the respective e-value, percent identity, and percent coverage for each gene. The highest-scoring allele (using its blast bitscore) from its respective gene cluster/location (depending on pruning method) is reported. Additionally, if a gene is detected multiple times in a genome, BTyper reports only the highest-scoring hit based on its BLAST bit score when -\-prune cluster is used. When -\-prune location is used, alleles of the same gene that appear in different locations in the genome (i.e. multiple copies) will be reported here. Note: in BTyper version 2.2.0, if the PlasmidFinder database is being used to detect plasmid replicons, they will be reported here in lieu of AMR genes.

* **If average nucleotide identity BLAST (ANIb) is being performed (-\-anib True):**
A tab-separated line, containing the *B. cereus* group species with the highest ANIb value (one of the 18 published *B. cereus* group species if the default published ANIb database is used; if the effective ANIb database is used, the closest "species" may be denoted by a RefSeq accession number rather than a species name), the ANIb value, and the percen coverage. If no tested *B. cereus* group species has an ANIb value > 95 (often regarded as a cutoff for bacterial species), a species of "Unknown" is reported, with the highest-ANIb-producing species in parentheses plus \* (e.g., "Unknown (Bacillus pseudomycoides)\*")

* **If *panC* clade typing is being performed (-\-panC True):**
A tab-separated line, containing the closest-matching *panC* clade (clade1, clade2, ... clade7 if the type strain 7-clade typing scheme is used, or cladeAlbus, cladeAnthracis, clade Cereus, ... cladeWiedmannii if the latest 18-species typing scheme is used), the closest-matching *B. cereus* group genome, percent identity, and percent coverage. A *panC* gene that does not match any gene in the database at &#8805; 75\% identity gives a clade designation of "None" (your isolate may not be a member of the *B. cereus* group), while a *panC* gene that is present at &#8805; 75\% identity but &#8804; 90\% identity gives a clade designation of "?" (a *panC* clade could not be determined for your isolate).

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
A tab-separated line containing strain information of the best-matching of 18 *B. cereus* group type strains (9 if using the legacy 16S rDNA database), percent identity, and percent coverage. Interpret this information at your own risk, as 16S rDNA sequencing is NOT recommended for typing *B. cereus* group isolates.

**genefiles**
*directory*
Directory in which BTyper deposits genefiles, (multi)fasta files which contain the sequences of all genes detected in a run. BTyper creates this directory within the btyper_final_results directory within your specified output directory (output_directory/btyper_final_results/genefiles).

***some_gene_genefile.fasta***
*file*
BTyper genefiles, (multi)fasta files which contain the sequences of all genes detected in a run. For virulence gene typing and antimicrobial resistance (AMR) gene detection, a file will be created for each virulence/AMR gene detected in a genome that meet your specified thresholds. The sequence of the best-matching gene/allele using its BLAST bit score is printed to the genefile. For virulence typing, if -\-virulence_database aa is selected (the default for BTyper), the amino acid sequence is printed. If -\-virulence_database nuc, is selected, the nucleotide sequence is reported. For MLST, 7 files are created (one for each of the 7 alleles). *rpoB*, *panC*, and 16S typing each produce one genefile. If BTyper is run using more than 1 genome as input (either in multifasta format, or if BTyper is run in a loop), genes from each genome are aggregated together in each genefile. These files are formatted so you can easily input them into your favorite aligner, phylogenetic tree construction program, the NCBI BLAST server, etc.

**isolatefiles**
*directory*
Directory in which BTyper deposits results directories for individual genomes. BTyper creates this directory within the btyper_final_results directory within your specified output directory (output_directory/btyper_final_results/isolatefiles).

***your_genome_results***
*directory*
Directory in which BTyper deposits additional results files for each input genome. BTyper creates this directory within the isolatefiles directory (output_directory/btyper_final_results/isolatefiles/*your_genome*_results). Within this directory, you'll find detailed tab-separated results files for each typing analysis performed, as well as fasta files containing genes extracted from the genome in question. If you're interested in virulence/AMR genes present in multiple copies in a genome, the location of each gene in a genome, sequences of all virulence genes detected in a particular isolate, alleles other than the best-matching one, etc., they will be deposited here.


If average nucleotide identity BLAST (ANIb) is being performed, an additional output directory with additional files is created:

**anib_blastfiles**
*directory*
Directory in which BTyper deposits directories containing average nucleotide identiy BLAST (ANIb) calculation files for each input genome; most users can delete this directory, as it contains intermediate files used to calculate ANIb relative to all reference genomes. BTyper creates this directory within the btyper_final_results directory within your specified output directory (output_directory/btyper_final_results/anib_blastfiles).

***your_genome***
*directory*
Directory in which BTyper deposits ANIb calculation files for each input genome. BTyper creates this directory within the anib_blastfiles directory (output_directory/btyper_final_results/anib_blastfiles/*your_genome*). Within this directory, you will find raw BLAST results (ends in .txt) and the hits used to calculate ANIb/coverage (ends in finalfragments.txt) for each input genome/reference genome combination; for example, if you have an input genome named mygenome.fasta, you will have a file called mygenome_vs_anthracis.txt (raw BLAST results of mygenome relative to *B. anthracis* str. Ames genome), a file called mygenome_vs_anthracis_finalfragments.txt (values used to calculate ANIb between mygenome and the *B. anthracis* str. Ames genome), mygenome_vs_cereus.txt (raw blast results of mygenome relative to *B. cereus* str. ATCC 14579 genome), a file called mygenome_vs_cereus_finalfragments.txt (values used to calculate ANIb between mygenome and *B. cereus* str. ATCC 14579 genome), and so on...)


------------------------------------------------------------------------

## Additional BTyper Scripts (BTyper version 2.3.0 and up)

### btyper2matrix.py

* Purpose: aggregates multiple BTyper file results files and produce a single matrix/text file (easier to read/interpret for users with more than one genome, easier to parse)

* Input: path to directory containing btyper final results files; these files have the suffix '_final_results.txt' and should be in a directory called 'btyper_final_results', provided the directory was not renamed

* Output: matrix file called "btyper2matrix_output.txt", deposited in an output directory of your choice; matrix is tab-separated and contains 1 row per final results file (genome). Includes a column for each virulence gene detected (denoted in column header by "vir|gene_name", with 1 for present in a genome and 0 for absent in a genome), a column for each AMR gene detected (denoted in column header by "amr|gene_name", with 1 for present in a genome and 0 for absent in a genome), and one column each for ANI species, *panC* clade, MLST sequence type, *rpoB* allelic type, and 16S rDNA species; cells with "NA" correspond to the absence of results for a particular analysis.

* Command structure:
```
btyper2matrix.py -i </path/to/directory/btyper_final_results/> -o </path/to/output/directory/>
```
For help, type ```btyper2matrix.py -h``` or ```btyper2matrix.py --help```

### build_btyper_anib_db.py

* Purpose: download databases(s) to be used with BTyper's ANIb option (-\-anib True); must be run before running ANIb

* Input: published or effective; specify the ANIb database to download for use with BTyper's -b/-\-anib option; published for 118M database with 18 published *Bacillus cereus* group species, or effective for 237M database, which includes published database plus 23 effective *Bacillus cereus* group species that have been proposed in the literature but not published as a novel species (41 species total; for most users' purposes, the published database is recommended)

* Output: 18 genomes stored in BTyper's seq_anib_db/published/ directory (published) or 18 genomes stored in BTyper's seq_anib_db/published/ directory, plus 21 genomes stored in BTyper's seq_anib_db/effective/ directory (effective)

* Command structure:
```
build_btyper_anib_db.py -db [published or effective]
```
Users will then be prompted in the terminal to type "yes" and press ENTER to confirm the download.

For help, type ```build_btyper_anib_db.py -h``` or ```build_btyper_anib_db.py --help```

------------------------------------------------------------------------

## Average Nucleotide Identity BLAST (ANIb) Q&A

### What is ANI?

Average nucleotide identity, or ANI, is a metric that can be used to <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC549018/">assign a bacterial genome to a bacterial species</a>. Very briefly, it is performed by fragmenting the bacterial genome in question (e.g., into fragments of a specific length, into coding regions), aligning those fragments against a reference genome (e.g., using BLAST, MUMmer), and then taking the average percentage of identical base pair matches of the aligned regions. 

In the past, DNA–DNA hybridization (DDH) has been used to delineate prokaryotic species, with memebers of a species <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2776425/">generally sharing DDH values of greater than 70% similarity</a>. ANI is a rapid and scalable *in silico* method that has been shown to correlate with DDH values, so it can be used to evaluate whether two genomes belong to the same species or different species. Currently, a <a href="https://www.ncbi.nlm.nih.gov/pubmed/17220447">species cutoff of 95 ANI is typically used for assigning bacterial species</a> (i.e., 2 genomes belong to the same species if they share &#8805; 95 ANI and different species if they share < 95 ANI), although a range of cutoffs have been proposed (e.g., 94 ANI, proposed by <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC549018/">Konstantinidis and Tiedje in 2005</a>, 95-96 ANI, proposed by <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC549018/">Richter and Rosselló-Móra (2009)</a> and later supported by <"<a href="https://www.ncbi.nlm.nih.gov/pubmed/24505072">Kim, et al. (2014)</a>). However, <a href="https://www.biorxiv.org/content/early/2017/11/27/225342">recent findings by Jain, et al. (2017) suggest that a species cutoff of 95 is adequate for most bacterial species</a>.

### Why use ANI?

Assuming you have a good quality genome assembly, ANI-based approaches should be able to tell you the species to which your *B. cereus* group genome belongs. While locus-based typing approaches (e.g. *panC*) are extremely valuable for *B. cereus* group isolate characterization, they can be incongruent with whole-genome phylogenetic clade. Also, if you went through the trouble of sequencing an entire *B. cereus* group genome, why throw away data by looking at just a few loci when you can use the WHOLE genome?

### Can I use BTyper to calculate ANI?

As of BTyper version 2.3.0, you can! Just add -\-anib True to your command (by default, BTyper is set to -\-anib False). You can compare your input genome(s) to:

* A database consisting of genomes of 18 published *Bacillus cereus* group species (referred to as the "published" database)

* The published database plus a database consisting of 21 effective *Bacillus* cereus group species that have been proposed in the literature but not published as a novel species (referred to as the "effective" database; 39 species total).

### Which ANIb genome database should I use?

**In nearly all cases, users should use the published genome database (used by default).** This database contains the genomes of all 18 published *B. cereus* group species and is the most accurate representation of our current knowledge of the *B. cereus* group, making it much easier to interpret the results. Furthermore, the published database takes up about half the space of the effective genome database (about 118M, compared to 237M for the effective database) and is much more stable (i.e., it takes a long time to publish a new *B. cereus* group species, so genome additions/removals to/from this database are rare).

A situation in which a user might select the effective species database over the published one is if she/he/they had a genome that was thought to belong to a novel species, and she/he/they wanted to compare the genome to putative novel species proposed in the literature to see if the species has been proposed before.

### If ANI is so great for determining bacterial species, why doesn't BTyper use it by default?

The goal of BTyper is to serve as a **rapid, high-throughput** tool for characterizing *Bacillus cereus* group species *in silico*, and running ANIb on a single draft *B. cereus* group genome using all other default settings increases the analysis time from about 9 seconds to about 1 minute (using the default published ANIb database; using the effective ANIb database, this increases the analysis time to about 2 minutes). For users that don't mind waiting a minute for their results, we recommend using ANIb, as it really is the best way to assign a *B. cereus* group genome to a *B. cereus* group species. However, for some users, this increase may not be trivial. Furthermore, the published and effective ANIb databases take up about 118M and 237M of disk space, respectively, and some users may not have that much to spare. As a result, we have elected to not include it as a default analysis method at this time.

### How does BTyper calculate ANI?

BTyper uses BLAST to calculate ANI; hence, the name ANIb. Because different tools calculate ANI differently (and because we love transparency), the ANIb algorithm implemented in BTyper can be outlined as follows:

0. If -\-draft_genome option is used, concatenate contigs/scaffolds into a single pseudochromosome with a spacer sequence of "NNnnNNnnNNnnNNnn" inserted in between each contig/scaffold (this is done so that the ANIb method is compatible with the other typing methods implemented in BTyper).

1. Fragment the input genome into 1020 bp fragments (see <a href="https://www.ncbi.nlm.nih.gov/pubmed/17220447">Goris, et al. 2007)</a>).

2. BLAST the fragments against each of the reference genomes (18 reference genomes if using the default published database, 39 reference genomes if using the effective database), using the following command structure:

```
NcbiblastnCommandline(query = fragments, db = reference_genome, out = fragments_vs_reference.txt, xdrop_gap_final = 150, evalue = 1e-15, max_target_seqs = 1, dust = "no", outfmt = '"6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qlen"')
```

3. For each hit in the BLAST outfmt 6 output file, do the following:

    a. Calculate a fragment alignment length by subtracting the number of gaps from the alignment length:
    
    frag_alnlen = float(hit.length) - float(hit.gaps)
    
    b. Substract number of mismatches from the fragment alignment length:
    
    frag_alnids = float(frag_alnlen) - float(hit.mismatch)
    
    c. Calculate the fragment coverage by taking the fragment alignment length and dividing by the query length:
    
    frag_anicoverage = (float(frag_alnlen) / float(hit.qlen)) * float(100)
    
    d. Calculate percent nucleotide identity by dividing the number of identical bases by the query length:
    
    frag_anipid = (float(frag_alnids) / float(hit.qlen) * float(100)
    
    e. If this meets cutoff thresholds AND is the top hit for that particular fragment, store it and print it to a final fragments file:
    
    if float(frag_anipid) > float(30) and float(frag_anicoverage) > float(70) and fragment not in used_frags:
        # store fragments and print the following tab-separated line to the final fragments file:
        query_fragment  float(hit.pident) float(hit.qlen) float(hit.bits) float(hit.length) float(hit.gaps)

4. For each final fragments file (i.e. reference genome):

     a. Sum up the hit.length column of the final fragments file and divide by the input genome length to get the ANIb coverage
     
     b. Sum up the hit.pident column of the final fragements file and divide by the number of final fragments (filtered BLAST hits) to get the ANIb value for that input genome/reference genome combination
     
5. Report the reference genome that yields the highest ANIb value for that particular query genome (print to BTyper final results file for that particular query genome)

### I ran BTyper's ANIb method, and my genome was assigned to species "Unknown"; have I discovered a new *B. cereus* group species?

A species assingment of "Unknown" only means that your genome did not share at least 95 ANIb with any *B. cereus* group species in the selected ANIb database (published or effective). How you interpret this requires some thought. Some possible next steps might be:

* Check the quality of your assembly (genome size, coverage, number of contigs, N50, contamination, etc.); a poor-quality assembly can produce artificially low ANIb values

* Check the 16S rDNA gene (-s True) to make sure that your isolate is a member of the *B. cereus* group, i.e. shares at least 97% identity with a member of the *B. cereus* group

* Try using a different tool to calculate ANIb, particularly one that will calculate pairwise ANIb between your genome and a reference genome (to increase speed, BTyper only calculates ANIb in one direction, with the input genome as the query and each database genome as a reference); different tools can produce slightly different ANIb results (e.g. 94.8 ANIb with Tool 1 vs. 95.2 ANI with Tool 2), and BTyper uses 95 ANIb as a hard cutoff for classifying a genome as an "Unknown" species. Try using <a href="http://jspecies.ribohost.com/jspeciesws/">JSpeciesWS</a> and/or <a href="https://github.com/widdowquinn/pyani">pyani</a> to perform pairwise ANIb calculations for your genome relative to *B. cereus* group species type strains. Do all tools produce ANIb values < 95? < 94? Do all tools report adequate coverage values?

* Compare your genome to the genomes of known *B. cereus* group species using the <a href="https://ggdc.dsmz.de/">Genome-to-Genome Distance Calculator (GGDC)</a> to obtain *in silico* DDH values; is this value below 70? What is the probability that it is > 70?

* Check your *panC*, *rpoB*, and MLST results; can your isolate be assigned to a known MLST sequence type? A known *panC* clade? If not, do all of the loci match known *B. cereus* group loci with high identity and coverage? Try constructing a phylogeny using the extracted *panC*, *rpoB*, and MLST loci of your isolate and the *panC*, *rpoB*, and MLST loci of *B. cereus* group type strains; does your isolate cluster among the *B. cereus* group type strains?

* Try running ANIb using BTyper again, this time with the effective ANIb database (-b_db effective); does your isolate share > 95 ANI with a predicted, putative *B. cereus* group species that has been reported in the literature before?


------------------------------------------------------------------------
  
  
## Frequently Asked Questions
  
* **Can I use partial nucleotide sequences (plasmid sequences, MLST genes, *rpoB* alleles, etc.) as input for BTyper?**
  
Sure! You don't have to use whole-genome sequencing data as input; you can technically use any nucleotide sequencing data and treat it as an assembly, as long as it is in fasta or multifasta format (any of the options that require assembly with SPAdes are designed for bacterial genomes). Although it's not necessary, you may want to adjust the options to only perform typing methods you're interested in to make your output easier to read (e.g. if you're interested in detecting virulence and antimicrobial resistance genes in a plasmid sequence, don't waste your time performing other typing methods; just set -\-mlst, -\-rpoB, and -\-panC to False).

* **Can I use whole-genome sequencing data from organisms that don't belong to the *Bacillus cereus* group?**
  
Yes! You can use whole-genome sequencing data from any bacterial species as input; in fact, BTyper's *rpoB* allelic type database provided by Cornell's Food Safety Lab actually contains allelic types for non-*Bacillus* species. Additionally, we've implemented antimicrobial resistance (AMR) gene detection using the ARG-ANNOT database in version 2.0.0 of BTyper, which may be used with any bacterial species. However, your results from certain typing analyses (i.e. MLST using the *B. cereus* group typing scheme) may be completely meaningless. Be extra cautious when interpreting them, and always take identity and coverage values for detected genes into consideration.


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

6. Now, let's run BTyper on our contigs, directly from our Downloads directory. Because this genome is made up of multiple contigs in multifasta format, rather than a single chromosome, we want to make sure we include the -\-draft_genome option in our command. Because we don't know much about this genome, and we want as much information as possible, let's perform all of the default typing methods with their default settings (virulence gene detection using an amino acid database, antimicrobial resistance gene detection using the ARG-ANNOT database, MLST, *rpoB* allelic typing, and *panC* clade typing, all using default thresholds), as well as 16S gene detection, just for fun. To run BTyper, type the following command, and press **Enter**:

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

* **Predicted AMR Genes***
This is a list of antimicrobial resistance (AMR) genes detected in our genome at 50\% coverage and 75\% identity (BTyper's default settings for AMR gene detection). It looks like we have a couple of AMR genes detected at various coverage/identity thresholds.

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
btyper -t seq -i bacillus.fasta -o . -m False -r False -p False -a False -v_db nuc -nuc_p 0 -nuc_q 0
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

* **-a False**
We only want to perform virulence typing, so we are setting antimicrobial resistance gene detection to False.

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

#### ANI/DDH Methods

Goris, Johan, et al. DNA–DNA hybridization values and their relationship to whole-genome sequence similarities. *International Journal of Systematic and Evolutionary Microbiology* 2007 Jan;57(Pt 1):81-91.

Jain, Chirag, et al. High-throughput ANI Analysis of 90K Prokaryotic Genomes Reveals Clear Species Boundaries. bioRxiv 225342; doi: https://doi.org/10.1101/225342.

Kim, Mincheol, et al. Towards a taxonomic coherence between average nucleotide identity and 16S rRNA gene sequence similarity for species demarcation of prokaryotes. *International Journal of Systematic and Evolutionary Microbiology* 2014 Feb;64(Pt 2):346-51.

Konstantinidis, Konstantinos T., and James M. Tiedje. Genomic insights that advance the species definition for prokaryotes. *Proceedings of the National Academy of Sciences of the United States of America* 2005 102(7): 2567–2572.

Meier-Kolthoff, J.P., et al. Genome sequence-based species delimitation with confidence intervals and improved distance functions. *BMC Bioinformatics* 2013 14:60. 

Pritchard, Leighton, et al. Genomics and taxonomy in diagnostics for food security: soft-rotting enterobacterial plant pathogens. *Analytical Methods* 2016 8(1): 12-24.

Richter, Michael and Ramon Rosselló-Móra. Shifting the genomic gold standard for the prokaryotic species definition. *Proceedings of the National Academy of Sciences of the United States of America* 2009 106(45): 19126–19131.

Richter, Michael, et al. JSpeciesWS: a web server for prokaryotic species circumscription based on pairwise genome comparison. *Bioinformatics* 2015 Nov 16. pii: btv681.

#### Antimicrobial Resistance (AMR) Gene Detection

Carroll, Laura M., et al. Whole-Genome Sequencing of Drug-Resistant *Salmonella enterica* Isolates from Dairy Cattle and Humans in New York and Washington States Reveals Source and Geographic Associations. *Applied and Environmental Microbiology* 2017 May 31;83(12).

Fu, L., et al. CD-HIT: accelerated for clustering the next-generation sequencing data. *Bioinformatics* 2012 Dec 1;28(23):3150-2.

Gupta, SK, et al. ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes. *Antimicrobial Agents and Chemotherapy* 2014;58(1):212-20.

Inouye, M., Harriet Dashnow, Lesley-Ann Raven, Mark B Schultz, Bernard J Pope, Takehiro Tomita, Justin Zobel and Kathryn E Holt. SRST2: Rapid genomic surveillance for public health and hospital microbiology labs. *Genome Medicine* 2014 Nov 20;6(11):90.

Lakin, S.M., et al. MEGARes: an antimicrobial resistance database for high throughput sequencing. *Nucleic Acids Research* 2017 Jan 4; 45(Database issue): D574–D580.

#### Plasmid Replicon Detection

Carattoli, Alessandra, et al. *In Silico* Detection and Typing of Plasmids using PlasmidFinder and Plasmid Multilocus Sequence Typing. *Antimicrob Agents and Chemotherapy* 2014 Jul; 58(7): 3895–3903.


#### Tutorial Genomes

Gee, JE, et al. Draft Genome Sequence of *Bacillus cereus* Strain BcFL2013, a Clinical Isolate Similar to G9241. *Genome Announcements* 2014 May 29;2(3).


------------------------------------------------------------------------


Disclaimer: BTyper and BMiner are pretty neat! However, no tool is perfect, and BTyper and BMiner cannot definitively prove whether an isolate is pathogenic or resistant to a particular antimicrobial. As always, interpret your results with caution. We are not responsible for taxonomic misclassifications, misclassifications of an isolate's pathogenic or antimicrobial resistance potential, and/or misinterpretations (biological, statistical, or otherwise) of BTyper and/or BMiner results.











