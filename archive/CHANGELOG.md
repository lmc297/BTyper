# BTyper CHANGELOG

All noteable changes to BTyper will be documented in this file

## [2.3.0] 2018-09-18
### Added
- Average nucleotide identity blast (ANIb) function (-\-anib/-b; set to False by default)
- seq_anib_db directory, which contains a list of genomes of 18 pubished *B. cereus* group species (published.txt), a list of genomes of 21 effective *B. cereus* group species (effective.txt; also lists two genomes that were thought to be potential putative new species but were recently removed from/not included in RefSeq due to small assembly size, which are denoted by a #), and empty published and effective directories for storing genomes to be used with the -\-anib option
- build_btyper_anib_db.py script, which downloads published and/or effective *B. cereus* group genomes listed in the published.txt and/or effective.txt files to be used with -\-anib option
- btyper2matrix.py script, which can be used to aggregate BTyper final results files into a single matrix/text file
- Added *sph* (sphingomyelinase C [Bacillus anthracis str. Ames]) gene to virulence gene databases

### Changed
- Updated ARG-ANNOT AMR database to most recent database (ARG-ANNOT version 4/May 2018)
- Updated MLST database to reflect current PubMLST *B. cereus* database (1795 STs)
- Updated README to reflect addition of ANIb and btyper2matrix.py

## [2.2.2] 2018-05-03
### Added
- Added new default 7-clade *panC* database containing only type strain *panC* genes for 9 *B. cereus* group species; we have found that this correlates better with whole-genome phylogenetic clade based on the current 7-clade system than the legacy database, which includes non-type strains

### Changed
- Changed default *panC* database from "legacy" to "type"

## [2.2.1] 2018-04-21
### Changed
- Forgot to upload virulence gene nucleotide sequence database for use with -v_db nuc option in version 2.2.0; fixed this
- Updated MLST database to reflect current PubMLST *B. cereus* database

## [2.2.0] 2018-03-21
### Added
- Added new default, location-based antimicrobial resistance (AMR) gene pruning method, which can take into account the position of a gene in the genoe: genes overlapping at more than some proportion threshold (default = 0.7) are treated as alleles of the same gene and only the top hit is reported. This is more accurate than the previous cluster-based method, which incorrectly reported tetracycline resistance genes at the same position in several *Salmonella* genomes as being different genes, and did not work with plasmid replicons
- Added most-recent PlasmidFinder plasmid replicon database as an option to be used with --amr method

### Changed
- Updated MLST database to reflect current PubMLST *B. cereus* database
- Updated README to account for location-based AMR gene detection parameters and PlasmidFinder replicon detection options; updated Tutorial 2 to account for inclusion of AMR gene detection

## [2.1.0] 2018-02-15
### Added
- Added new panC database containing panC genes from 18 *B. cereus* group species type strains: 9 from the original 9 species, and 9 from proposed new species in Liu, et al. 2017 (https://www.ncbi.nlm.nih.gov/pubmed/28792367); rather than numbered clades (clade1, clade2, ...clade7), this database uses species (cladeAnthracis, cladeCereus, etc.)
- Added new 16S rDNA database containing 16S rDNA genes from 18 *B. cereus* group species type strains: 9 from the original 9 species, and 9 from proposed new species in Liu, et al. 2017 (https://www.ncbi.nlm.nih.gov/pubmed/28792367)
- Added -panC_db/--panC_database option to specify which panC database BTyper should use: -panC_db "latest" for the new panC database described above, or -panC_db "legacy" for the original database and 7-clade typing scheme
- Added -s_db/--s_database option to specify which 16S rDNA database BTyper should use: -s_db "latest" for the 18-species database described above, or -s_db "legacy" for the original 9-species database

### Changed
- Updated MLST database to reflect current PubMLST *B. cereus* database
- Fixed typo in help message; help message displayed a default value of "50" for a minimum antimicrobial resistance gene identity value, but the default value was 75 as stated in the manual (BTyper was using a minimum identity of 75 as default, but the help message displayed 50 to the user when "btyper -h" was called)
- Updated homebrew formula depends_on "python" line (depends_on :python is no longer used)
- Updated README to account for new species (*panC* and 16S rDNA database changes), as well as an update to Tutorial 1 to account for the addition of AMR gene detection in BTyper v. 2.0.0

## [2.0.3] 2017-11-14

### Changed
- Fixed ARG-ANNOT antimicrobial resistance gene database (clusters 61 and 64) to have parentheses around (Phe)
- Replaced “(“ and “)” in virulence database with “[“ and “]” for *bpsB* gene
- Updated MLST database to reflect current PubMLST *B. cereus* database
- Fixed btyper_path variable bug so users who run from source file (rather than install using Homebrew) do not need to run the program within the btyper-2.0.3 directory
- Updated README to include instructions on how to install pip and biopython

### Removed
- Removed biopython installation from brew formula (lmc297/homebrew-btyper; see README on how to install using pip)

## [2.0.2] 2017-08-16

### Added
- Included formatted MEGARes database (current version 1.0.1) for antimicrobial resistance gene detection as an option
- Included tblastx as an option for antimicrobial resistance gene detection

### Changed
- Updated ARG-ANNOT database to include all alleles present in the database version 3 (Note: in BTyper version 2.0.1 and 2.0.0, only representative alleles selected by clustering using cd-hit-est were used)
- Changed "/" symbols in ARG-ANNOT fasta headers to "_"
- Updated README to reflect addition of MEGARes database and the addition of tblastx as an AMR gene detection algorithm; updated references and disclaimer

## [2.0.1] 2017-08-09

### Added
- Updated formatted (see version 2.0.0) ARG-ANNOT database to reflect ARG-ANNOT version 3 (March 2017)
- Added bpsXABCDEFGH operon (alternative capsule protein-encoding genes for some anthrax-causing B. cereus s.s. strains) to virulence gene databases

## [2.0.0] 2017-06-29

### Added
- Antimicrobial resistance gene detection function (--amr/-a) using the ARG-ANNOT database and nucleotide BLAST (blastn)
- ARG-ANNOT antimicrobial resistance gene database, clustered using cd-hit-est and a sequence identity cutoff of 0.8

### Changed
- Re-formatted virulence gene database prefixes from "<cluster#>___" to "<cluster#>vir___"
- Added 3 plcR alleles to virulence gene database (plcR from B. anthracis str. Ames, B. cereus ATCC 14579, and B. thuringiensis serovar konkukian str. 97-27)
- Updated MLST database to reflect current PubMLST B. cereus database
- Updated BTyper citation information to reflect paper's acceptance to AEM
- Updated README to reflect addition of antimicrobial resistance gene detection & fixed typos

## [1.0.0] 2017-05-14

### Added
- BTyper initial commit (virulence gene detection, panC clade typing, rpoB allelic typing, MLST, and 16S typing)
