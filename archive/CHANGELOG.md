# BTyper CHANGELOG

All noteable changes to BTyper will be documented in this file

## [2.0.2] 2017-08-16

### Added
- Included formatted MEGARes database (current version 1.0.1) for antimicrobial resistance gene detection as an option
- Included tblastx as an option for antimicrobial resistance gene detection

### Changed
- Updated ARG-ANNOT database to include all alleles present in the database version 3 (Note: in BTyper version 2.0.1 and 2.0.0, only representative alleles selected by clustering using cd-hit-est were used)

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
