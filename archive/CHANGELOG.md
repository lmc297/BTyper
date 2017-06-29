# BTyper CHANGELOG

All noteable changes to BTyper will be documented in this file

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
