# MINTIA - Metagenomic INsertT BIoinformatic Annotation

## Descrition
###### Functional metagenomics is used to understand who is doing what in microbial ecosystems. DNA sequencing can be prioritized by activity-based screening of libraries obtained by cloning and expressing metagenomic DNA fragments in an heterologous host. When large insert libraries are used, allowing a direct access to the functions encoded by entire metagenomic loci sizing several dozens of kbp, NGS is required to identify the genes that are responsible for the screened function. MINTIA is an easy to use pipeline assembling and annotating metagenomic inserts.

###### The assembly module (assemble) assembles, cleans and extracts the longest and most covered contigs for each DNA insert. It handles reads (454, ion torrent,...) or read pairs (Illumina,...). This tools is not able to process PacBio or Oxford Nanopore reads. It sub-selects by default 300X read coverage depending on the sequencing platform and assembles them, removes cloning vector and selects best contigs. Only contigs with length and average depth over given thresholds are kept. The produced HTML report includes a dynamic graphic with contig length and coverage for each sample.

###### The annotation module (annotate) aims at obtaining main gene functions and a functional classification. The pipeline launches at least prokka for ORF detection, generating fasta files for genes and proteins as well as a tabular file containing ORFs description. Depending on additional options selected, contigs and ORFs are aligned against NCBI NR (Non Redundant) as well as SP (SwissProt) and COGss databases. The produced HTML report includes all results and allow to explore annotations based on [igv.js](https://github.com/igvteam/igv.js) an embeddable interactive genome visualization component.


## Table of content
- [Installation](#installation)
	- [Tools dependancies](#tools-dependancies)
	- [Databanks](#databanks)
- [License](#license)
- [Copyright](#copyright)
- [Contact](#contact)

## Installation
This MINITA repository is for command line user.

### Tools dependancies

| Bioinformatics tools | Version |
| ------------- |-------------|
| cross_match | 1.090518 |
| diamond | 0.8.24 |
| Emboss (transeq) | 6.4.0.0 |
| metagene | |
| megan | 5.10.6 |
| rpsblast | 2.2.26 |
| samtools | 1.3.1 |
| spades | |

|Unix tools|Version|
| ------------- |-------------|
|bgzip||
|file|file-5.04|
|grep|GNU grep 2.6.3|
|tabix| 0.2.5 (r964)|
|which||


### Databanks
Reference database are needed to...

## License
GNU GPL v3

## Copyright
2019 INRA

## Contact
support.sigenae@inra.fr
