# MINTIA - Metagenomic INsertT BIoinformatic Annotation

## Description
> Functional metagenomics is used to understand who is doing what in microbial ecosystems. DNA sequencing can be prioritized by activity-based screening of libraries obtained by cloning and expressing metagenomic DNA fragments in an heterologous host. When large insert libraries are used, allowing a direct access to the functions encoded by entire metagenomic loci sizing several dozens of kbp, NGS is required to identify the genes that are responsible for the screened function. MINTIA is an easy to use pipeline assembling and annotating metagenomic inserts.

> The assembly module (assemble) assembles, cleans and extracts the longest and most covered contigs for each DNA insert. It handles reads (454, ion torrent,...) or read pairs (Illumina,...). This tools is not able to process PacBio or Oxford Nanopore reads. It sub-selects by default 300X read coverage depending on the sequencing platform and assembles them, removes cloning vector and selects best contigs. Only contigs with length and average depth over given thresholds are kept. The produced HTML report includes a dynamic graphic with contig length and coverage for each sample.

> The annotation module (annotate) aims at obtaining main gene functions and a functional classification. The pipeline launches at least prokka for ORF detection, generating fasta files for genes and proteins as well as a tabular file containing ORFs description. Depending on additional options selected, contigs and ORFs are aligned against NCBI NR (Non Redundant) as well as SP (SwissProt) and COGss databases. The produced HTML report includes all results and allow to explore annotations based on [igv.js](https://github.com/igvteam/igv.js) an embeddable interactive genome visualization component.


## Table of content
- [Installation](#installation)
	- [Tools dependancies](#tools-dependancies)
	- [Databanks](#databanks)
	- [Install](#install)
- [Run MINTIA](#run-mintia)
	- [Check tools dependancies](#check-tools-dependancies)
	- [Assemble](#assemble)
	- [Annotate](#annotate)
- [License](#license)
- [Copyright](#copyright)
- [Contact](#contact)

## Installation
This MINITA repository is for command line user.

#### Tools dependancies

| Bioinformatics tools | Version || Unix tools | Version |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| **`cross_match`** | `1.090518` || **`bgzip`** | |
| **`spades`** | `v3.13.0` || **`file`** | `file-5.04` |
| **`prokka`** | `1.13.3` || **`grep`** | `GNU grep 2.6.3` |
| **`diamond`** | `0.8.24` || **`which`** |  |
| **`megan`** | `5.10.6` || **`xvfb-run`** | |
| **`samtools`** | `1.3.1` || **`tabix`** | `0.2.5 (r964)` |
| **`rpsblast`** | `2.2.26` || | |


#### Databanks
Reference databases are needed to the "annotate" module.
- NR 
- Uniprot/Swissprot
- COGs 
> Create COGs DB
```sh
$ wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd//cdd.tar.gz
$ tar -xvzf cdd.tar.gz
$ makeprofiledb -title COG.3-28-17 -in Cog.pn -out Cog.v3-28-17 -threshold 9.82 -scale 100.0 -dbtype rps -index true
```

#### Install
Clone this repository:
```sh
$ git clone --recursive https://github.com/Bios4Biol/MINTIA.git
```
Use conda to install the third party software:
```sh
$ cd MINTIA
$ conda env create --file environment.yaml
```
Two dependencies will not be installed by conda and must be installed "manually":
- cross_match required for the assemble module (step1): [cross_match](http://www.phrap.org/consed/consed.html#howToGet)
- MEGAN5 (optional) for the annotate module (step2): [megan5](https://software-ab.informatik.uni-tuebingen.de/download/megan5/welcome.html)

## Run MINTIA

#### Check tools dependancies

```
$ ./mintia.pl check -h
Name:
     mintia.pl - Fosmid assembly and annotation pipeline.

Check Synopsis:
     mintia.pl check

Check Options:
    -h, --help
             Print help
```

#### Assemble

```
$ ./mintia.pl assemble -h
Name:
     mintia.pl - Fosmid assembly and annotation pipeline.

Assemble Synopsis:
     mintia.pl assemble --input FILE[S] --vectorSeq FILE --dirOutputs STR

Assemble Options:
    -i, --input FILE[S]
             Fastq(.gz) file(s)
             For each sample one OR two fastq file must be provided:
             - Paired data must contain R1/R2: R[12].f[ast]q[.gz]
               Ex: sampleName1_R1.fastq sampleName1_R2.fastq
                   sampleName2_R1.fq.gz sampleName2_R2.fq.gz
             - Single data
               Ex: sampleName.f[ast]q[.gz]

    -v, --vectorSeq FILE
             Path to the vector fasta file

    --length INT
             Fosmid's expected length [40000]

    --minimalContigLength INT
             Contig's minimum length [1000]

    --minimalContigDepth INT
             Contig's minimum depth [8]

    -c, --maxDepth INT
             Coverage, maximum depth use to filter input reads [300]

    -d, --dirOutputs STR
             Path to the outputs directory

    -Z, --zipOutput STR
             Zip output name [mintia_assemble.zip]

    -H, --htmlOutput STR
             HTML output name [mintia_assemble.html]

    -L, --logOutput STR
             Log output file name [mintia_assemble.log]

    -t, --threads
             number of threads for SPADES [8]

    -h, --help
             Print help
```
Example on test data:
```
./mintia.pl assemble -t 1 --input Data/Input/Assemble/BifidoAdolescentis.simul*gz --vectorSeq Data/Input/Assemble/pCC1FOS.fasta -len 40000 -c 300 -d Data/Output/Assemble/
```

#### Annotate

```
$ ./mintia.pl annotate -h
Name:
     mintia.pl - Fosmid assembly and annotation pipeline.

Annotate Synopsis:
     mintia.pl annotate --input FILE --dirOutputs STR

Annotate Options:
    -i, --input FILE
             Fasta(.gz) file

    -s, --separator CHAR [#]
             Which separator allows retreiving the fosmid name
             This separator will be also use to create ORF id
             Example: >fosmidName1#contig1... |
                      >fosmidName1#contig2... | => fosmidName1
                      >fosmidName1#contig3... |

    -F, --FunctionalAndTaxonomic
             Run functional and taxonomic annotations

    -e, --evalue FLOAT
             Maximum diamond e-value to report alignments [10e-8]

    --query-cover INT
             Minimum diamond query cover% to report an alignment [50]

    -M, --Megan FILE
             Run MEGAN - A license file must be provided

    -C, --Cog FILE
             Run annotations with COGs, DB COGs path file

    -c, --cMaxEvalue FLOAT
             Max Evalue for Rps-Blast COGs filtering [10e-8]

    -S, --SubmissionFiles
             Build submission files

    -D, --DiamondAgainstPrivateDB FILE
             Run diamond against your own protein reference fasta file

    -t, --threads INT
             Number of threads for Blast [8]

    -d, --dirOutputs STR
             Path to the outputs directory

    -H, --htmlOutput STR
             HTML output name [mintia_annotate.html]

    -L, --logOutput STR
             Log output file name [mintia_annotate.log]

    -k, --keepTmpFiles
             Keep temporary files

    -h, --help
             Print help
```

## License
GNU GPL v3

## Copyright
2019 INRA

## Contact
support.sigenae@inra.fr
