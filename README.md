# MINTIA
     Metagenomic INsertT BIoinformatic Annotation in Galaxy

## Usage
     mintia_assembly.pl -h

## Commands:
     check    - step 0 to check the dependencies
     assemble - step 1 to assemble raw reads...
     annotate - step 2 to annotate filtered and cleaned scaffold(s)


## Description
     Step 0: check the dependencies

     Step 1: assembles raw  reads, looks  for and removes the  cloning vector,  and
     extracts the longest and the most covered contigs. It has been build to handle
     two  types of  raw  reads  as inputs:  single (454, ion torrent reads, ...) or
     paired (Illumina,...) reads.
     This tools is not able to process PacBio or Oxford Nanopore reads.

     Step 2: annotate filtered and cleaned scaffold(s) provided by the step 1.

## Options

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

    -l, --length INT
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
             Zip output name [mintia_assembly.zip]

    -H, --htmlOutput STR
             HTML output name [mintia_assembly.html]

    -L, --logOutput STR
             Log output file name [mintia.log]

    -t, --threads
             number of threads for SPADES [8]

    -h, --help
             Print help
