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
