# fortuna Fly Atlas

fortuna Fly Atlas pipeline is used to catalog novel AS events in scRNASeq data. It proceeds in the steps shown in the figure below.

<p align="center">
<img src="https://user-images.githubusercontent.com/37735817/216773111-e62af577-becd-43fc-aa8d-3e2843a9ca04.png" width=50% height=50%>
</p>

During step 1, we use fortuna with *-bam* and *-alt* flags to identify and quantify novel alternative splicing events and ```cellranger count``` pipeline to perform barcode correction and UMI collapsing. Both tools use the same annotation (.gtf), genome (.fasta) and reads (.fastq) files as input. Fortuna outputs novel AS event catalog and a pseudoalignment .bam file. Cell Ranger will output several files, out of which we will only use the alignment .bam file.  

During step 2, we convert the binary .bam file produced by fortuna into a .sam file using the provided p2g (*pseudo-to-genome-alignment*) tool. Similarly, we convert the binary .bam file produced by Cellranger using ```samtools view``` command. The final part of our pipeline, barcodeAS tool, adjusts fortuna’s catalog of novel AS events according to the information obtained from the two recently created .sam files and saves the result in file .tsv format. 

## Dependencies
fortuna pipeline requires the following software: 
* [fortuna](https://github.com/canzarlab/fortuna)
* [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)
* [samtools](https://www.htslib.org/)

## Installation
After cloning the repository, the p2g and barcodeAS tools have to be compiled using the following commands:

* g++ -std=c++11 -O3 -Wall p2g.cpp -o p2g -I\<path-to-fortuna\>/ext/htslib -L\<path-to-fortuna\>/ext/htslib -lhts -lz -lpthread
* g++ -std=c++11 -O3 -Wall barcodeAS.cpp -o barcodeAS

Both p2g and barcodeAS tools require C++11 to compile and run. 

## Fly Cell Atlas
We evaluated our pipeline on Drosophila single cell data samples from the Fly Cell Atlas downloaded from [ArrayExpress](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10519/sdrf). Reference package was created by following the guidelines from the [FCA github page](https://github.com/FlyCellAtlas/genome_references). 
