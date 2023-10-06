# Angua_Luggage is an HTS Bioinformatics pipeline and toolkit

## Angua

### Installing

I recommend using mamba, just install mamba to your base environment and run:
```
mamba install -c mwodring angua-luggage
```
You will need a local copy of the following to run the main Angua pipeline:

- An NCBI protein database.
- An NCBI nucleotide database.
- The Megan na2t and pa2t databases (the most up-to-date one is important). [Need help?](#Megan)

This toolkit includes a [script](#ICTV) to generate a viral database from ICTV accessions.

### Quick-start

To run Angua with its default settings (qc > bbduk > qc > Trinity > mmseqs2 > Blastn & Blastx > Megan):
```
Angua main [RAW_READS] [OUTPUT_DIR] -pa2t [MEGAN PROTEIN DB] -na2t [MEGAN NUC DB] -nt-db [NUCLEOTIDE BLAST DB] -nr-db [PROTEIN BLAST DB] --cluster -bba [BBDUK ADAPTER FILE]
```
You can do this from the directory containing the raw directory or using absolute paths to the raw and output directory; both should work. 

Make sure paired reads end with R1_001.fastq.gz and R2_001.fastq.gz (Illumina default). Single-ended reads (R1_001.fastq.gz only) will also work, but for short (<50 bp) reads, please use -a spades or you may not get any contigs.

Angua creates .finished files to track its progress and allow you to pick up where you left off. Remove these if you want to repeat a step for whatever reason.

### Megan dbs

Go to the Megan 6 [downloads](https://software-ab.cs.uni-tuebingen.de/download/megan6/welcome.html). You'll want the files starting with megan-map (pa2t) and megan-nucl (na2t). 

### Back-mapper

Angua back-mapper maps a directory of reads to a directory of fasta files using bwa-mem2. The output is histograms of coverage, indexed bam files, and .tsvs.

Output is one histogram/bam/etc. file for each sample/fasta combination. In this case, fasta means one >, not the .fasta. You can input several individual sequences per fasta, and each .fasta becomes a subdirectory in the resulting folder.

For example, inputting a directory with two samples and a directory with two .fastas, one with one sequence, one with two sequences.

```
out/
├─ Histograms/
├─ MonopartiteVirus/
│  ├─ Sample1_Genome
│  ├─ Sample2_Genome
├─ BipartiteVirus/
│  ├─ Sample1_Segment1
│  ├─ Sample1_Segment2
│  ├─ Sample2_Segment1
│  ├─ Sample2_Segment2
in/
├─ raw_reads/
│  ├─ Sample2_R1.fq
│  ├─ Sample1_R1.fq
├─ in_fastas/
│  ├─ BipartiteVirus.fasta
│  ├─ MonopartiteVirus.fasta
```

### Taxa-finder

Uses Nhmmer to map contigs to plant hosts. You will need databases for these. [NOTE: Mog isn't too familiar with this, so details will come soon.]

## ICTV

You will need to download the [latest ICTV VMR database](https://ictv.global/vmr) file as an input. There is a link: 'Download current virus metadata resource'.

**Script last run by Mog**: 6/8/2023. If it breaks for you, please contact me and I'll fix it; the spreadsheet probably changed format.

Place it in a folder and run:
```
makeICTVdb [FOLDER] [ENTREZ email] 
```
Run --help for details. It will default to plant hosts only, you may restrict it with other criteria if you wish, or provide an api key for faster retrieval.

## Luggage

### Use cases

Angua_Luggage is a Bioinformatics tool bringing together a few useful pieces of software to analyse the output of the Angua pipeline (other pipeline outputs can be used in theory). If you use another pipeline, Luggage might still work for you; as long as you have contigs and Blast files in XML format/.rma6 format Megan files, Luggage should be of use to you.

Luggage has two main functions. 

- One (**parseBlast** and **parseBlastMegan**) is to quickly summarise pipeline (Blastn/X/Megan) output in .csv format (and output contigs matching desired species, if possible). 
- The other (**Annotatr**) is to automate some basic annotations of contigs: pfam domains and ORFs, alongside coverage. This is to aid in triage in case of several novel viruses, or just a quick way of looking at coverage for diagnostic purposes.

### Additional databases

For Annotatr, you will need a local copy of the current pfam database. 

### Inputs to Luggage

In all cases Luggage will need a directory. If you just have one file, please put it in a directory by itself first.