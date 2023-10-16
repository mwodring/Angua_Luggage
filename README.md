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

### Parameters 

You may also use trimmed reads as an input, Angua will refer to them as 'raw' as it doesn't trim them, but QC can still be performed as normal.

-a trinity assembles using Trinity, -a spades uses spades (it will move and rename scaffolds.fa files). Single or paired reads are inferred from filenames (R1, R2).

-sort [blastn length] [blastx length] overrides the automatic 200 and 1000 contig sorting, and clusters the longer set of reads as usual. 

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

- One (**parseBlast** and **parseMegan**) is to quickly summarise pipeline (Blastn/X/Megan) output in .csv format (and output contigs matching desired species, if possible). 
- The other (**Annotatr**) is to automate some basic annotations of contigs: pfam domains, cytoplasmic/non-cytoplasmic/transmembrane domains and ORFs, alongside coverage. This is to aid in triage in case of several novel viruses, or just a quick way of looking at coverage for diagnostic purposes.

### Additional databases

#### Annotatr 

##### Pfam

For Annotatr, you will need a local copy of the [current pfam database](https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/) HMM files.

Follow the [pfam_scan docs](https://github.com/aziele/pfam_scan) if you get lost. 

Broadly:

```
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
```

```
mkdir pfamdb
gunzip -c Pfam-A.hmm.dat.gz > pfamdb/Pfam-A.hmm.dat
gunzip -c Pfam-A.hmm.gz > pfamdb/Pfam-A.hmm
rm Pfam-A.hmm.gz Pfam-A.hmm.dat.gz
```

```
hmmpress pfamdb/Pfam-A.hmm
```

Use this pfamdb/ folder as input to annotatr db_dir.

##### Phobius

Download [Phobius](http://software.sbc.su.se/cgi-bin/request.cgi?project=phobius) phobius101_linux.tar.gz.

With your Angua-Luggage environment active:

```
phobius-register phobius101_linux.tar.gz
```

To add it to the present conda environment.

### Inputs to Luggage

In all cases Luggage will need a directory. If you just have one file, please put it in a directory by itself first.

### parseBlast

```
parseBlast [blast_dir] [out_dir] -wl [whitelist.txt] -bl [blacklist.txt] -r [trimmed_reads] -bt [N/X/P] -c [contigs] -e [NCBI email]
```

whitelist.txt can be a single word or a .txt file of words, one line per word. The same goes for blacklist. Luggage will take any items matching at least one of the whitelist terms, and exclude it if any of the blacklist terms are present.

By giving parseBlast trimmed reads, it will map reads to NCBI accessions found in the Blast file. 

By giving it contigs, it will label contigs by their hits and create .fasta files for you.

If you used the database generated by fetchICTV, use the --ictv flag as the .xml is different.

### parseMegan

To create a .csv file of Megan hits (virus only) and output contigs matching those hits:

```
parseMegan [in_dir] [out_dir] [contigs] -r [trimmed_reads] -a2t [(n/p)a2t)] -bt [N/P/X] -co
```

in_dir can be blast files (use -m to run Megan on them) or a directory of .rma6 files.

### Annotatr (ORFs and protein families)

```
annotatr [in_dir] [out_dir] [pfam_db] -t [trimmed_reads]
```

Will generate graphs of ORFs and pfam hits, along with Phobius-predicted cytoplasmic, non-cytoplasmic and transmembrane domains, next to coverage determined by the same pipeline as back-mapper. 

I would recommend using directories output by the parseBlast or parseMegan pipeline rather than contigs directly from e.g. Angua.