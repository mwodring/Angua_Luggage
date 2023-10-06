#Angua_Luggage is a post-pipeline processing helper

##Installing

conda install -c mwodring angua-luggage

I recommend using mamba for this, just install mamba to your base environment and replace 'conda' with 'mamba'. 

You will need a local copy of the following:

- An NCBI protein database.
- An NCBI nucleotide database.
- The Megan na2t and pa2t databases (the most up-to-date one is important).

This toolkit includes a [script](#ICTV) to generate a viral database from ICTV accessions.
If you would like to use Luggage Annotatr, you'll also need a local copy of the pfam database.

##Quick-start

To run Angua with its default settings (bbduk > Trinity > mmseqs2 > Blastn &| Blastx > Megan):

Angua main [RAW_READS] [OUTPUT_DIR] -pa2t [MEGAN PROTEIN DB] -na2t [MEGAN NUC DIR] -nt-db [NUCLEOTIDE BLAST DB] -nr-db [PROTEIN BLAST DB] --cluster -bba [BBDUK ADAPTER FILE]

You can do this from the directory containing the raw directory or using absolute paths to the raw and output directory; both should work.

Angua automatically creates .finished files to track its progress and allow you to pick up where you left off. Remove these if you want to repeat a step for whatever reason.

##Luggage use cases

Angua_Luggage is a Bioinformatics tool bringing together a few useful pieces of software to analyse the output of the Angua pipeline (other pipeline outputs can be used in theory). If you use another pipeline, Luggage might still work for you; as long as you have contigs and Blast files in XML format/.rma6 format Megan files, Luggage should be of use to you.

Luggage has two main functions. One is to quickly summarise pipeline (Blastn/X/Megan) output in .csv format (and output contigs matching desired species, if possible). The other is to automate some basic annotations of contigs: pfam domains and ORFs, alongside coverage. This is to aid in triage in case of several novel viruses, or just a quick way of looking at coverage for diagnostic purposes.

##Inputs to Luggage

In all cases Luggage will need a directory. If you just have one file, please put it in a directory by itself first.

##ICTV

You will need to download the [latest ICTV VMR database](https://ictv.global/vmr) file as an input. There is a link: 'Download current virus metadata resource'.
Place it in a folder and run:
makeICTVdb [FOLDER] [ENTREZ email] 
Run --help for details. It will default to plant hosts only, you may restrict it with other criteria if you wish, or provide an api key for faster retrieval.