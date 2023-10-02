#Angua_Luggage is a post-pipeline processing helper

##Luggage use cases

Angua_Luggage is a Bioinformatics tool bringing together a few useful pieces of software to analyse the output of the Angua pipeline. If you use another pipeline, Luggage might still work for you; as long as you have contigs and Blast files in XML format/.rma6 format Megan files, Luggage should be of use to you.

Luggage has two main functions. One is to quickly summarise pipeline output in .csv format (and output contigs matching desired species, if possible). The other is to automate some basic annotations: pfam domains and ORFs, alongside coverage. This is to aid in triage in case of several novel viruses, or just a quick way of looking at coverage for diagnostic purposes.

##Inputs to Luggage

In all cases Luggage will need a directory. If you just have one file, please put it in a directory by itself first.

