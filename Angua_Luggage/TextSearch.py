# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 09:55:02 2022

@author: mwodring
"""
import argparse
import sys
import os
from .FastaKit import seqHandler, logFileOutput
from .utils import new_logger
            
def parseArguments():
    parser = argparse.ArgumentParser(description = "Runs 'text search'.")
    #Consider .add_mutually_exclusive_group for fastas and nofull.
    #Probably worth breaking into subcommands by now nd running those from main.
    
    parser.add_argument("file",
                       help = "Blast .xml file to search.")
    parser.add_argument("output",
                        help = "Output folder.")
    parser.add_argument("-a", "--all", 
                        help = "Give all hits, not just the top hit for each query.", 
                        action = "store_true")
    parser.add_argument("-st", "--searchterm",
                        help = "Text to look for in the Blast output. Default VIRUS. Use a .txt file, one per line, for a whitelist.",
                        default = "virus")
    parser.add_argument("-ml", "--minlength",
                        help = "Minimum contig length to check. Default 200.",
                        type = int, default = 200)
    parser.add_argument("-b", "--bitscore",
                        help = "Minimum bitscore to filter on. Default 0 i.e. returns all hits.",
                        type = int, default = 0)
    parser.add_argument("-csv", "--outputcsv",
                        help = "Output findings as a .csv file.",
                        action = "store_true")
    parser.add_argument("-bl", "--blacklist",
                        help = "Text to exclude from in the Blast output. Input a .txt file one item per line to exclude multiple terms.",
                        default = "phage")
    parser.add_argument("-c", "--contigs",
                        help = ".fasta file containing the contigs used for the Blast query, if you'd like the reads extracted.",
                        default = False)
    parser.add_argument("--ictv",
                        help = "ICTV db?",
                        action = "store_true")
    parser.add_argument("-nf", "--nofull",
                        help = "Suppress output of a combined csv of all the samples in the input folder.",
                        action = "store_true")
    parser.add_argument("-f", "--fastas",
                        help = "Output NCBI matches as fastas (for bwa). Overrides csv setting to True.",
                        action = "store_true")
    parser.add_argument("-e", "--email",
                        help = "Entrez email for NCBI fetching. Required if using NCBI to get accessions.")
    parser.add_argument("-r", "--raw",
                        help = "Directory of raw reads if bwa is desired. (For single-ended reads blasts.)")
    return parser.parse_args()

def unpackTxt(text_file: str) -> list:
    with open(text_file, "r") as txt:
        data = txt.read()
        new_list = [term.upper() for term in data.split("\n") if term != ""]
        return new_list
    
def runTextSearch(in_file: str, output: str, **kwargs):
    handler = None
    if kwargs["search_term"].endswith(".txt"):
        search_term = unpackTxt(kwargs["search_term"])
            
    if kwargs["blacklist"].endswith(".txt"):
        blacklist = unpackTxt(kwargs["blacklist"])

    handler = seqHandler(xml = in_file)
    handler.addFolder("out", output)
    handler.get_BlastRecord(list(search_term), list(blacklist), kwargs["minlength"], kwargs["get_all"], kwargs["bitscore"], kwargs["ictv"])
    
    if kwargs["contigs"]:
        sample_name = os.path.splitext(os.path.basename(in_file))[0]
        contig_file = os.path.join(kwargs["contigs"], f"{sample_name}.fasta")
        handler.addFasta(contig_file)
        handler.unpackBlastRecord()
        #Need to change this to work with the individual search terms.
        handler.outputFasta("vir", by_species = True)
    
    if kwargs["output_csv"] or kwargs["fastas"]:
        handler.outputCSV(os.path.splitext(os.path.basename(in_file))[0], "vir")
    #handler.outputFullLog(["Alignments"])
    return handler
    
#Should modify this to search a whole folder or else mergeCSVOutput will run repeatedly.
def runTextSearchOnFolder(in_dir: str, output: str, **kwargs):
    handler = None
    if not os.path.exists(output):
        os.mkdir(output)
    logFileOutput(output, "Luggage")
    
    for file in os.scandir(in_dir):
        if file.name.endswith(".xml"):
            handler = runTextSearch(in_file = os.path.join(in_dir, file.name), output = output, **kwargs)
    
    if not handler:
        #Do this as an exception instead.
        print("No .xml found. Wrong folder?")
        quit()
    
    if not kwargs["no_full"]:
        merged_csv = handler.mergeCSVOutput()
    
    if kwargs["fastas"]:
        handler.TextSearchToFasta(output, kwargs["email"]) 
        if kwargs["raw"]:
            bwa_outdir = os.path.join(output, "bwa")
            handler.runBwaTS(raw_dir = kwargs["raw"], bwa_dir = bwa_outdir)
            handler.addToCSV(merged_csv)
        #Add reads mapped to .csv file using pandas.
        
#Allows running as standalone if need be - it will unpack the arguments itself.
#We're looking for a way to avoid the repetition. 
if __name__ == '__main__':
    options = parseArguments()
    logger = new_logger(name = __name__, output_dir = options.output)
    sys.exit(runTextSearchOnFolder(options.file, options.output, 
                                   search_term = options.searchterm, get_all = options.all, 
                                   minlength = options.minlength, output_csv = options.outputcsv, 
                                   bitscore = options.bitscore, ictv = options.ictv, 
                                   no_full = options.nofull, fastas = options.fastas, 
                                   email = options.email, blacklist = options.blacklist, 
                                   raw = options.raw, contigs = options.contigs))
