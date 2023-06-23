# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 09:55:02 2022

@author: mwodring
"""
import argparse
import sys
from FastaKit import seqHandler
            
def parseArguments():
    parser = argparse.ArgumentParser(description = "Runs 'text search'.")
    parser.add_argument("file",
                       help = "Blast .xml file to search.")
    parser.add_argument("output",
                        help = "Output folder.")
    parser.add_argument("-a", "--all", 
                        help = "Give all hits, not just the top hit for each query.", 
                        action = "store_true")
    parser.add_argument("-st", "--searchterm",
                        help = "Text to look for in the Blast output. Default VIRUS.",
                        default = "VIRUS")
    parser.add_argument("-ml", "--minlength",
                        help = "Minimum contig length to check. Default 200.",
                        default = 200)
    parser.add_argument("-csv", "--outputcsv",
                        help = "Output findings as a .csv file.",
                        action = "store_true")
    #Not yet implemented.
    parser.add_argument("-bl", "--blacklist",
                        help = "Term to exclude. Default 'retrovirus'.",
                        default = "retrovirus")
    parser.add_argument("-c", "--contigs",
                        help = ".fasta file containing the contigs used for the Blast query.",
                        default = False)
    return parser.parse_args()

def runTextSearch(search_term: str, in_file: str, output: str, get_all: bool, minlength: int, output_csv: str, contigs = False):
    if not in_file.endswith(".xml"):
        print("File needs to be in xml. Please check the Blast documentation for outputting this format.")
        quit()

    handler = seqHandler(xml = in_file)
    handler.get_BlastRecord(search_term, minlength, get_all = get_all)
    
    if contigs:
        handler.addFasta(contigs)
        handler.unpackBlastRecord()
        handler.outputFasta(output, search_term)
    
    if output_csv:
        handler.outputCSV(output, "".join(in_file.split("/"))[-1], search_term)
    seqHandler._logger.info("TextSearch complete.")
    handler.outputFullLog(["Fastas", "Alignments", "CSVs"])
 
#Allows running as standalone if need be - it will unpack the arguments itself.
#We're looking for a way to avoid the repetition. 
if __name__ == '__main__':
    options = parseArguments()
    options.searchterm = options.searchterm.upper()
    sys.exit(runTextSearch(options.searchterm, options.file, options.output, options.all, options.minlength, options.outputcsv, options.contigs))
