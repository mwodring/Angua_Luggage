# -*- coding: utf-8 -*-
"""
@author: mwodring
"""

import argparse, sys, os, logging, logging.config
from ..utils import SearchParams

from ..LuggageInterface import blastParser

from pkg_resources import resource_filename

#logging_conf = resource_filename("Angua_Luggage", "data/logging.conf")
#logging.config.fileConfig(logging_conf)
logging.basicConfig(stream = sys.stdout)
LOG = logging.getLogger(__name__)

def parseArguments():
    parser = argparse.ArgumentParser(description = "Runs 'text search'.")
    #Consider .add_mutually_exclusive_group for fastas and nofull.
    #Probably worth breaking into subcommands by now nd running those from main.
    
    #REQUIRED
    parser.add_argument("in_dir",
                       help = "Folder containing .xml file(s).")
    parser.add_argument("out_dir",
                        help = "Output folder.")
    
    #INPUT_FILES
    parser.add_argument("-c", "--contigs",
                        help = ".fasta file containing the contigs used for the Blast query, if you'd like the reads extracted.",
                        default = False)
    parser.add_argument("-r", "--raw",
                        help = "Directory of raw reads if bwa is desired. (For single-ended reads blasts.)")
    
    #SWITCHES
    parser.add_argument("-a", "--get_all", 
                        help = "Give all hits, not just the top hit for each query.", 
                        action = "store_true")
    parser.add_argument("--ictv",
                        help = "ICTV db?",
                        action = "store_true")
    parser.add_argument("-atf", "--acc_to_fa",
                        help = "Output NCBI matches as fastas (for bwa).",
                        action = "store_true")
                        
    #SEARCH_PARAMS
    search_params = parser.add_argument_group("search_params")
    search_params.add_argument("-st", "--searchterm",
                               help = "Text to look for in the Blast output. Default VIRUS. Use a .txt file, one per line, for a whitelist.",
                               default = "virus")
    search_params.add_argument("-ml", "--minlen",
                               help = "Minimum contig length to check. Default 200.",
                               type = int, default = 200)
    search_params.add_argument("-b", "--bitscore",
                                help = "Minimum bitscore to filter on. Default 0 i.e. returns all hits.",
                                type = int, default = 50)
    search_params.add_argument("-bl", "--blacklist",
                                help = "Text to exclude from in the Blast output. Input a .txt file one item per line to exclude multiple terms.",
                                default = "phage")
                        
    parser.add_argument("-e", "--email",
                        help = "Entrez email for NCBI fetching. Required if using NCBI to get accessions.")
    return parser.parse_args()

def getTerms(text_file: str) -> list:
    with open(text_file, "r") as txt:
        data = txt.read()
        return [term.upper() for term in data.split("\n") if term != ""]

def runTextSearch(handler, args):
    whl = getTerms(args.searchterm) if args.searchterm.endswith(".txt") else list(
                   args.searchterm)
    bl = getTerms(args.blacklist) if args.blacklist.endswith(".txt") else list(
                  args.blacklist)
    handler.findBlastFiles(ictv = args.ictv)
        
    queries_parsed, hits = handler.parseAlignments(search_params = 
                                                   SearchParams(whl,
                                                                args.minlen,
                                                                args.bitscore,
                                                                bl),
                                                   get_all = args.get_all)
    handler.hitsToCSV(os.path.splitext(os.path.basename(args.searchterm))[0])
    handler.mergeCSVOutput()
    return f"Total queries checked: {queries_parsed} Total hits found: {hits}"

def getEmail():
    email = input()
    #Now validate.
    return email
    
def main():
    args = parseArguments()
    
    handler = blastParser("xml", args.in_dir)
    handler.addFolder("out", args.out_dir)
    
    LOG.info(runTextSearch(handler, args))
    
    #I could make this a function but for now meh.
    if args.contigs:
        handler.addFolder("contigs", args.contigs)
        handler.findFastaFiles("contigs")
        handler.hitContigsToFasta(by_species = True)
    
    if args.raw and not args.acc_to_fa:
        args.acc_to_fa == True
    
    if args.acc_to_fa:
        while not args.email:
            print("Need an NCBI email for accessions:")
            args.email = getEmail()
        handler.hitAccessionsToFasta(args.email)
    
    if args.raw:
        handler.findFastaFiles("raw")
        tsvs = handler.runBwaTS(args.raw)
        handler.appendMappedToCSV()
        
if __name__ == "__main__":
    sys.exit(main())