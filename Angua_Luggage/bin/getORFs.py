# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 10:19:35 2022

@author: mwodring
"""

import argparse, sys, os, logging, logging.config
from ..LuggageInterface import Annotatr

def parseArguments():
    parser = argparse.ArgumentParser(description = 
    """Runs pfam on ORFs from contigs and generates plots and/or GFF3 from the 
    results.""")
    
    #REQUIRED
    parser.add_argument("in_dir",
                        help = "Directory containing contigs to scan.")
    parser.add_argument("out_dir",
                        help = """Directory to place resulting files in. 
                                If it doesn't exist, it will be made for you.""")
    parser.add_argument("db_dir",
                        help = "Directory containing the local pfam database.")
          
    #OPTIONS
    parser.add_argument("-t", "--trimmed",
                        help = "Directory of trimmed reads to map to the contigs.",
                        required = False)
    parser.add_argument("-np", "--no_plot",
                        help = "Turn off outputting a plot.",
                        action = "store_true")
    parser.add_argument("--gff3",
                        help = "Output annotations as a GFF3 file.",
                        action = "store_true")
    parser.add_argument("-ex", "--extend",
                        help = "Number of underscore extensions to the file. -Angua specific-",
                        required = False,
                        type = int)
                        
    return parser.parse_args()
    
def main():
    logging.basicConfig(stream = sys.stdout)
    LOG = logging.getLogger(__name__) 

    args = parseArguments()
    plotter = Annotatr("contigs", args.in_dir, args.extend)
    plotter.getORFs(args.out_dir)
    plotter.runPfam(args.db_dir)
    plotter.getAnnotations(no_plot = args.no_plot, 
                           gff3 = args.gff3, 
                           trimmed_dir=args.trimmed)
    
if __name__ == "__main__":
    sys.exit(main())
