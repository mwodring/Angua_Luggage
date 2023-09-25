from .parseBlastXML import getTerms
from ..LuggageInterface import rmaHandler
import argparse

def parseArguments():
    parser = argparse.ArgumentParser(description = "Runs Megan on .xml files and generates a report.")

    #REQUIRED
    parser.add_argument("in_dir",
                       help = "Folder containing .xml file(s).")
    parser.add_argument("out_dir",
                        help = "Output folder.")
    parser.add_argument("contigs",
                        help = ".fasta file containing the contigs used for the Blast query.")
    
    #INPUT_FILES
    parser.add_argument("-r", "--raw",
                        help = "Directory of raw reads if bwa is desired. (For single-ended reads blasts.)")
    parser.add_argument("-a2t",
                        help = "Location of database file for Megan accessions to taxa.")
    
    #OTHER
    parser.add_argument("--ictv",
                        help = "ICTV db?",
                        action = "store_true")
    parser.add_argument("-atf", "--acc_to_fa",
                        help = "Output NCBI matches as fastas (for bwa etc.).",
                        action = "store_true")
    parser.add_argument("-bt", "--blast_type",
                        help = "Type of blast used. N, P or X.",
                        default = "N")
    parser.add_argument("-m", "--runmegan",
                        help = "Start from blast files by running Megan.",
                        action = "store_true")                    
    parser.add_argument("-ex", "--extend",
                        help = "Number of underscores to remove from sample names.",
                        type = int, default = 1)
    return parser.parse_args()

def main():
    args = parseArguments()
    handler = rmaHandler("out", args.out_dir, extend = args.extend)
    if args.runmegan:
        handler.addFolder("xml", args.in_dir)
        blast_kind = "Blast" + args.blast_type.upper()
        handler.addFolder("contigs", args.contigs)
        handler.findFastaFiles("contigs")
        handler.blast2Rma(db = args.a2t, blast_kind = blast_kind)
    else:
        handler.addFolder("megan", args.in_dir)
    
    handler.getMeganReport()

	#Hook into the existing csv stuff and continue from there.
	
if __name__ == "__main__":
    sys.exit(main())