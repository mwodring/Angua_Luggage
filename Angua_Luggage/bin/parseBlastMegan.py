from .parseBlastXML import getTerms
from ..LuggageInterface import rmaHandler
import argparse, logging, sys, os

logging.basicConfig(stream = sys.stdout, level=logging.DEBUG)
LOG = logging.getLogger(__name__)

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
    parser.add_argument("-co", "--contigsout",
                        help = "Output contigs matching virus species. Default false.",
                        action = "store_true")
    return parser.parse_args()

def main():
    args = parseArguments()
    
    handler = rmaHandler("out", os.path.abspath(args.out_dir), extend = args.extend)
    handler.addFolder("contigs", os.path.abspath(args.contigs))
    handler.findFastaFiles("contigs")
    blast_kind = "Blast" + args.blast_type.upper()
    
    if args.runmegan:
        handler.addFolder("xml", os.path.abspath(args.in_dir))
        handler.blast2Rma(db = args.a2t, blast_kind = blast_kind)
    else:
        handler.addFolder("megan", args.in_dir)
        handler.findRmas(db = args.a2t, blast_kind = blast_kind)
    
    header = ["sample", "contig", "rank", "species"]
    handler.getMeganReport()
    handler.hitsToCSV(header)
    handler.mergeCSVOutput(header)
    
    if args.contigsout:
        handler.hitContigsToFasta()
	
if __name__ == "__main__":
    sys.exit(main())