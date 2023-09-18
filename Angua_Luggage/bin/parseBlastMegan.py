from .parseBlastXML import getTerms

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
    
    #SWITCHES
    parser.add_argument("--ictv",
                        help = "ICTV db?",
                        action = "store_true")
    parser.add_argument("-atf", "--acc_to_fa",
                        help = "Output NCBI matches as fastas (for bwa etc.).",
                        action = "store_true")
    parser.add_argument("-bt", "--blast_type",
                        help = "Type of blast used. N, P or X.",
                        default = "N")
                        
    return parser.parse_args()

def main():
	args = parseArguments()
    handler = rmaHandler("xml", args.in_dir)
	blast_kind = "Blast" + args.blast_type.upper()
	handler.blast2Rma()
	
	

if __name__ == "__main__":
    sys.exit(main())