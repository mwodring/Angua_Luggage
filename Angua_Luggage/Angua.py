import os, sys, shutil, multiprocessing
import .exec_Angua
import argparse
import datetime
import inspect
import pysam

from pathlib import Path
from Bio import SeqIO
from Bio.Blast import NCBIXML

from .Luggage import fileHandler, toolBelt
from .exec_utils import *

#Doesn't use a file right this second but it will.
logging.basicConfig(stream = sys.stdout)
LOG = logging.getLogger(__name__)

class Angua(fileHandler):
	def __init__(self):
		super().__init__()
		#This is hard-coded for now but I'll change it soon.
		self.initPipelineFolders()
	
	def optional_folders(func):
		#This function looks at the input args and checks if it needs to override
		#any defaults.
        @functools.wraps(func)
        def wrapper_optional_folders(self, *args, **kwargs):
            for arg_name, arg_contents in kwargs.items():
                if arg_name.endswith("dir") and arg_contents:
                    dir_kind = "_".join(arg_name.split("_")[:-1])
                    self.addFolder(dir_kind, arg_contents)
            return func(self, *args, **kwargs)
        return wrapper_optional_folders
        
	def extendFolderMultiple(self, orig_dir_kind: str, 
								   dir_kinds: list, dir_names: list):
		#Zip?
		for dir_kind, dir_name in zip(dir_kinds, dir_names):
			self.extendFolder(orig_dir_kind, dir_kind, dir_name)
			
	def initPipelineFolders(self):
		self.extendFolder("out", "results", "Results")
		self.extendFolder("out", "qc", "QC")
		self.extendFolderMultiple("QC", ["QC_raw_F", "QC_raw_M"],
										["FastQC", "MultiQC"])
		self.extendFolderMultiple("QC", ["QC_trimmed_F", "QC_trimmed_M"],
										["FastQC", "MultiQC"])
		self.extendFolder("out", "trimmed", "Bbduk")
		self.extendFolder("out", "assembly", options.assembler)
		self.extendFolder("out", "contigs", "Contigs")
		self.extendFolder("out", "clustered", "Mmseqs2")
		self.extendFolder("out", "unmapped", "Unmapped")
		self.extendFolderMultiple("unmapped", 
								 ["unm_bwa", "unm_reads",  "unm_clust"
								  "unm_blastn"],
								 ["BwA", "Reads", "Mmseqs2", "Blastn"])
		self.extendFolderMultiple("out", ["blastn", "blastx"], 
										 ["Blastn", "Blastx"])
		self.extendFolder["out", "megan", "Megan"]
		self.extendFolderMultiple("megan", ["meg_blastn", "meg_blastx"],
										   ["blastn", "blastx"])
										   
	@optional_folders
	def pre_qc_checks(self, qc_threads: int, raw_dir = "", 
							QC_raw_F_dir = "", QC_raw_M_dir = ""):
		# Run FastQC and MultiQC on the supplied directory
		fastqc_dir = self.getFolder("QC_raw_F")
		exec_Angua.fastQC(qc_threads, self.getFolder("raw"), fastqc_dir)
		exec_Angua.multiQC(fastqc_dir, self.getFolder("QC_raw_M"))
	
	@optional_folders
	def post_qc_checks(self, qc_threads: int, 
							 QC_trimmed_F_dir = "", QC_trimmed_M_dir = ""):
		exec_Angua.fastQC(qc_threads, self.getFolder("trimmed"), 
									  self.getFolder("QC_trimmed_F")
		exec_Angua.multiQC(self.getFolder("QC_trimmed_F"),
						   self.getFolder("QC_trimmed_M")	
	
	@optional_folders
	def run_bbduk(self, min_len: int, adapters: str, min_q: int, 
						raw_dir = "", trimmed_dir = ""): 
		trimmed_dir = self.getFolder("trimmed")
		out_dir = self.getFolder("raw")
		for file_R1 in self.getFiles("raw", "_R1.fasta"):
			file_R2 = file.replace("_R1", "_R2")
			sample_name_R1 = getSampleName(file_R1)
			sample_name_R2 = sample_name_R1.replace("_R1", "_R2")

			trimmer_input_R1 = os.path.join(raw_dir, file_R1)
			trimmer_input_R2 = os.path.join(raw_dir, file_R2)
			trimmer_output_R1 = os.path.join(trimmed_dir, file.replace('_L001_R1_001', '_R1'))
			trimmer_output_R2 = os.path.join(trimmed_dir, file_R2.replace('_L001_R2_001', '_R2'))

			exec_Angua.runBbduk(trimmed_input_R1, trimmer_input_R2
								trimmer_output_R1, trimmer_output_R2
								min_len, adapters, min_q)
	
	#I can probably use a function to set up the _R2 files since this repeatd the above.
	@optional_folders
	def run_trinity(mem: str, threads: int,
					trimmed_dir = "", trinity_dir = "")y
		for file_R1 in self.getFiles("trimmed", "_R1.fasta"):
			file_R2 = file.replace("_R1", "_R2")
			sample_name_R1 = getSampleName(file_R1)
			sample_name_R2 = sample_name_R1.replace("_R1", "_R2")
			
			trimmed_dir, trinity_dir = self.getFolder("trimmed"), 
								       self.getFolder("assembly")
			trinity_input_R1 = os.path.join(trimmed_dir, file)
			trinity_input_R2 = os.path.join(trimmed_dir, file)
			trinity_output = os.path.join(
							 trinity_dir, 
							 f"{sample_name_R1.replace('_R1', '_trinity')}")
			trinity_log = os.path.join(
						  trinity_dir,
						  f"{sample_name_R1.replace('_R1', '.log')}")

			log = exec_Angua.runTrinity(trinity_input_R1, trinity_input_R2,
									    trinity_output, mem, threads)
			LOG.info(log)
	
	### Sort and rename contigs
	def sort_fasta_by_length(self, min_len: int):
		out = self.extendFolder("contigs", f"sorted_{min_len}", str(min_len))
		for file in self.getFiles("contigs", ".fasta"):
			sample_name = file.split("_trinity")[0]
			output_file = os.path.join(out, 
									   f"{sample_name}_sorted_{min_length}.fasta"
			self._toolBelt.filterFasta(file, output_file, "len", min_len)
	
	@optional_folders
	def cluster_reads(self, dir_kind: str, perc, threads: int, 
							cluster_dir = "", cluster-in_dir = ""):
		if not self.getFolder(dir_kind):
			dir_kind = "cluster-in"
		if not self.getFolder("cluster"):
			out_dir = self.extendFolder(dir_kind, "cluster", "Mmseqs2")
		else:
			out_dir = self.getFolder("cluster")
		for file in self.getFiles(dir_kind):
			sample_name = getSampleName(file)
			exec_Angua.cluster(file, out_dir, perc, threads) 
			# Remove and rename intermediate files
			os.remove(os.path.join(out_dir, f"{file}_all_seqs.fasta"))
			before_rep_seq, after_rep_seq = os.path.join(
											out_dir, 
											f"{file}_rep_seq.fasta"),
											os.path.join(out_dir,
											f"{sample_name}_rep_seq.fasta")
			os.rename(before_rep_seq, after_rep_seq)
			before_tsv, after_tsv = os.path.join(out_dir,
												 f"{file}_cluster.tsv"),
												os.path.join(out_dir,
												 f"{sample_name}_cluster.tsv")
			os.rename(before_tsv, after_tsv)

	# Remove the tmp dir
	shutil.rmtree(f"{output_dir}/tmp/")

	LOG.info("Clustering complete.")
		
	def get_unmapped(min_len: int, min_alignment_len: int, threads: int,
					 trimmed_dir: ""):
	for file in self.getFiles("trimmed", ".fasta"):
		sample_name = getSampleName(file, extend=1)
		self.extendFolder("out", "unmapped", "Unmapped")
		current_bwa_dir = self.extendFolder("unmapped", "bwa_current", sample_name)
		shutil.copy(file, current_bwa_dir)
		input_file = os.path.join(current_bwa_dir, file)
		sam_file = os.path.join(current_bwa_dir, f"{sample_name}_sort.sam"})

		raw_reads = findRawBySample(sample_name)

		# Index reference
		runBwa(input_file, raw_reads, sam_file)

		# BWA and samtools
		sort_file = samToIndexedBam(in_sam = sam_file, 
									out_sam = os.path.join(current_bwa_dir, 
									f"{sample_name}.bam"))

		# Pysam
		reads_dir = self.extendFolder("unmapped", "reads_out", "Reads")
		bamfile = pysam.AlignmentFile(sort_file, "rb")
		with open(os.path.join(reads_dir, f"{sample_name}.fasta", "w") as outfile:
			for read in bamfile.fetch(until_eof = True):
				try:
					if read.query_alignment_length >= (int(read.infer_read_length()) / 3):
						if read.query_alignment_length >= int(min_alignment_len):
							good += 1
						else:
							outfile.write(f">{read.query_name}\n")
							outfile.write(f"{read.query_sequence}\n")
					else:
						outfile.write(f">{read.query_name}\n")
						outfile.write(f"{read.query_sequence}\n")
				except:
					outfile.write(f">{read.query_name}\n")
					outfile.write(f"{read.query_sequence}\n")

		# Close file
		bamfile.close()
		
def main():
	angua_version = 3
	options = parse_arguments()

	### Main Pipeline
	if(sys.argv[1] == "main"):
		angua = Angua("out", options.output)
		if not options.no_qc:
			angua.pre_qc_checks(options.qc_threads, raw_dir = options.in_dir)
		
		if options.bbduk_adapters:
			run_bbduk(options.bbduk_minl, options.bbduk_adapters, options.bbduk_q
					  raw_dir = options.in_dir)
			if not options.no_qc:
				angua.post_qc_checks(options.qc_threads)
		
		if options.assembler:
			if options.assembler.upper() == "TRINITY":
				angua.run_trinity(options.trinity_mem, options.trinity_cpu)
		
		if not options.sort:
			options.sort = [200, 1000]
		
		for num in options.sort:
			angua.sort_fasta_by_length(num)
		
		if options.cluster:
			angua.cluster_fasta(f"sorted_{max(options.sort)}", options.cluster_perc, options.cluster_threads)
		
		if options.unmapped:
			angua.get_unmapped(max(options.sort), options.bwa_threads, options.min_alignment)
			angua.cluster_fasta("reads_out", options.cluster_perc, options.cluster_threads)
			blastn_unmapped = Blast(f"{options.output}/Unmapped/Mmseqs2/", f"{options.output}/Unmapped/Blastn/", "blastn", options.blast_pool, "blastn", "megablast", options.ictv_db, options.blastn_threads, options.blast_descriptions, options.blast_alignments, "qualified", ".fasta")
			blastn_unmapped.run_blast_parallel()
			text_search(f"{options.output}/Unmapped/Blastn/", f"{options.output}/Results/", "Unmapped_reads_viruses", "", options.min_alignment)
		if options.nt_db:
			blastn = Blast(f"{options.output}/Contigs/200/", f"{options.output}/Blastn/", "blastn", options.blast_pool, "blastn", "megablast", options.nt_db, options.blastn_threads, options.blast_descriptions, options.blast_alignments, "single", ".fasta")
			blastn.run_blast_parallel()
		if options.megan_blastn == "Y":
			run_megan(f"{options.output}/Blastn/", f"{options.output}/Megan/Blastn/", "BlastN", f"{options.output}/Contigs/200/", options.megan_na2t)
		if options.blastx == "Y":
			blastx = Blast(f"{options.output}/Mmseqs2/", f"{options.output}/Blastx/", "blastx", options.blast_pool, "blastx", "blastx", options.nr_db, options.blastx_threads, options.blast_descriptions, options.blast_alignments, "qualified", ".fasta")
			blastx.run_blast_parallel()
		if options.megan_blastx == "Y":
			run_megan(f"{options.output}/Blastx/", f"{options.output}/Megan/Blastx/", "BlastX", f"{options.output}/Mmseqs2/", options.megan_pa2t)
		if options.stats == "Y":
			stats(f"{options.output}/QC/raw/MultiQC/multiqc_data/multiqc_general_stats.txt", f"{options.output}/QC/trimmed/MultiQC/multiqc_data/multiqc_general_stats.txt", ".", f"{options.output}/Trinity/", f"{options.output}/Results/")
		if options.doc_env == "Y":
			document_env("Angua", angua_version, options.output, options)

	### Plant-Finder Pipeline
	elif(sys.argv[1] == "taxa-finder"):
		
		# Directory list to create
		dirs = ["Taxa-Finder"]

		if options.create_dirs == "Y":
			create_dirs(options.output, dirs)

		if options.find_taxa == "Y":
			find_taxa(options.input, options.output, options.hmm, options.hmm_min_length)

	### Back-Mapper Pipeline
	elif(sys.argv[1] == "back-mapper"):
		
		# Directory list to create
		dirs = ["Back-Mapper"]

		if options.create_dirs == "Y":
			create_dirs(options.output, dirs)

		if options.back_mapper == "Y":
			back_mapper(options.input, f"{options.output}/Back-Mapper/", options.reference, options.threads, options.delim, options.mapq, options.flag, options.coverage)

	### Taxa-Finder Pipeline
	elif(sys.argv[1] == "text-searcher"):

		# Directory list to create
		dirs = ["Text-Searcher"]

		if options.create_dirs == "Y":
			create_dirs(options.output, dirs)

		if options.text_searcher == "Y":
			text_search(options.input, f"{options.output}/Text-Searcher/", options.output_filename ,options.search_term, options.min_alignment)
			
################################################################################
def run_megan(input_dir, output_dir, mode, reads, a2t):

	input_files = looper(input_dir, "single", ".fasta")

	for file in input_files:
		input_file = f"{input_dir}/{file}"
		input_reads = f"{reads}/{file.split('.')[0]}.fasta"

		subprocess.call(f"blast2rma -i {input_file} -f BlastXML -bm {mode} -r {input_reads} -ms 75 -sup 1 -a2t {a2t} -o {output_dir}", shell = True)

		print(f"Megan complete for sample: {file.split('.')[0]}")

################################################################################
def find_taxa(input_dir, output_dir, hmm, size):

	input_files = looper(input_dir, "qualified", ".fasta")

	for file in input_files:
		sample_name = file.split(".")[0]
		sample_in = f"{input_dir}/{file}"
		sample_out = f"{output_dir}/Taxa-Finder/{sample_name}"

		# Run nhmmer
		subprocess.call(f"nhmmer -o {sample_out}.hmm --tblout {sample_out}.tbl {hmm} {sample_in}", shell = True)

		# Parse table output
		with open(f"{sample_out}.tbl") as nhmmer_tbl:
			with open(f"{sample_out}.bed", "w") as bed_out:
				for line in nhmmer_tbl:
					if(line.startswith("#")):
						next(nhmmer_tbl)
					else:
						contig_ID = line.split()[0]
						start = line.split()[6]
						end = line.split()[7]

						if(int(line.split()[6]) > int(line.split()[7])):
							start = line.split()[7]
							end = line.split()[6]

						# Write output
						if(int(end) - int(start) > int(size)):
							bed_out.write(f"{contig_ID}	{start}	{end}\n")

			# Write fasta
			subprocess.call(f"bedtools getfasta -s -fo {sample_out}.fasta -fi {sample_in} -bed {sample_out}.bed", shell = True)

			# Move index file
			os.rename(f"{sample_in}.fai", f"{sample_out}.fasta.fai")

################################################################################
def back_mapper(input_dir, output_dir, reference, threads, delim, mapq, flag, coverage):

	with open(reference) as ref_fasta:
		seq_count = 0
		for line in ref_fasta:
			if(">" in line):
				seq_count += 1

		if(seq_count > 1):
			print("WARNING")
			print("Multiple sequences detected in the reference fasta file. Please ensure that each reference file represents ONE DISTINCT genome.")
			print("If multiple genomes are used in a single reference file, then reads that align equally well will ONLY be assigned to a single reference.")
			print("This can and will skew the mapping results.")
			print("Consider using the --flag 0 option to correct this.")


	input_files = looper(input_dir, "qualified", "_R1")

	# Index reference
	# subprocess.call(f"bwa-mem2 index {reference}", shell = True)

	for file in input_files:
		sample_ID = file.split(delim)[0]
		Path(output_dir, sample_ID).mkdir(exist_ok = True, parents = True)
		fileR2 = file.replace("_R1", "_R2")

		# Index reference
		subprocess.call(f"bwa-mem2 index {reference}", shell = True)

		if os.path.isfile(f"{input_dir}/{file}"):
			r1_in = f"{input_dir}/{file}"
			r2_in = f"{input_dir}/{fileR2}"
			sample_out = f"{output_dir}/{sample_ID}/"

			# Map raw reads
			subprocess.call(f"bwa-mem2 mem -t {threads} {reference} {r1_in} {r2_in} > {sample_out}{sample_ID}.sam", shell = True)
			subprocess.call(f"samtools view -q {mapq} -F {flag} -bS {sample_out}/{sample_ID}.sam > {sample_out}{sample_ID}.bam", shell = True)
			subprocess.call(f"samtools sort {sample_out}/{sample_ID}.bam > {sample_out}{sample_ID}_sort.bam", shell = True)
			subprocess.call(f"samtools index {sample_out}{sample_ID}_sort.bam", shell = True)
			subprocess.call(f"samtools idxstats {sample_out}{sample_ID}_sort.bam > {sample_out}{sample_ID}_stats.txt", shell = True)

			# Calculate Coverage
			if coverage == "Y":
				subprocess.call(f"average-coverage.py {sample_out}{sample_ID}_sort.bam -o {sample_out}{sample_ID}_coverage.tsv", shell = True)

	print("Back mapping completed.")

################################################################################
def stats(input_raw_stats, input_trimmed_stats, qualifier, normalised_reads_dir, output_dir):

	sample_dict = {}

	# Raw data stats
	with open(input_raw_stats) as raw_in:
		next(raw_in)
		for line in raw_in:
			line = line.strip()
			sample_id = line.split("\t")[0]
			reads = line.split("\t")[5]
			if "_R1" in sample_id:
				sample_name = sample_id.split("_L")[0]
				sample_dict[sample_name] = []
				sample_dict[sample_name].append(reads)

	# Trimmed data stats
	with open(input_trimmed_stats) as trimmed_in:
		next(trimmed_in)
		for line in trimmed_in:
			line = line.strip()
			sample_id = line.split("\t")[0]
			reads = line.split("\t")[5]
			if "_R1" in sample_id:
				sample_name = sample_id.split("_R")[0]
				sample_dict[sample_name].append(reads)

	# Normalisation stats
	input_files = looper(normalised_reads_dir, "qualified", ".log")

	for file in input_files:
		input_file = f"{normalised_reads_dir}/{file}"
		sample_name = file.split(".")[0]

		with open(input_file) as trinity_log_in:
			for line in trinity_log_in:
				if "reads selected during normalization" in line:
					line = line.strip()
					norm_reads = line.split(" ")[0]
					sample_dict[sample_name].append(norm_reads)
	
	# Output stats
	with open(f"{output_dir}/Angua_stats.tsv", "w") as stats_out:
		stats_out.write("Sample	Raw_reads	Trimmed_reads	Normalised_reads	Normalised_reads_percentage\n")
		for sample in sample_dict:
			stats_out.write(f"{sample}	{'	'.join(sample_dict[sample])}	{(float(sample_dict[sample][2]) / float(sample_dict[sample][1])) * 100}\n")

	# Output results template
	with open(f"{output_dir}/Angua_results.tsv", "w") as results_out:
		results_out.write("Sample	Blastn	Blastn_notes	Blastx	Blastx_notes\n")
		for sample in sample_dict:
			results_out.write(f"{sample}\n")

################################################################################
def document_env(script_name, script_version, output_dir, input_params):
	# Report the arguments used to run the program
	# Report the environemnt the program was run in

	print(f"Printing {script_name} pipeline version information")
	with open(f"{output_dir}/{script_name}Pipeline_params.txt", "w") as log_output:
		log_output.write(f"{script_name} Pipeline Version: {script_version}\n")
		log_output.write(f"Datetime: {datetime.datetime.now()}\n")
		log_output.write(f"Parameters:\n")
		for arg in vars(input_params):
			log_output.write(f"{arg} {getattr(input_params, arg)}\n")
	subprocess.call(f"conda list > {output_dir}/{script_name}_env.txt", shell = True)

################################################################################

################################################################################
def parse_arguments():
	parser = argparse.ArgumentParser(prog = "Angua", 
									 description = "Runs the Angua pipeline.")
	subparsers = parser.add_subparsers(help = "sub-command help")
	main = subparsers.add_parser("main", 
								 help = "Runs the Angua pipeline.")
	taxa_finder = subparsers.add_parser("taxa-finder", 
										help = "Runs the taxa-finder pipeline. Requires a HMM file and directory of fasta files. nhmmer > bedtools")
	back_mapper = subparsers.add_parser("back-mapper", 
										help = "Runs the back-mapper pipeline. bwa-mem2 > samtools > bamtocov")

	################################################################################
	### Main Pipeline

	# Key arguments
	main.add_argument("--input", 
					  help = "Path to raw data directory.", required = True)
	main.add_argument("--output", 
					  help = "Directory where output data will be generated.", required = True)
	main.add_argument("--nt_db", 
					  help = "Path to the nt database.")
	main.add_argument("--nr_db", 
					 help = "Path to the nr database.")
	main.add_argument("--ictv_db", 
					  help = "Path to the ICTV nucleotide database generated by makeICTVDB.")
	
	main.add_argument("-mn2t", "--megan_na2t", 
					  help = "Path to the megan nucl_acc2tax file.")
	main.add_argument("-mp2t", "--megan_pa2t", 
					  help = "Path to the megan prot_acc2tax file.")

	# Extra arguments, useful for if a specific job has failed and you don't want to start from scratch

	main.add_argument("-nodir", "--create_dirs_off", 
					  help = "Suppress creating directory structure.",
					  action = "store_true")
	main.add_argument("--noqc", 
					  help = "Do not run FastQC and MultiQC.",
					  action = "store_true")
	main.add_argument("-a", "--assembler", 
					  help = "Choice of assembler. 'N' skips assembly.", 
					  choices = ["trinity"]
					  default = "trinity")
	main.add_argument("-s", "--sort",
					  help = "Bins to sort contigs into. The highest will be used for Blastx, and the rest for Blastn, if these flags are set. Defaults to 200 and 1000."
					  nargs = "*",
					  type = int)
	main.add_argument("--cluster", 
					  help = "Clusters the highest bin before Blastx.",
					  action = "store_true")
	main.add_argument("-um", "--unmapped", 
					  help = "Extract unmapped and weak alignments, and looks for viruses.", 
					  action = "store_true")

	main.add_argument("-ns", "--no_stats", 
					 help = "Do not generate read stats and results template.",
					 action = "store_true")
	main.add_argument("-nde", "--no_doc_env", 
					  help = "Do not log environment and parameter details."
					  action = "store_true")

	# Tool specific parameters

	# FastQc and MultiQC
	main.add_argument("-qct", "--qc_threads", 
					  help = "Number of threads to use to generate QC plots. Default round(0.9*total threads).", 
					  default = round(os.cpu_count() * 0.9))

	# Bbduk
	main.add_argument("-bba", "--bbduk_adapters", 
					  help = "Bbduk adapter references.")
	main.add_argument("-bbq", "--bbduk_q", 
					  help = "Bbduk phred quality trim parameter. Default 10", 
					  default = 10)
	main.add_argument("-bbml", "--bbduk_minl", 
					  help = "Bbduk minimum length. Default 50", 
					  default = 50)

	# Trinity
	main.add_argument("-tcpu", "--trinity_cpu", 
					  help = "Trinity CPU parameter. Default 60.", 
					  default = 60)
	main.add_argument("-tmem", "--trinity_mem", 
					  help = "Trinity max memory parameter. Default 200G", 
					  default = "200G")

	# MMseq2
	main.add_argument("-clp", "--cluster_perc", 
					  help = "What percentage identity to cluster at. Default 0.95.", default = 0.95)
	main.add_argument("-clt", "--cluster_threads", 
					  help = "Number of threads to run mmseq2 with. Default round(0.9*total threads).", 
					  default = round(os.cpu_count() * 0.9))

	# Unmapped
	main.add_argument("-bwt", "--bwa_threads", 
					  help = "Number of threads to use with BWA. Default round(0.9*total threads).", 
					  default = round(os.cpu_count() * 0.9))
	main.add_argument("-mina", "--min_alignment", 
					  help = "Minimum alignment length for a 'good' alignment. Default 50.", 
					  default = 50)

	# Blast
	main.add_argument("-blp", "--blast_pool", 
					  help = "TMaximum number of blast processes allowed in the pool at any one time. Default 8.",
					  default = 8)
	main.add_argument("blt", "--blastn_threads", 
					  help = "Number of threads used for each blastn process. Default 16.", 
					  default = 16)
	main.add_argument("-blxt", "--blastx_threads", 
					  help = "Number of threads used for running blastx. Default 130.", 
					  default = 130)
	main.add_argument("-bld", "--blast_descriptions", 
					  help = "Number of descriptions shown. Default 25.", 
					  default = 25)
	main.add_argument("-bla", "--blast_alignments", 
					  help = "Number of alignments shown. Default 25.", 
					  default = 25)

	################################################################################
	### Taxa-finder Pipeline

	# Key arguments
	taxa_finder.add_argument("--input", 
							 help = "Location of the contig directory.", 
							 required = True)
	taxa_finder.add_argument("--output", 
							help = "Location of the output directory.", 
							 required = True)
	taxa_finder.add_argument("--hmm", 
							 help = "Location of input hmm file.",
							 required = True)

	taxa_finder.add_argument("-nodir", "--create_dirs_off", 
							 help = "Suppress creation of directory structure.", 
							 action = "store_true")
	taxa_finder.add_argument("--find_taxa", 
							 help = "Runs the specified HMM and extracts fasta files. Default true.",
							 action = "store_true")

	# Tool specific parameters

	taxa_finder.add_argument("--hmm_min_length", 
							 help = "Determine the minimum size for a sequence to be retrieved from the HMM search. Default 500.", 
							 default = 500)

	################################################################################
	### Back-Mapper Pipeline

	# Key arguments
	back_mapper.add_argument("--input", 
							help = "Location of the input directory.", 
							required = True)
	back_mapper.add_argument("--output", 
							 help = "Location of the output directory.", 
							 required = True)
	back_mapper.add_argument("--reference", 
							 help = "Location of reference fasta file.", 
							 required = True)

	back_mapper.add_argument("-nodir", "--create_dirs_off", 
							 help = "Suppress creation of directory structure.", 
							 action = "store_true")
	back_mapper.add_argument("-nocov", "--no_coverage", 
							 help = "Suppress generating coverage output.", 
							 action = "store_true")

	# Extra arguments
	back_mapper.add_argument("--threads", help = "Number of threads to use. Default is 23.", default = round(os.cpu_count() * 0.9))
	back_mapper.add_argument("--delim", help = "What chatacter separates the ID from the rest of the filename. Default is an underscore.", default = "_")
	back_mapper.add_argument("--mapq", help = "Filter reads with a MAPQ of >= X. Default is 0", default = 0)
	back_mapper.add_argument("--flag", help = "Filter reads with the specified samflag. Default is 2304", default = 2304)

	################################################################################
	return parser.parse_args()

################################################################################

if __name__ == '__main__':
	main()
