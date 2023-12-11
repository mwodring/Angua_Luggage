import os, sys, shutil, multiprocessing
import argparse
import datetime, inspect, pysam
import functools
import pandas as pd

from pathlib import Path

import dataclasses 
from dataclasses import dataclass, field

from pathlib import Path
from Bio import SeqIO
from Bio.Blast import NCBIXML

from .Luggage import fileHandler, toolBelt
from .exec_utils import *
from .utils import *
from .LuggageInterface import blastParser, spadesTidy

from collections.abc import Generator

#Doesn't use a file right this second but it will.
logging.basicConfig(stream = sys.stdout, level=logging.INFO)
LOG = logging.getLogger(__name__)

@dataclass
class Blast:
	input_files: list[str]
	output_dir: str
	mode: str
	blast_pool: int
	blast_type: str
	blast_task: str
	database_path: str
	threads: int
	alignments: int

	def run_blast_parallel(self):
		# Get pool size
		pool = multiprocessing.Pool(processes = 1)
		if self.mode == "blastn":
			pool = multiprocessing.Pool(processes = int(self.blast_pool))
		elif self.mode == "blastx":
			pool = multiprocessing.Pool(processes = 1)
		
		results = [
			pool.apply_async(
			self.blast_query, 
			args = (file, 
			f"{self.output_dir}/{os.path.basename(file)}.{self.blast_task}.xml")
			)
			for file in self.input_files
		]

		for p in results:
			p.get()

	def blast_query(self, input_file, output_file):
		LOG.info(f"Starting blast for {os.path.basename(input_file)}.")
		# Run blast query
		blast_child = subprocess.Popen(
			[self.blast_type, 
			 "-task", self.blast_task,
			 "-db", self.database_path,
			 "-query", input_file,
			 "-num_threads", str(self.threads),
			 "-outfmt", "5",
			 "-out", output_file],
			stdout = subprocess.PIPE,
			stderr = subprocess.PIPE,
			universal_newlines = True
		)
		
		#Should process the stderr before marking something as complete.
		blast_output, blast_error = blast_child.communicate()
		
		if blast_error:
			LOG.error(blast_error)
		LOG.info(f"Blast complete for query: {os.path.basename(input_file)}.")
		
class Angua(fileHandler):
	#DECORATOR
	def check_complete(func):
		@functools.wraps(func)
		def wrapper_check_complete(self, *args, **kwargs):
			stage = "_".join(func.__name__.split("_")[1:])
			finished_dir = self.getFolder(kwargs.get("out_dir", stage))
			finished_file = os.path.join(finished_dir, f"{stage}.finished")
			if os.path.exists(finished_file):
				LOG.info(f"{stage} already complete, skipping!")
			else:
				output = func(self, *args, **kwargs)
				if output:
					Path(finished_file).touch()
					LOG.info(f"{stage} completed successfully.")
				else: 
					LOG.error(f"Something went wrong with {stage}.")
		return wrapper_check_complete
		
	def extendFolderMultiple(self, orig_dir_kind: str, 
								   dir_kinds: list, dir_names: list):
		for dir_kind, dir_name in zip(dir_kinds, dir_names):
			self.extendFolder(orig_dir_kind, dir_kind, dir_name)
	
	def getFoldersIn(self, dir_kind: str):
		yield from [f.path for f in os.scandir(self.getFolder(dir_kind)) if f.is_dir()]
			
	def getContigsSorted(self, length = "min"):
		sorted_folders = [f for f in self.getFoldersIn("contigs")]
		nums = [int(os.path.basename(f)) for f in sorted_folders]
		if length == "min":
			target = min(nums)
		elif length == "max":
			target = max(nums)
		else:
			target = int(length)
		for f, n in zip(sorted_folders, nums):
			if n == target:
				return f
	
	@check_complete	
	def run_pre_qc(self, qc_threads: int) -> int:
		# Run FastQC and MultiQC on the supplied directory
		fastqc_dir = self.getFolder("QC_raw_F")
		for file in self.getFiles("raw", [".fastq.gz", ".fastq"]):
			fastQC(qc_threads, file, fastqc_dir)
		multiQC(fastqc_dir, self.getFolder("QC_raw_M"))
		return 1
	
	#Definitely repeating myself here.
	@check_complete
	def run_post_qc(self, qc_threads: int) -> int:
		fastqc_dir = self.getFolder("QC_trimmed_F")
		for file in self.getFiles("bbduk", [".fastq", ".fastq.gz"]):
			fastQC(qc_threads, file, self.getFolder("QC_trimmed_F"))
		multiQC(self.getFolder("QC_trimmed_F"),
				self.getFolder("QC_trimmed_M"))	
		return 1
		
	@check_complete
	def run_bbduk(self, min_len: int, adapters: str, min_q: int): 
		trimmed_dir = self.getFolder("bbduk")
		raw_dir = self.getFolder("raw")
		processed = 0
		for file_R1 in self.getFiles("raw", "_R1_001.fastq.gz"):
			file_R2 = file_R1.replace("_R1", "_R2")
			base_R1, base_R2 = os.path.basename(file_R1), os.path.basename(file_R2)
			sample_R1 = base_R1.replace("_L001_R1_001", "_R1")
			sample_R2 = base_R2.replace("_L001_R2_001", "_R2")
			out_R1 = os.path.join(trimmed_dir, sample_R1)
			out_R2 = os.path.join(trimmed_dir, sample_R2)
			runBbduk(file_R1, file_R2,
					 out_R1, out_R2,
					 min_len, adapters, min_q)
			processed += 1
		if processed > 0:
			return 1
		else:
			return
		
	def run_trinity(self, mem: str, threads: int):
		self.assembler = "Trinity"
		finished_file = os.path.join(self.getFolder("contigs"), 
													"Trinity.finished")
		if os.path.exists(finished_file):
			return ("Trinity already completed, skipping.")
		processed = 0
		for file_R1 in self.getFiles("bbduk", "_R1.fastq.gz"):
			file_R2 = file_R1.replace("_R1", "_R2")
			base_R1 =  os.path.basename(file_R1).split(".")[0]
			trimmed_dir, trinity_dir = self.getFolder("bbduk"), self.getFolder("contigs")
			trinity_output = os.path.join(trinity_dir, 
										  base_R1.replace("_R1", "_trinity"))
			trinity_log = trinity_output.replace("_trinity", ".log")
			LOG.info(f"Running Trinity on {getSampleName(file_R1, extend=1)}.")
			in_reads = [file_R1]
			if os.path.exists(file_R2):
				in_reads.append(file_R2)
			log = runTrinity(in_reads,
							 trinity_output, mem, threads)
			with open(trinity_log, "wb") as tlog:
				tlog.write(log)
			shutil.rmtree(trinity_output, ignore_errors=True)
			processed += 1
		if processed <= 0:
			return "Trinity failed!"
		else:
			Path(finished_file).touch()
			return "Trinity complete."
	
	def run_spades(self):
		self.assembler = "Spades"
		finished_file = os.path.join(self.getFolder("contigs"), 
													"Spades.finished")
		if os.path.exists(finished_file):
			return ("Spades already completed, skipping.")
		for file_R1 in self.getFiles("bbduk", ["_R1.fastq.gz",
											   "_R1.fastq"]):
			in_reads = [file_R1]
			file_R2 = file_R1.replace("_R1", "_R2")
			if os.path.exists(file_R2):
				in_reads.append(file_R2)
			sample_name = getSampleName(file_R1)
			out_dir = self.extendFolder("contigs", sample_name, sample_name)
			runSpades(in_reads, out_dir)
		tidy = spadesTidy("in", self.getFolder("contigs"))
		tidy.spadesToDir(self.getFolder("contigs"), cleanup = True)
		Path(finished_file).touch()
		return "Spades complete."
		
	### Sort and rename contigs - uses SeqIO directly until I can get seqtk to rename like this.
	def sort_fasta_by_length(self, min_len: int):
		output_dir = self.extendFolder("contigs", f"sorted_{min_len}", str(min_len))
		for input_file in self.getFiles("contigs", ".fasta"):
			sample_name = getSampleName(input_file, extend=1)
			output_file = os.path.join(output_dir, 
									   f"{sample_name}_sorted_{min_len}.fasta")
			with open(output_file, "w+") as contigs_out:
				for seq_record in SeqIO.parse(open(input_file, mode = "r"), "fasta"):
					if(len(seq_record.seq) >= int(min_len)):
						seq_record.id = f"{sample_name}_{seq_record.id}"
						seq_record.description = f"{sample_name}_{seq_record.description}"
						SeqIO.write(seq_record, contigs_out, "fasta")
					
	@check_complete
	def run_mmseqs2(self, perc, threads: int, in_dir = None,
											  out_dir = None) -> int:
		dir_kind = "cluster_in" if not in_dir else in_dir
		if not self.getFolder("mmseqs2"):
			out_dir = self.extendFolder(dir_kind, "mmseqs2", "Mmseqs2")
		else:
			out_dir = self.getFolder("mmseqs2")
		tmp = f"{out_dir}/tmp/"
		for file in self.getFiles(dir_kind):
			sample_name = getSampleName(file)
			out_file = os.path.join(out_dir, sample_name)
			mmseqs2(file, out_file, perc, threads, tmp) 
			# Remove and rename intermediate files
			to_remove = [file for file in self.getFiles("mmseqs2") if "all_seqs" in file]
			for file in to_remove:
				os.remove(file)
		# Remove the tmp dir
		shutil.rmtree(tmp)

		LOG.info("Clustering complete.")
		return 1
		
	@check_complete
	def run_blastn(self, options, in_dir = None, out_dir = None, ictv = False) -> int:
		in_dir = self.getContigsSorted("min") if not in_dir else self.getFolder(in_dir)
		out_dir = self.getFolder("blastn") if not out_dir else self.getFolder(out_dir)
		db = options.nt_db if not ictv else options.ictv
		blastn = Blast([os.path.join(in_dir, file) for file in os.listdir(in_dir) if file.endswith(".fasta")], 
					   out_dir, 
					   "blastn", options.blast_pool, 
					   "blastn", "megablast", 
					   db, 
					   options.blastn_threads,
					   options.blast_alignments)
		blastn.run_blast_parallel()
		return 1
			
	@check_complete
	def run_blastx(self, options, in_dir = None, out_dir = None) -> int:
		mmseqs2 = self.getFolder("mmseqs2")
		if mmseqs2 and not in_dir:
			in_files = [file for file in self.getFiles("mmseqs2", ".fasta")]
		elif not mmseqs2 and not in_dir:
			in_dir = self.getContigsSorted("max")
			in_files = [os.path.join(in_dir, file) for file in os.listdir(in_dir) 
						if file.endswith(".fasta")]
		elif in_dir:
			in_files = [file for file in self.getFolder(in_dir)]
		out_dir = self.getFolder("blastx") if not out_dir else self.getFolder(out_dir)
		blastx = Blast(in_files, out_dir, 
					   "blastx", options.blast_pool, 
					   "blastx", "blastx", 
					   options.nr_db, 
					   options.blastx_threads, 
					   options.blast_alignments)
		blastx.run_blast_parallel()
		return 1
	
	@check_complete
	def run_megan_blastn(self, na_db: str) -> int:
		self.run_megan("BlastN", na_db)
		return 1
		
	@check_complete
	def run_megan_blastx(self, pa_db: str) -> int: 
		self.run_megan("BlastX", pa_db)
		return 1
	
	def run_megan(self, blast_type: str, megan_db: str, 
						in_dir = None, out_dir = None, contigs_dir = None) -> int:
		blast_lower = blast_type.lower()
		in_dir = self.getFolder(blast_lower) if not in_dir else self.getFolder(in_dir)
		out_dir = self.getFolder(f"megan_{blast_lower}") if not out_dir else self.getFolder(out_dir)
		if not contigs_dir:
			contigs = self.getContigsSorted("min") if blast_type == "BlastN" else self.getFolder("mmseqs2")
			if not contigs:
				contigs = self.getContigsSorted("max")
		else:
			contigs = self.getFolder(in_dir)
		self.addFolder("megan_in", contigs)
		self.findFastaFiles("megan_in")
		files = (os.path.join(in_dir, file) for file in os.listdir(in_dir) if file.endswith(".xml"))
		for file in files:
			sample_name = getSampleName(file)
			current_contigs = self._toolBelt.getToolsByName("fasta", sample_name)[0].filename
			self._toolBelt.blast2Rma(file, out_dir, megan_db, 
									 current_contigs, blast_type,
									 sample_name)
		return 1

	def stats(self, options):
		def renameReads(dataframe, new_name: str):
			dataframe = dataframe[is_R1_mask]
			dataframe = dataframe[["Sample", 
								   "FastQC_mqc-generalstats-fastqc-percent_fails"]]
			dataframe = dataframe.rename(columns = 
				      {"FastQC_mqc-generalstats-fastqc-percent_fails": 
				       new_name})
			return dataframe
		# Raw data stats
		input_raw_stats = os.path.join(self.getFolder("QC_raw_M"),
									   "multiqc_data",
									   "multiqc_general_stats.txt")
		df = pd.read_csv(input_raw_stats, sep = "\t")
	
		is_R1_mask = df['Sample'].str.contains('R1') 
		df = renameReads(df, "raw_reads_failed")
		cut_samples = {sample : 
					  sample.replace("_L001_R1_001", "_R1") for sample in df["Sample"]}
		df["Sample"].replace(cut_samples, inplace = True)
		
		if options.bbduk_adapters:
			# Trimmed data stats
			input_trimmed_stats = os.path.join(self.getFolder("QC_trimmed_M"),
											   "multiqc_data",
											   "multiqc_general_stats.txt") 
			t_df = pd.read_csv(input_trimmed_stats, sep = "\t")
			t_df = renameReads(t_df, "trimmed_reads_failed")
			df = pd.merge(df, t_df, on="Sample")

			assem_data = []
		
		#Trinity stats:
		if self.assembler == "Trinity":
			for sample_name, norm_reads in self.getTrinityNormReads():
					assem_data.append({"Sample": f"{sample_name}_R1",
									   "normalised_reads" : norm_reads})
			a_df = pd.DataFrame(assem_data)
			df = pd.merge(df, a_df, on="Sample")

		# Output stats
		out_dir = self.getFolder("results")
		out_stats = os.path.join(out_dir, "Angua_stats.tsv")
		with open(out_stats, "w") as stats_out:
			df.to_csv(stats_out, sep = "\t")
				
	def getTrinityNormReads(self):
			for file in self.getFiles("contigs", ".log"):
				sample_name = getSampleName(file)
				with open(file, "r") as trinity_log_in:
					for line in trinity_log_in:
						if "reads selected during normalization" in line:
							line = line.strip()
							norm_reads = line.split(" ")[0]
							yield sample_name, norm_reads
	
	def document_env(self, script_name: str, script_version: float, input_params):
		# Report the arguments used to run the program
		# Report the environemnt the program was run in

		LOG.info(f"Printing {script_name} pipeline version information")
		out_log = os.path.join(self.getFolder("out"),
								f"{script_name}Pipeline_params.txt")
		env_out = os.path.join(self.getFolder("out"),
							   f"{script_name}_env.txt")
		with open(out_log, "w") as log_output:
			log_output.write(f"{script_name} Pipeline Version: {script_version}\n")
			log_output.write(f"Datetime: {datetime.datetime.now()}\n")
			log_output.write(f"Parameters:\n")
			for arg in vars(input_params):
				log_output.write(f"{arg} {getattr(input_params, arg)}\n")
		with open(env_out, "w") as txt:
			subprocess.run(["conda", "list"], stdout = txt)

class taxaFinder(fileHandler):
	def find_taxa(self, hmm, size):
		for file in self.getFiles("in", ".fasta"):
			sample_name = getSampleName(file)
			sample_out = os.path.join(self.getFolder("out"),
									  sample_name)

			# Run nhmmer
			runNhmmer(sample_out, hmm, file)
			
			# Parse table output
			with open(f"{sample_out}.tbl") as nhmmer_tbl:
				with open(f"{sample_out}.bed", "w+") as bed_out:
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
							if(int(end) - int(start) > int(size)):
								bed_out.write(f"{contig_ID}	{start}	{end}\n")
					bedtoolsWriteFasta(sample_out, file)
				os.rename(f"{file}.fai", f"{sample_out}.fasta.fai")

class backMapper(blastParser):
	def runBwaTS(self, threads, mapq, flag):
		all_R1 = (file for file in os.listdir(self.getFolder("trimmed"))
                       if "R1" in file)
		for R1 in all_R1:	
			sample_name = getSampleName(R1, extend=1)
			both_reads = self.findFastaBySample(sample_name, dir_kind = "trimmed")
			tmp_dir = self.extendFolder("trimmed", "tmp_dir", "tmp")
			tmp_out = [os.path.join(tmp_dir, os.path.basename(read)) for read in both_reads]
			for r_in, r_out in zip(both_reads, tmp_out):
				shutil.copyfile(r_in, r_out)
			super().runBwaTS(tmp_dir, "ref", 0, threads, mapq, flag, text_search = False)
			shutil.rmtree(tmp_dir)

def main():
	angua_version = 4
	options = parse_arguments()

	### Main Pipeline
	if(sys.argv[1] == "main"):
		angua = Angua("out", os.path.abspath(options.output))
		angua.addFolder("raw", os.path.abspath(options.input))
		if not options.noqc:
			angua.extendFolder("out", "QC", "QC")
			angua.extendFolderMultiple("QC", ["pre_qc", "post_qc"],
											 ["Raw", "Trimmed"])
			angua.extendFolderMultiple("pre_qc", ["QC_raw_F", "QC_raw_M"],
												 ["FastQC", "MultiQC"])
			angua.run_pre_qc(options.qc_threads)
		
		if options.bbduk_adapters:
			angua.extendFolder("out", "bbduk", "Bbduk")
			angua.run_bbduk(options.bbduk_minl, options.bbduk_adapters, 
							options.bbduk_q)
			if not options.noqc:
				angua.extendFolderMultiple("post_qc", ["QC_trimmed_F", "QC_trimmed_M"],
													  ["FastQC", "MultiQC"])
				angua.run_post_qc(options.qc_threads)
		else:
			angua.addFolder("bbduk", os.path.abspath(options.input))
		
		if options.assembler != "N":
			if options.assembler == "trinity":
				angua.extendFolder("out", "contigs", "Trinity")
				LOG.info(angua.run_trinity(options.trinity_mem, 
										   options.trinity_cpu))
			if options.assembler == "spades":
				angua.extendFolder("out", "contigs", "Spades")
				LOG.info(angua.run_spades())
		
		if not options.sort:
			options.sort = [200, 1000]
		
		for num in options.sort:
			angua.sort_fasta_by_length(num)
		
		if options.cluster:
			angua.extendFolder("out", "mmseqs2", "mmseqs2")
			angua.run_mmseqs2(options.cluster_perc, 
							  options.cluster_threads,
							  f"sorted_{max(options.sort)}")
		if options.nt_db:
			out_dir = angua.extendFolder("out", "blastn", "BlastN")
			angua.run_blastn(options)
		
		if options.megan_na2t:
			angua.extendFolder("out", "megan", "Megan")
			blastn = angua.extendFolder("megan", "megan_blastn", "BlastN")
			angua.run_megan_blastn(options.megan_na2t)
							
		if options.nr_db:
			out_dir = angua.extendFolder("out", "blastx", "BlastX")
			angua.run_blastx(options)
			
		if options.megan_pa2t:
			angua.extendFolder("out", "megan", "Megan")
			blastx = angua.extendFolder("megan", "megan_blastx", "BlastX")
			angua.run_megan_blastx(options.megan_pa2t)
			
		if not options.no_stats:
			angua.extendFolder("out", "results", "Results")
			angua.stats(options)
		if not options.no_doc_env:
			angua.document_env("Angua", angua_version, options)

	### Plant-Finder Pipeline
	elif(sys.argv[1] == "taxa-finder"):
		
		tf = taxaFinder("in", os.path.abspath(options.input))
		tf.addFolder("out", os.path.abspath(options.output))
		tf.find_taxa(options.hmm, options.hmm_min_length)

	### Back-Mapper Pipeline
	elif(sys.argv[1] == "back-mapper"):
		#Atm this hooks into some functions blastParser has. Which feels hacky but it works.
		bm = backMapper("out", os.path.abspath(options.output))
		bm.addFolder("ref", os.path.abspath(options.input))
		bm.addFolder("trimmed", os.path.abspath(options.trimmed))
		bm.runBwaTS(options.threads, options.mapq, options.flag)
						  
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
										help = "Runs the back-mapper pipeline.")

	################################################################################
	### Main Pipeline

	# Key arguments
	main.add_argument("input", 
					  help = "Path to raw data directory.")
	main.add_argument("output", 
					  help = "Directory where output data will be generated.")
	main.add_argument("--nt_db", 
					  help = "Path to the nt database.")
	main.add_argument("--nr_db", 
					 help = "Path to the nr database.")

	main.add_argument("-na2t", "--megan_na2t", 
					  help = "Path to the megan nucl_acc2tax file.")
	main.add_argument("-pa2t", "--megan_pa2t", 
					  help = "Path to the megan prot_acc2tax file.")

	# Extra arguments, useful for if a specific job has failed and you don't want to start from scratch

	main.add_argument("--noqc", 
					  help = "Do not run FastQC and MultiQC.",
					  action = "store_true")
	main.add_argument("-a", "--assembler", 
					  help = "Choice of assembler. 'N' skips assembly.", 
					  choices = ["trinity", "spades"],
					  default = "trinity")
	main.add_argument("-s", "--sort",
					  help = "Bins to sort contigs into. The highest will be used for Blastx, and the rest for Blastn, if these flags are set. Defaults to 200 and 1000.",
					  nargs = "*",
					  type = int)
	main.add_argument("--cluster", 
					  help = "Clusters the highest bin before Blastx.",
					  action = "store_true")

	main.add_argument("-ns", "--no_stats", 
					 help = "Do not generate read stats and results template.",
					 action = "store_true")
	main.add_argument("-nde", "--no_doc_env", 
					  help = "Do not log environment and parameter details.",
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
	main.add_argument("-blt", "--blastn_threads", 
					  help = "Number of threads used for each blastn process. Default 16.", 
					  default = 16)
	main.add_argument("-blxt", "--blastx_threads", 
					  help = "Number of threads used for running blastx. Default 130.", 
					  default = 130)
	main.add_argument("-bla", "--blast_alignments", 
					  help = "Number of alignments shown. Default 25.", 
					  default = 25)

	################################################################################
	### Taxa-finder Pipeline

	# Key arguments
	taxa_finder.add_argument("input", 
							 help = "Location of the contig directory.")
	taxa_finder.add_argument("output", 
							help = "Location of the output directory.")
	taxa_finder.add_argument("--hmm", 
							 help = "Location of input hmm file.",
							 required = True)
							 
	# Tool specific parameters

	taxa_finder.add_argument("-hmmml", "--hmm_min_length", 
							 help = "Determine the minimum size for a sequence to be retrieved from the HMM search. Default 500.", 
							 default = 500)

	################################################################################
	### Back-Mapper Pipeline

	# Key arguments
	back_mapper.add_argument("input", 
							help = "Location of the input directory of reference files.")
	back_mapper.add_argument("trimmed", 
							 help = "Location of trimmed reads to map.")
	back_mapper.add_argument("output", 
							 help = "Location of the output directory.")

	# Extra arguments
	back_mapper.add_argument("-t", "--threads", 
							 help = "Number of threads to use. Default is 23.", 
							 default = round(os.cpu_count() * 0.9))
	back_mapper.add_argument("-mq", "--mapq", 
							 help = "Filter reads with a MAPQ of >= X. Default is 0", 
							 default = 0)
	back_mapper.add_argument("-f", "--flag", 
							 help = "Filter reads with the specified samflag. Default is 2304", 
							 default = 2304)

	return parser.parse_args()

################################################################################

if __name__ == '__main__':
	sys.exit(main())