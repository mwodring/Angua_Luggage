import os
import sys
import subprocess
import multiprocessing
import argparse
import shutil
import datetime
import inspect
import pysam

from pathlib import Path
from Bio import SeqIO
from Bio.Blast import NCBIXML

def main():

	# Angua script version
	angua_version = 1

	# Input arguments
	options = parse_arguments()

	### Main Pipeline
	if(sys.argv[1] == "main"):

		# Directory list to create
		dirs = ["Results",
			"QC",
			"QC/raw/FastQC",
			"QC/raw/MultiQC",
			"QC/trimmed/FastQC",
			"QC/trimmed/MultiQC",
			"Bbduk",
			"Trinity",
			"Contigs",
			"Contigs/200",
			"Contigs/1000",
			"Mmseqs2",
			"Unmapped",
			"Unmapped/BWA",
			"Unmapped/Reads",
			"Unmapped/Mmseqs2",
			"Unmapped/Blastn",
			"Blastn",
			"Blastx",
			"Megan",
			"Megan/Blastn",
			"Megan/Blastx",
			"Results"
			]

		if options.create_dirs == "Y":
			create_dirs(options.output, dirs)
		if options.qc == "Y":
			qc_checks(options.input, f"{options.output}/QC/raw/", options.qc_threads)
		if options.qc == "Y":
			run_bbduk(options.input, f"{options.output}/Bbduk/", options.bbduk_minl, options.bbduk_adapters, options.bbduk_q)
		if options.qc == "Y":
			qc_checks(f"{options.output}/Bbduk/", f"{options.output}/QC/trimmed/", options.qc_threads)
		if options.trinity == "Y":
			run_trinity(f"{options.output}/Bbduk/", f"{options.output}/Trinity/", options.trinity_mem, options.trinity_cpu)
		if options.sort == "Y":
			sort_fasta(f"{options.output}/Trinity/", f"{options.output}/Contigs/200/", 200)
			sort_fasta(f"{options.output}/Trinity/", f"{options.output}/Contigs/1000/", 1000)
		if options.cluster == "Y":
			cluster_fasta(f"{options.output}/Contigs/1000/", f"{options.output}/Mmseqs2/", options.cluster_perc, options.cluster_threads)
		if options.unmapped == "Y":
			get_unmapped(f"{options.output}/Contigs/200/", f"{options.output}/Bbduk/", f"{options.output}/Unmapped/", options.bwa_threads, options.min_alignment)
			cluster_fasta(f"{options.output}/Unmapped/Reads/", f"{options.output}/Unmapped/Mmseqs2/", options.cluster_perc, options.cluster_threads)
			blastn_unmapped = Blast(f"{options.output}/Unmapped/Mmseqs2/", f"{options.output}/Unmapped/Blastn/", "blastn", options.blast_pool, "blastn", "megablast", options.ictv_db, options.blastn_threads, options.blast_descriptions, options.blast_alignments, "qualified", ".fasta")
			blastn_unmapped.run_blast_parallel()
			text_search(f"{options.output}/Unmapped/Blastn/", f"{options.output}/Results/", "Unmapped_reads_viruses", "", options.min_alignment)
		if options.blastn == "Y":
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
def create_dirs(output_dir, dirs):
	# Create the directories if they do not exist

	for directory in dirs:
		Path(output_dir, directory).mkdir(exist_ok = True, parents = True)

	print("Directory creation completed.")

################################################################################
def qc_checks(input_dir, output_dir, threads):
	# Run FastQC and MultiQC on the supplied directory

	subprocess.call(f"fastqc -t {threads} {input_dir}/* -o {output_dir}/FastQC/", shell = True)
	subprocess.call	(f"multiqc {output_dir}/FastQC/ -o {output_dir}/MultiQC/", shell = True)

	print("Quality control checks completed.")

################################################################################
def run_bbduk(input_dir, output_dir, min_len, adapters, min_q):
	# Call the looper function to get a list of files in the directory
	# Run Bbduk on the returned list of flies

	input_files = looper(input_dir, "qualified", "_R1")

	for file in input_files:
		file_R2 = file.replace("_R1", "_R2")
		sample_name_R1 = file.split(".")[0]
		sample_name_R2 = sample_name_R1.replace("_R1", "_R2")

		trimmer_input_R1 = f"{input_dir}/{file}"
		trimmer_input_R2 = f"{input_dir}/{file_R2}"
		trimmer_output_R1 = f"{output_dir}/{file.replace('_L001_R1_001', '_R1')}"
		trimmer_output_R2 =f"{output_dir}/{file_R2.replace('_L001_R2_001', '_R2')}"

		subprocess.call(f"bbduk.sh in1={trimmer_input_R1} in2={trimmer_input_R2} out1={trimmer_output_R1} out2={trimmer_output_R2} minlen={min_len} ktrim=r k=23 mink=11 hdist=1 ref={adapters} qtrim=r trimq={min_q}", shell = True)

	print("Bbduk complete.")

################################################################################
def run_trinity(input_dir, output_dir, mem, threads):
	# Call the looper function to get a list of files in the directory
	# Run Trinity on the returned list of flies

	input_files = looper(input_dir, "qualified", "_R1")

	for file in input_files:
		file_R2 = file.replace("_R1", "_R2")
		sample_name_R1 = file.split(".")[0]
		sample_name_R2 = sample_name_R1.replace("_R1", "_R2")

		trinity_input_R1 = f"{input_dir}/{file}"
		trinity_input_R2 = f"{input_dir}/{file_R2}"
		trinity_output = f"{output_dir}/{sample_name_R1.replace('_R1', '_trinity')}"
		trinity_log = f"{output_dir}/{sample_name_R1.replace('_R1', '.log')}"

		subprocess.call(f"Trinity --seqType fq --max_memory {mem} --left {trinity_input_R1} --right {trinity_input_R2} --CPU {threads} --full_cleanup --output {trinity_output} > {trinity_log}", shell = True)

	print("Trinity complete.")

################################################################################
def sort_fasta(input_dir, output_dir, min_length):
	### Sort and rename contigs

	input_files = looper(input_dir, "single", ".fasta")

	for file in input_files:
		if(file.endswith(".fasta")):
			sample_name = file.split("_trinity")[0]
			input_file = f"{input_dir}/{file}"
			output_file = f"{output_dir}/{sample_name}_sorted_{min_length}.fasta"

			with open(output_file, "w") as contigs_out:
				for seq_record in SeqIO.parse(open(input_file, mode = "r"), "fasta"):
					seq_record.id = f"{output_dir.split('/')[0]}_{sample_name.split('_')[-1]}_{seq_record.id}"
					seq_record.description = f"{output_dir.split('/')[0]}_{sample_name.split('_')[-1]}_{seq_record.description}"

					if(len(seq_record.seq) >= int(min_length)):
						SeqIO.write(seq_record, contigs_out, "fasta")

	print("Sorting complete.")

################################################################################
def cluster_fasta(input_dir, output_dir, cluster_perc, threads):

	input_files = looper(input_dir, "single", ".fasta")

	for file in input_files:
		sample_name = file.split(".fasta")[0]

		subprocess.call(f"mmseqs easy-cluster -c {cluster_perc} --threads {threads} -v 0 {input_dir}/{file} {output_dir}/{sample_name}.fasta {output_dir}/tmp", shell = True)

		# Remove and rename intermediate files
		os.remove(f"{output_dir}/{file}_all_seqs.fasta")
		os.rename(f"{output_dir}/{file}_rep_seq.fasta", f"{output_dir}/{sample_name}_rep_seq.fasta")
		os.rename(f"{output_dir}/{file}_cluster.tsv", f"{output_dir}/{sample_name}_cluster.tsv")

	# Remove the tmp dir
	shutil.rmtree(f"{output_dir}/tmp/")

	print("Clustering complete.")

################################################################################
def get_unmapped(input_assembly_dir, input_reads_dir, output_dir, min_alignment_len, threads):
	
	input_files = looper(input_assembly_dir, "single", ".fasta")

	for file in input_files:
		sample_name = file.split("_sorted")[0]

		Path(f"{output_dir}/BWA/", sample_name).mkdir(exist_ok = True, parents = True)
		subprocess.call(f"cp {input_assembly_dir}/{file} {output_dir}/BWA/{sample_name}/{file}", shell = True)
		input_file = f"{output_dir}/BWA/{sample_name}/{file}"

		r1_in = f"{input_reads_dir}/{sample_name}_R1.fastq.gz"
		r2_in = f"{input_reads_dir}/{sample_name}_R2.fastq.gz"

		# Index reference
		subprocess.call(f"bwa-mem2 index {input_file}", shell = True)

		# BWA and samtools
		subprocess.call(f"bwa-mem2 mem -t {threads} {input_file} {r1_in} {r2_in} > {output_dir}/BWA/{sample_name}/{sample_name}.sam", shell = True)
		subprocess.call(f"samtools view -q 0 -F 2304 -bS {output_dir}/BWA/{sample_name}/{sample_name}.sam > {output_dir}/BWA/{sample_name}/{sample_name}.bam", shell = True)
		subprocess.call(f"samtools sort {output_dir}/BWA/{sample_name}/{sample_name}.bam > {output_dir}/BWA/{sample_name}/{sample_name}_sort.bam", shell = True)
		subprocess.call(f"samtools index {output_dir}/BWA/{sample_name}/{sample_name}_sort.bam", shell = True)

		# Pysam
		bamfile = pysam.AlignmentFile(f"{output_dir}/BWA/{sample_name}/{sample_name}_sort.bam", "rb")
		with open(f"{output_dir}/Reads/{sample_name}.fasta", "w") as outfile:
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

################################################################################
def run_megan(input_dir, output_dir, mode, reads, a2t):

	input_files = looper(input_dir, "single", ".fasta")

	for file in input_files:
		input_file = f"{input_dir}/{file}"
		input_reads = f"{reads}/{file.split('.')[0]}.fasta"

		subprocess.call(f"blast2rma -i {input_file} -f BlastXML -bm {mode} -r {input_reads} -ms 75 -sup 1 -a2t {a2t} -o {output_dir}", shell = True)

		print(f"Megan complete for sample: {file.split('.')[0]}")

################################################################################
def looper(input_dir, mode, qualifier):
	# Loop through files in a directory and return a list of contents

	files = []
	for file in os.listdir(input_dir):
		if mode == "qualified":
			if(qualifier in file):
				files.append(file)
		elif mode == "single":
			files.append(file)

	return files

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
def text_search(input_dir, output_dir, output_filename, search_term, min_alignment_len):

	input_files = looper(input_dir, "qualified", ".xml")

	all_virus_dict = {}
	for file in input_files:
		input_file = f"{input_dir}/{file}"
		sample_ID = file.split(".")[0]

		with open(input_file) as xml_in:
			blast_records = NCBIXML.parse(xml_in)

			virus_dict = {}
			for record in blast_records:
				for alignment in record.alignments:
					if search_term in alignment.hit_def.upper():

						# Check alignment length is not shorter than 1/3 of the query length
						if len(alignment.hsps[0].match) >= (record.query_length / 3):
							# Check the alignment length is not shorter than 50nts
							if len(alignment.hsps[0].match) >= int(min_alignment_len):

								# Assign variables
								hit = alignment.hit_def
								pid = alignment.hsps[0].identities / len(alignment.hsps[0].match) * 100
								bitscore = alignment.hsps[0].bits
								coverage = len(alignment.hsps[0].match) / record.query_length * 100
								alignment_len = len(alignment.hsps[0].match)
								query_len = record.query_length
								aligned_seq = alignment.hsps[0].query

								blast_vars = (hit, pid, bitscore, coverage, alignment_len, query_len, aligned_seq)

								# Add hit to dict
								if record.query in virus_dict:
									if alignment.hsps[0].bits >= list(virus_dict[record.query].items())[0][1][2] * 0.9:
										virus_dict[record.query][alignment.hit_id] = []

										# Append to dict
										for i in blast_vars:
											virus_dict[record.query][alignment.hit_id].append(i)

								else:
									virus_dict[record.query] = {}
									virus_dict[record.query][alignment.hit_id] = []

									# Append to dict
									for i in blast_vars:
										virus_dict[record.query][alignment.hit_id].append(i)

			all_virus_dict[sample_ID] = virus_dict

	# Write output
	with open(f"{output_dir}/{output_filename}.tsv", "w") as unmapped_out:
		unmapped_out.write("Sample	Read_ID	Database_hit	Percentage_identity	Coverage	Query_length	Aligned_sequence\n")
		for sample in all_virus_dict:
			for read in all_virus_dict[sample]:
				for hit in all_virus_dict[sample][read]:
					taxa = all_virus_dict[sample][read][hit][0]
					pid = all_virus_dict[sample][read][hit][1]
					cov = all_virus_dict[sample][read][hit][3]
					query_len = all_virus_dict[sample][read][hit][5]
					aligned_seq = all_virus_dict[sample][read][hit][6]

					unmapped_out.write(f"{sample}	{read}	{taxa}	{pid}	{cov}	{query_len}	{aligned_seq}\n")

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
class Blast:
	def __init__(
		self,
		input_dir,
		output_dir,
		mode,
		blast_pool,
		blast_type,
		blast_task,
		database_path,
		threads,
		descriptions,
		alignments,
		loop_mode,
		loop_qualifier
	):
		self.input_dir = input_dir
		self.output_dir = output_dir
		self.mode = mode
		self.blast_pool = blast_pool
		self.blast_type = blast_type
		self.blast_task = blast_task
		self.database_path = database_path
		self.threads = threads
		self.descriptions = descriptions
		self.alignments = alignments
		self.loop_mode = loop_mode
		self.loop_qualifier = loop_qualifier

	############################################################################
	def run_blast_parallel(self):

		input_files = looper(self.input_dir, self.loop_mode, self.loop_qualifier)

		# Get pool size
		pool = multiprocessing.Pool(processes = 1)
		if self.mode == "blastn":
			pool = multiprocessing.Pool(processes = int(self.blast_pool))
		elif self.mode == "blastx":
			pool = multiprocessing.Pool(processes = 1)

		results = [
			pool.apply_async(self.blast_query, args = (f"{self.input_dir}/{file}", f"{self.output_dir}/{file.split('.')[0]}.{self.blast_task}.xml")
			)
			for file in input_files
		]

		for p in results:
			p.get()

	############################################################################
	def blast_query(self, input_file, output_file):

		# Blast query
		blast = (
			f"{self.blast_type} "
			f"-task {self.blast_task} "
			f"-db {self.database_path} "
			f"-query {input_file} "
			f"-num_threads {self.threads} "
			f"-num_descriptions {self.descriptions} "
			f"-num_alignments {self.alignments} "
			f"-outfmt 5 "
			f"-out {output_file}"
		)

		# Run blast query
		blast_child = subprocess.Popen(
			str(blast),
			stdout = subprocess.PIPE,
			stderr = subprocess.PIPE,
			universal_newlines = True,
			shell = (sys.platform != "win32")
		)

		blast_output, blast_error = blast_child.communicate()

		print(f'Blast complete for query: {input_file.split("/")[-1]}')

################################################################################
def parse_arguments():
	parser = argparse.ArgumentParser(prog = "Angua", description = "Runs the Angua pipeline.")
	subparsers = parser.add_subparsers(help = "sub-command help")
	main = subparsers.add_parser("main", help = "Runs the Angua pipeline.")
	taxa_finder = subparsers.add_parser("taxa-finder", help = "Runs the taxa-finder pipeline. Requires a HMM file and directory of fasta files. nhmmer > bedtools")
	back_mapper = subparsers.add_parser("back-mapper", help = "Runs the back-mapper pipeline. bwa-mem2 > samtools > bamtocov")
	text_searcher = subparsers.add_parser("text-searcher", help = "Runs the taxa-finder pipeline. Searches xml files for search terms.")

	################################################################################
	### Main Pipeline

	# Key arguments
	main.add_argument("--input", help = "This is the location of the raw data directory.", required = True)
	main.add_argument("--output", help = "This is where the output data will be generated.", required = True)

	main.add_argument("--nt_db", help = "This is the path to the nt database.", default = "/data/bigbio_00/smcgreig/Blast_databases/nt/nt_db_28012022/nt")
	main.add_argument("--nr_db", help = "This is the path to the nr database.", default = "/data/bigbio_00/smcgreig/Blast_databases/nr/nr_db_27012022/nr")
	main.add_argument("--ictv_db", help = "This is the path to the ictv database.", default = "/data/bigbio_00/smcgreig/Blast_databases/ICTV/ICTV_viruses")
	main.add_argument("--megan_na2t", help = "This is the path to the megan nucl_acc2tax file.", default = "/home/smcgreig/miniconda2/envs/Angua3/opt/megan-6.12.3/megan-nucl-Jan201.db")
	main.add_argument("--megan_pa2t", help = "This is the path to the megan prot_acc2tax file.", default = "/home/smcgreig/miniconda2/envs/Angua3/opt/megan-6.12.3/megan-map-Jan2021.db")

	# Extra arguments, useful for if a specific job has failed and you don't want to start from scratch

	main.add_argument("--create_dirs", help = "Creates the directory structure. Default Y.", choices = ["Y", "N"], default = "Y")
	main.add_argument("--bbduk", help = "Runs bbduk. Default Y.", choices = ["Y", "N"], default = "Y")
	main.add_argument("--qc", help = "Runs fastQC and multiQC before and after trimming. Default Y.", choices = ["Y", "N"], default = "Y")
	main.add_argument("--trinity", help = "Run trinity. Default Y.", choices = ["Y", "N"], default = "Y")
	main.add_argument("--sort", help = "Sort contigs, based on length, into >=200 and >=1000. Default Y.", choices = ["Y", "N"], default = "Y")
	main.add_argument("--cluster", help = "Clusters contigs >= 1000. Default Y.", choices = ["Y", "N"], default = "Y")
	main.add_argument("--unmapped", help = "Extracts unmapped and weak alignments, and looks for viruses", choices = ["Y", "N"], default = "Y")
	main.add_argument("--blastn", help = "Run blastn. Default Y.", choices = ["Y", "N"], default = "Y")
	main.add_argument("--blastx", help = "Run blastx. Default Y.", choices = ["Y", "N"], default = "Y")
	main.add_argument("--megan_blastn", help = "Run megan for blastn. Default Y.", choices = ["Y", "N"], default = "Y")
	main.add_argument("--megan_blastx", help = "Run megan for blastx. Default Y.", choices = ["Y", "N"], default = "Y")
	main.add_argument("--stats", help = "Generate read stats and results template. Default Y.", choices = ["Y", "N"], default = "Y")
	main.add_argument("--doc_env", help = "Creates environment and parameter details. Default Y.", choices = ["Y", "N"], default = "Y")

	# Tool specific parameters

	# FastQc and MultiQC
	main.add_argument("--qc_threads", help = "The number of threads to use to generate QC plots. Default round(0.9*total threads).", default = round(os.cpu_count() * 0.9))

	# Bbduk
	main.add_argument("--bbduk_adapters", help = "Bbduk adapter references.", default = "/home/smcgreig/Scripts/Angua_test/RNA_adapters.fasta")
	main.add_argument("--bbduk_q", help = "Bbduk phred quality trim parameter. Default 10", default = 10)
	main.add_argument("--bbduk_minl", help = "Bbduk minimum length. Default 50", default = 50)

	# Trinity
	main.add_argument("--trinity_cpu", help = "Trinity CPU parameter. Default 60.", default = 60)
	main.add_argument("--trinity_mem", help = "Trinity max memory parameter. Default 200G", default = "200G")

	# MMseq2
	main.add_argument("--cluster_perc", help = "What percentage identity to cluster at. Default 0.95.", default = 0.95)
	main.add_argument("--cluster_threads", help = "Number of threads to run mmseq2 with. Default round(0.9*total threads).", default = round(os.cpu_count() * 0.9))

	# Unmapped
	main.add_argument("--bwa_threads", help = "The number of threads to use with BWA. Default round(0.9*total threads).", default = round(os.cpu_count() * 0.9))
	main.add_argument("--min_alignment", help = "The minimum alignment length for a 'good' alignment. Default 50.", default = 50)

	# Blast
	main.add_argument("--blast_pool", help = "This is the maximum number of blast processes allowed in the pool at any one time. Default 8.", default = 8)
	main.add_argument("--blastn_threads", help = "This is the number of threads used for each blastn process. Default 16.", default = 16)
	main.add_argument("--blastx_threads", help = "This is the number of threads used for running blastx. Default 130.", default = 130)
	main.add_argument("--blast_descriptions", help = "This is the number of descriptions shown. Default 25.", default = 25)
	main.add_argument("--blast_alignments", help = "This is the number of alignments shown. Default 25.", default = 25)

	################################################################################
	### Taxa-finder Pipeline

	# Key arguments
	taxa_finder.add_argument("--input", help = "This is the location of the contig directory.", required = True)
	taxa_finder.add_argument("--output", help = "This is the location of the output directory.", required = True)
	taxa_finder.add_argument("--hmm", help = "This is the location of input hmm file. Default /biostore/bigbio_00/smcgreig/rbcL_databases/rbcl_hmm/rbcl.hmm", default = "/biostore/bigbio_00/smcgreig/rbcL_databases/rbcl_hmm/rbcl.hmm")

	taxa_finder.add_argument("--create_dirs", help = "Creates the directory structure. Default Y.", choices = ["Y", "N"], default = "Y")
	taxa_finder.add_argument("--find_taxa", help = "Runs the specified HMM and extracts fasta files. Default Y.", choices = ["Y", "N"], default = "Y")

	# Tool specific parameters

	taxa_finder.add_argument("--hmm_min_length", help = "Determine the minimum size for a sequence to be retrieved from the HMM search. Default 500.", default = 500)

	################################################################################
	### Back-Mapper Pipeline

	# Key arguments
	back_mapper.add_argument("--input", help = "This is the location of the input directory.", required = True)
	back_mapper.add_argument("--output", help = "This is the location of the output directory.", required = True)
	back_mapper.add_argument("--reference", help = "This is the location of reference fasta file.", required = True)

	back_mapper.add_argument("--create_dirs", help = "Creates the directory structure. Default Y.", choices = ["Y", "N"], default = "Y")
	back_mapper.add_argument("--back_mapper", help = "Runs bwa and samtools. Default Y.", choices = ["Y", "N"], default = "Y")
	back_mapper.add_argument("--coverage", help = "Generate coverage output. Default Y.", choices = ["Y", "N"], default = "Y")

	# Extra arguments
	back_mapper.add_argument("--threads", help = "Number of threads to use. Default is 23.", default = round(os.cpu_count() * 0.9))
	back_mapper.add_argument("--delim", help = "What chatacter separates the ID from the rest of the filename. Default is an underscore.", default = "_")
	back_mapper.add_argument("--mapq", help = "Filter reads with a MAPQ of >= X. Default is 0", default = 0)
	back_mapper.add_argument("--flag", help = "Filter reads with the specified samflag. Default is 2304", default = 2304)

	################################################################################
	### Text-Search Pipeline

	# Key arguments
	text_searcher.add_argument("--input", help = "This is the location of the input directory.", required = True)
	text_searcher.add_argument("--output", help = "This is the location of the output directory.", required = True)
	text_searcher.add_argument("--search_term", help = "Term to search for. Default VIRUS.", default = "VIRUS")

	text_searcher.add_argument("--create_dirs", help = "Creates the directory structure. Default Y.", choices = ["Y", "N"], default = "Y")
	text_searcher.add_argument("--text_searcher", help = "Inspects blast XML files for a specified search term. Default Y.", choices = ["Y", "N"], default = "Y")

	# Extra arguments
	text_searcher.add_argument("--min_alignment", help = "The minimum alignment length for a 'good' alignment. Default 50.", default = 50)
	text_searcher.add_argument("--output_filename", help = "The output filename. Default Blastn_text_search.", default = "Blastn_text_search")
	
	
	return parser.parse_args()

################################################################################

if __name__ == '__main__':
	main()
