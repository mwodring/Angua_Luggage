import subprocess
from .AnguaUtils import back_mapper
from Bio import Entrez
import os

def backMapToBedGraph(trimmed_dir: str, outputdir: str, ref_file: str):
    just_filename = os.path.basename(ref_file)
    ref_sample_name = "_".join(just_filename.split("_")[:4])
    seqHandler.splitBbduk(trimmed_dir)
    for d in os.listdir(trimmed_dir):
        if d == ref_sample_name:
            back_mapper(input_dir = os.path.join(trimmed_dir, d), 
                        output_dir = outputdir, 
                        reference = ref_file, 
                        threads = 10, delim = "_", 
                        mapq = "0", flag = 2304, coverage = "Y", 
                        min_alignment_length = 50, softclip_perc = 1.0)  
    if os.path.exists(outputdir):
        for d in os.listdir(outputdir):
            for seq in os.scandir(os.path.join(outputdir, d)):
                if seq.name.endswith("_sort.bam"):
                    seqname = os.path.splitext(seq.name)[0]
                    #https://blog.liang2.tw/posts/2016/01/plot-seq-depth-gviz/#convert-sequencing-depth-to-bedgraph-format
                    #os.system(f"bedtools genomecov -bg -ibam {out_dir}/{d}/{seq.name} | gzip > {outputdir}/{d}/{seqname}.bedGraph.gz")
                    gz_outfile = os.path.join(outputdir, d, f"{seqname}.bedGraph.gz")
                    bg_outfile = os.path.join(out_dir, , d, seq.name)
                    bg_proc = subprocess.Popen(["bedtools", "genomecov", "-bg", "ibam", bg_outfile], stdout = PIPE)
                    gz_proc = subprocess.Popen(["gzip"], stdin = bg_proc, stdout = gz_outfile)

def fetchSRA(output_folder: str, accession: str):
    #Trying to move away from FastaKit logging to getting Angua to do it but unsure how to do that here.
    #self._logger.info(f"Fetching {accession}")
    #.strip is added due to trailing newlines.
    #https://blog.dalibo.com/2022/09/12/monitoring-python-subprocesses.html
    cmd = ["fasterq-dump", "-p", "-S", "-O", output_folder, accession.strip()]
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as proc:
        errs = []
        for line in proc.stderr:
            #self._logger.info(line.strip())
            errs.append(line)
        stdout, _ = proc.communicate()
    result = subprocess.CompletedProcess(cmd, proc.returncode, stdout, "\n".join(errs))
    #In future I would like to find a way for this to check the filesize of the accessions against the memory available.
    return result

 def getpercentage(a, b):
    return(a / b) * 100

def fetchEntrezFastas(id_list: list, email: str, outputdir: str, api = False, proxy = 3128, filename = "viruses"):
    #To help with FERA proxy shenanigans.
    os.environ["https_proxy"] = f"http://webcache:{proxy}"
    fasta_file = os.path.join(outputdir, filename)
    Entrez.email = email
    if api:
        api_key = api
    handle = Entrez.efetch(db = "nucleotide", 
                           id = set(id_list), 
                           rettype = "fasta")
    sequences = SeqIO.parse(handle, "fasta")
    with open(fasta_file, "w") as fasta_file:
        #SeqIO returns the count when it works, which is handy.
        count = SeqIO.write(sequences, fasta_file, "fasta")
    logger.info(f"{count} sequences found and written for {os.path.basename(fasta_file)}.")

def getNumMappedReads(bwa_file: str, sample_name: str, seq_name: str, bwa_folder: str) -> str:
    num_mapped_reads = subprocess.Popen(["samtools", "view", "-F", "0x04", "-c", f"{bwa_folder}/{bwa_file}"], stdout = subprocess.PIPE).communicate()[0]
    num_mapped_reads = num_mapped_reads.strip()
    num_mapped_reads = num_mapped_reads.decode("utf8")
    return num_mapped_reads
    
def runBwa(fa, bwa_reads, out_file):
    subprocess.run(["bwa", "index", fa])
    bwa_proc = subprocess.Popen(["bwa", "mem", "-v", "0", "-t", 
                                 "12", "-v", "3",
                                 fa, bwa_reads], stdout = PIPE)
    view_proc = subprocess.Popen(["view", "-S", "-b"], 
                                 stdin = bwa_proc, stdout = PIPE)
    sort_proc = subprocess.Popen(["samtools", "sort"], 
                                 stdout = out_file) 