import subprocess, os, logging, pysam
import urllib
from subprocess import PIPE
from Bio import Entrez
from shutil import move as shmove
#TODO: Move the generation to fasta tool / fileHandler.
from Bio import SeqIO

LOG = logging.getLogger(__name__)
LOG.addHandler(logging.NullHandler())

def runGzip(file: str):
    subprocess.run(["pigz", file])

def fetchSRA(output_folder: str, accession: str):
    LOG.info(f"Fetching {accession}")
    #.strip is added due to trailing newlines.
    #https://blog.dalibo.com/2022/09/12/monitoring-python-subprocesses.html
    cmd = ["fasterq-dump", "--seq-defline", "@$sn[_$rn]/$ri", "-S", "-O", output_folder, accession.strip()]
    with subprocess.Popen(cmd, stdout= PIPE, stderr=subprocess.PIPE, text=True) as proc:
        errs = []
        for line in proc.stderr:
            LOG.warn(line.strip())
        stdout, _ = proc.communicate()
    result = subprocess.CompletedProcess(cmd, proc.returncode, stdout, "\n".join(errs))
    #In future I would like to find a way for this to check the filesize of the accessions against the memory available.
    return result

def safeEntrez(db_type: str, rettype: str, id_list: list[str]):
    try:
        handle = Entrez.efetch(db = db_type, 
                               id = set(id_list), 
                               rettype = rettype)
    except urllib.error.HTTPError:
        LOG.error(f"Unable to find at least one accession in: {id_list}.")
        handle = None
    return handle
    
def fetchEntrez(id_list: list, email: str, 
                api = False, proxy = 3128,
                db_type = "nuc"):
    #To help with FERA proxy shenanigans.
    os.environ["https_proxy"] = f"http://webcache:{proxy}"
    Entrez.email = email
    if api:
        api_key = api
    if db_type == "nuc":
        handle = safeEntrez("nuccore", "fasta", id_list)
    else:
        handle = safeEntrez("protein", "fasta", id_list)
    return handle
 
def getNumMappedReads(bwa_file: str) -> str:
    num_mapped_reads = subprocess.Popen(["samtools", "view", "-F", "0x04", 
                                        "-c", f"{bwa_file}"], 
                                         stdout = PIPE).communicate()[0]
    num_mapped_reads = num_mapped_reads.strip()
    num_mapped_reads = num_mapped_reads.decode("utf8")
    return num_mapped_reads

def outputSamHist(sorted_file: str, out_file: str):
    subprocess.run(["samtools", "coverage", sorted_file, "-m", "-o", out_file])
    
def runBedtools(out_file: str, bam: str):
    with open(out_file, "wb") as bg:
        subprocess.run(["bedtools", "genomecov", "-bg", "-ibam", bam],
                          stdout = bg)
    subprocess.run(["pigz", out_file])

def runBwa(fa: str, bwa_reads: list[str], out_file: str, threads = 12):
    index_result = subprocess.run(["bwa-mem2", "index", fa], capture_output=True)
    LOG.info(index_result.stdout)
    proc_call = ["bwa-mem2", "mem", "-v", "3", "-t", str(threads), fa]
    proc_call.extend(bwa_reads)
    LOG.debug(proc_call)
    with open(out_file, "w+") as sam:
        subprocess.run(proc_call, stdout=sam)
    
def samSort(bam_file: str, sam_file: str, mapq: int, flag: int):
    with open(bam_file, "w+") as bam:
        subprocess.run(["samtools", "view", "-q", str(mapq), "-F", str(flag), 
                        "-bS", sam_file], stdout = bam)
    pysam.sort("-o", bam_file, bam_file)
    subprocess.run(["samtools", "index", bam_file])
    idx_txt_out = os.path.splitext(sam_file)[0] + "_stats.txt"
    with open(idx_txt_out, "w+") as txt:
        subprocess.run(["samtools", "idxstats", bam_file], stdout = txt)
                                 
def runPfam(fasta_file, outfile, db_dir):
    with open(outfile, "w") as output:
        subprocess.run(["pfam_scan.pl", "-fasta", fasta_file, "-dir", db_dir, 
        "-json", "pretty"], stdout = output)
        
def runBlast2Rma(file: str, outdir: str, db: str, reads: str, 
                 blast_kind = "BlastN"):
    subprocess.run(["blast2rma", "-i", file, "-f", "BlastXML", "-o", outdir, 
                    "-ms", "75", "-sup", "1", "-a2t", db, "-bm", blast_kind, 
                    "-r", reads])
                    
def runRma2Info(filename, outfile):
    with open(outfile, "w") as output:
            subprocess.run(["rma2info", "--in", filename, "-vo", "-n", "-r2c", 
                            "Taxonomy", "-r", "-u", "false", "-v"], stdout = output)

def fastQC(threads: int, in_file: str, out_dir: str):
    subprocess.run(["fastqc", "-t", str(threads), in_file, "-o", out_dir])
    
def multiQC(input_dir: str, output_dir: str):
    subprocess.run(["multiqc", input_dir, "-o", output_dir])

def runBbduk(in_R1: str, in_R2: str, out_R1: str, out_R2: str,
             min_len: int, adapters: str, min_q: int):
    subprocess.run(["bbduk.sh", 
    f"in1={in_R1}", f"in2={in_R2}", f"out1={out_R1}", f"out2={out_R2}",
    f"minlen={min_len}", "ktrim=r", "k=23", "mink=11", "hdist=1",
    f"ref={adapters}", "qtrim=r", f"trimq={min_q}"])
   
def runTrinity(in_reads: list[str], out_file: str, mem: str, threads: int):
    args = ["Trinity", "--seqType", "fq", "--max_memory", mem]
    if len(in_reads) > 1:
        args.extend(["--left", in_reads[0], "--right", in_reads[1]])
    else:
        args.extend(["--single", in_reads[0]])
    args.extend(["--CPU", str(threads), 
                "--full_cleanup", 
                "--output", out_file])
    proc = subprocess.run(args, stdout = subprocess.PIPE)
    return trinity_proc.stdout

def mmseqs2(in_file: str, output: str, perc, threads: int, tmp: str) -> None:
    subprocess.run(["mmseqs", "easy-cluster", in_file, output, tmp,
                    "-c", str(perc), "--threads", str(threads), 
                    "-v", "1"])

def runNhmmer(sample_out: str, hmm, sample_in: str):
    subprocess.run(["nhmmer", "-o", f"{sample_out}.hmm", 
                    "--tblout", f"{sample_out}.tbl", hmm, sample_in])
                    
def bedtoolsWriteFasta(sample_out: str, sample_in: str):
    subprocess.run(["bedtools", "getfasta", "-s", "-fo", 
                    f"{sample_out}.fasta", "-fi", sample_in,
                    "-bed", f"{sample_out}.bed"])

def runSpades(reads: list[str], out_dir: str):
    args = ["spades.py", "-o", out_dir, "--careful"]
    if len(reads) > 1:
        args.extend(["-1", reads[0], "-2", reads[1]])
    else:
        args.extend(["-s", reads[0]])
    subprocess.run(args)
    
def runMakeblastdb(in_fasta: str, db_name: str, db_type = "nucl"):
    subprocess.run(["makeblastdb", "-in", in_fasta, 
                    "-parse_seqids", "-out", db_name, "-dbtype", db_type])