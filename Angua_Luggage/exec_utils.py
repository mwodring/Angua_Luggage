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
    
def samToIndexedBam(in_sam: str, out_bam: str):
    view_proc = subprocess.Popen(["samtools", "view", "-q", "0", "-F", 
                                 "2304", "-bS", in_sam], stdout = subprocess.PIPE)
    with open(out_bam, "w+") as bam:
        bam.write(view_proc.stdout)
    sort_proc = subprocess.Popen(["samtools", "sort", 
                                "out_bam"])
    out_sorted = out_bam.replace(".bam", "_sort.bam")
    with open(out_sorted, "w+") as sort:
        sort.write(sort_proc.stdout)
    subprocess.run(["samtools", "index", out_sorted])
    return out_sorted
    
def runBedtools(out_file: str, bam: str):
    proc = subprocess.run(["bedtools", "genomecov", "-bg", "-ibam", bam],
                          stdout = PIPE)
    with open(out_file, "wb") as bg:
        bg.write(proc.stdout)
    subprocess.run(["pigz", out_file])

def runBwa(fa: str, bwa_reads: list[str], out_file: str):
    index_result = subprocess.run(["bwa-mem2", "index", fa], capture_output=True)
    LOG.info(index_result.stdout)
    proc_call = ["bwa-mem2", "mem", "-v", "0", "-t", 
                 "12", "-v", "3", 
                 fa]
    proc_call.extend(bwa_reads)
    bwa_proc = subprocess.run(proc_call, stdout = PIPE)
    with open(out_file, "wb") as sam:
        sam.write(bwa_proc.stdout)

def samSort(bam_file: str, sam_file: str):
    pysam.sort("-o", sam_file, bam_file)
                                 
def runPfam(fasta_file, outfile, db_dir):
    with open(outfile, "w") as output:
        subprocess.run(["pfam_scan.pl", "-fasta", fasta_file, "-dir", db_dir, 
        "-json", "pretty"], stdout = output)
        
def runBlast2Rma(file, outdir, db, reads, blast_kind = "BlastN"):
    subprocess.run(["blast2rma", "-i", file, "-f", "BlastXML", "-o", outdir, 
                    "-ms", "75", "-sup", "1", "-a2t", db, "-bm", blast_kind, 
                    "-r", reads])
                    
def runRma2Info(filename, outfile):
    with open(outfile, "w") as output:
            subprocess.run(["rma2info", "--in", filename, "-vo", "-n", "-r2c", 
                            "Taxonomy", "-r", "-u", "false", "-v"], stdout = output)