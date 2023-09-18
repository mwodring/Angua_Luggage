import subprocess

#Note to self - have all of these return their stdout?

def fastQC(threads: int, input_dir: str, output_dir: str):
    subprocess.run(["fastqc", "-t", {threads}, f"{input_dir}/*", "-o", output_dir])
    
def multiQC(input_dir: str, output_dir: str):
    subprocess.run(["multiqc", input_dir, "-o", output_dir])

def runBbduk(in_R1: str, in_R2: str, out_R1: str, out_R2: str
             min_len: int, adapters: str, min_q: int):
    subprocess.run(["bbduk.sh", 
    f"in1={trimmer_input_R1}", "in2={trimmer_input_R2}", 
    f"out1={trimmer_output_R1}", f"out2={trimmer_output_R2}",
    f"minlen={min_len}", "ktrim=r", "k=23", "mink=11", "hdist=1",
    f"ref={adapters}", "qtrim=r", f"trimq={min_q}"])
   
def runTrinity(in_R1: str, in_R2: str, out_file: str, mem: str, threads: int):
    trinity_proc = subprocess.run(["Trinity", "--seqType fq", "--max_memory", mem,
    "--left", in_R1, "--right", in_R2, "--CPU", threads, 
    "--full_cleanup", "--output", out_file], stdout = subprocess.PIPE)
    return trinity_proc.stdout

def mmseqs2(in_file: str, out_dir: str, perc, threads: int) -> None:
    tmp = os.path.join(os.path.dirname(in_file), "tmp")
    subprocess.run(["mmseqs", "easy-cluster", "-c", "perc", "--threads", threads, "-v", "0", in_file, out_file])