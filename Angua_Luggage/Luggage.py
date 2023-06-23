# -*- coding: utf-8 -*-
"""
Created on Fri May 19 15:51:33 2023

@author: mwodring
"""

import os
import re
import subprocess
import pandas as pd

class seqHandler:
    def __init__(self, folder: str, folder_type = "misc"):
        self.resetFolders()
        self.addFolder(folder_type, folder)
        self.toolBelt = toolBelt()

    def resetFolders(self):
        self._folders = {}

    def addFolder(self, folder_type: str, folder: str):
        self._folders[folder_type] = folder
    
    def getFolder(self, folder_type: str):
        if self._folders[folder_type]:
            return self._folders[folder_type]
        raise ValueError("No folder of this type.")
    
    def getFiles(self, folder_type: str, file_end = "*"):
        folder_name = self.getFolder(folder_type)
        for file in os.listdir(folder_name):
            if re.search(file, f".*[{file_end}]$"):
                yield f"{folder_name}/{file}"
                
    def mergeCSVOutput(self):
        out_csv = csvHandler(list(blastTool.header))
        for file in self.getFiles(folder_type = "csv", file_end = ".textsearch.csv"):
            out_csv.addCSVContents(file, sample = True)
        self.merged_csv = out_csv.mergeCSVOutput(self.getFolder("csv"))
       
    def addMappedToCSV(self, csv_file = ""):
        if not csv_file and self.merged_csv:
            csv_file = self.merged_csv
        out_csv = csvHandler(["sample", "species", "read_no"])
        for file in self.getFiles("bwa", ".tsv"):
            out_csv.addTSVContents(file)
        out_csv.outputMappedReads(folder_name = self.getFolder("csv"), csv_file = csv_file)
   
    def renameSRA():
        for file in self.getFiles(folder_type = "raw", file_end = "_1.fastq.gz"):
            entrezHandler.renameSRA(file, self.getFolder('raw'))
    
    def fetchSRAList(self, outout_folder: str, SRA_file: str):
        results = []
        self.addFolder(folder = output_folder, folder_type = "raw")
        with open(SRA_file, "r") as accessions:
            for accession in accessions.readlines():
                result = entrezHandler.fetchSRA(output_folder, accession)
                results.append(result)
                
    def textSearchToFasta(self, outputdir: str, email: str):
        hit_dir = os.path.join(outputdir, "hit_fastas")
        self.addFolder("ref", hit_dir)
        for file in self.getFiles("csv", ".textsearch.csv"):
            filename = os.path.splitext(os.path.basename(file))[0]
            sample_name = "_".join(filename.split("_")[:-3])
            accessions = csvHandler.getCSVAccesions(file)
            self.fasta_file = entrezHandler.fetchEntrezFastas(id_list = all_accessions, outputdir = hit_dir, email = email, filename = f"{os.path.splitext(sample_name)[0]}_accessions.fasta")
        #Might be worth naming these based on the samples they come from rather than the full csv.
    
    #This only works for a specific format atm. How to go about this?  
    def getSampleName(file: str):
        sample = os.path.splitext(os.path.basename(file))[0]
        sample = sample.split(".")[0]
        return sample
    
    @staticmethod
    def getpercentage(a, b):
        return(a / b) * 100

class bwaHandler():
    def __init__():
        pass
    
    def getNumMappedReads(bwa_file: str, sample_name: str, seq_name: str) -> str:
        num_mapped_reads = subprocess.Popen(["samtools", "view", "-F", "0x04", "-c", f"{self._folders['bwa']}/{bwa_file}"], stdout = subprocess.PIPE).communicate()[0]
        num_mapped_reads = num_mapped_reads.strip()
        num_mapped_reads = num_mapped_reads.decode("utf8")
        return num_mapped_reads
        
class entrezHandler():
    def __init__():
        pass
        
    @staticmethod    
    def renameSRA(file, raw_folder) -> int:
        sample_num = 1
        #This will need testing!
        whole_sample_name = seqHandler.getSampleName(file)
        samplename = whole_sample_name.split("_")[0]   
        #Uses the fact the SRR reads are stored under the 'raw' folder type in seqHandler. Extension is based on what Angua anticipates.
        new_filename = os.path.join(raw_folder, f"{samplename}_S{sample_num}_L001_R1_001.fastq.gz")
        os.rename(file, new_filename)
        #Same for second read. This just ensures the reads get the same sample name. 
        new_filename_R2 = new_filename.replace("R1", "R2")
        file_R2 = file.replace("_1", "_2")
        os.rename(file_R2, new_filename_R2)
        #Last of all increment the sample number.
        sample_num += 1
        #Returns the sample number for logging purposes if so desired.
        return sample_num
    
    @staticmethod
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
        
    @staticmethod
    def fetchEntrezFastas(id_list: list, email: str, outputdir: str, api = False, proxy = 3128, filename = "viruses"):
        #To help with FERA proxy shenanigans.
        os.environ["https_proxy"] = f"http://webcache:{proxy}"
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)
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
        return fasta_file
        
class csvHandler(seqHandler):
    def __init__(self, header_df: list):
        self.df_all = []
        self.header_df = header_df
    
    @staticmethod
    def getSampleName(file: str):
       sample = super().getSampleName(file)
       return sample
    
    @statismethod
    def getCSVAccessions(csv: str) -> DataFrame:
        df = pd.read_csv(file, sep = ",")
        all_accessions = df['NCBI accession'].to_list()
        return all_accessions
    
    def addCSVContents(self, file: str, sample = False):
        tmp_df = pd.read_csv(file, sep = ",")
        if sample:
            sample_name = self.getSampleName(file)
            tmp_df["sample"] = sample_name
        self.df_all.append(tmp_df)
     
    def mergeCSVOutput(self, folder: str) -> str:
        df_merged = pd.concat(self.df_all)
        self.header_df.append("sample")
        df_merged.columns = self.header_df
        all_csv = f"{folder}/all_samples.csv"
        df_merged.to_csv(all_csv)
        return all_csv
        
    def addTSVContents(self, file: str):
        tmp_df = pd.read_csv(file, sep = "\t", header = None)
        self.df_all.append(tmp_df)
    
    #This could probably be refactored to be more multi-use and less fragile tbh.
    def outputMappedReads(self, folder_name: str, csv_file: str):
        all_mapped_df = pd.concat(self.df_all)
        all_mapped_df.columns = self.header_df
        merged_df = pd.read_csv(csv_file, sep = ",")
        self.header_df = list(merged_df.columns).append("read_no")
        full_df = all_mapped_df.merge(merged_df, how = "outer")
        full_df.drop("Unnamed: 0", axis=1)
        test_csv_file = os.path.join(folder_name, "test_csv.csv")
        full_df.to_csv(test_csv_file, index = False)
   
class toolBelt:
    def __init__(self):
        #May be wiser to set up a way of doing (x)_tools on the fly? Or a YAML and a for loop?
        self.fasta_tools = []
        self.blast_tools = []
        self.pfam_tools = []
        self.bwa_tools = []
        self.rma_tools = []
        
    def addFastaTool(self, filename: str):
        self.fasta_tools.append(fastaTool())
        
    def addBlastTool(self, filename: str):
        self.blast_tools.append(blastTool())

class Tool:
    def __init__(self):
        pass

class fastaTool(Tool):
    def __init__(self):
        pass
    
class blastTool(Tool):
    def __init__(self):
        pass
    
class rmaTool(Tool):
    def __init__(self):
        pass

class bwaTool(Tool):
    def __init__(self):
        pass