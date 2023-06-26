# -*- coding: utf-8 -*-
"""
Created on Fri May 19 15:51:33 2023

@author: mwodring
"""

import os
import re
import subprocess
import pandas as pd
import csv
from shutil import move as 
from collections import defaultdict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .utils import new_logger, count_calls, Cleanup, add_logger_outdir
import .exec_utils

import importlib.resources
from . import data

from xml.parsers.expat import ExpatError
#Angua_test needs changing to Python 3.9+ to use resources.files.
fa_yaml = importlib.resources.open_binary(data, "fastaTool.yaml")
fastaTool_dict = yaml.full_load(fa_yaml)

#seqHandler deals with file and folder management and keeps tabs on which fasta reflects which 
#blast, etcetera. It also deals with file renaming. 
class seqHandler:
    def __init__(self, folder: str, folder_kind = "misc"):
        self.resetFolders()
        self.addFolder(folder_kind, folder)
        self.toolBelt = toolBelt()
    
    #Set folders to an empty dict.
    def resetFolders(self):
        self._folders = {}
    
    #Adds a folder to track, with its 
    def addFolder(self, folder_kind: str, folder: str):
        if not os.path.exists(folder)
            os.mkdir(folder)
        self._folders[folder_kind] = folder
    
    def getFolder(self, folder_kind: str):
        if self._folders[folder_kind]:
            return self._folders[folder_kind]
        raise ValueError("No folder of this type.")
    
    def getFiles(self, folder_kind: str, file_end = "*"):
        folder_name = self.getFolder(folder_kind)
        for file in os.listdir(folder_name):
            if re.search(file, f".*[{file_end}]$"):
                yield os.path.join(folder_name, file)
                
    def addFasta(self, filename: str, seqs = None, frame = 1, ID = "N/A", species = "N/A"):
        self.toolBelt.addFastaTool(filename, seqs, frame, ID, species)
                
    def mergeCSVOutput(self):
        out_csv = csvHandler(list(blastTool.header))
        for file in self.getFiles(folder_kind = "csv", file_end = ".textsearch.csv"):
            out_csv.appendCSVContents(file, sample = True)
        self.merged_csv = out_csv.mergeCSVOutput(self.getFolder("csv"))
       
    def appendMappedToCSV(self, csv_file = ""):
        if not csv_file and self.merged_csv:
            csv_file = self.merged_csv
        out_csv = csvHandler(["sample", "species", "read_no"])
        for file in self.getFiles("bwa", ".tsv"):
            out_csv.addTSVContents(file)
        out_csv.outputMappedReads(folder_name = self.getFolder("csv"), csv_file = csv_file)
   
    def renameSRA():
        for file in self.getFiles(folder_kind = "raw", file_end = "_1.fastq.gz"):
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
    
    def fetchSRAList(self, outout_folder: str, SRA_file: str):
        results = []
        self.addFolder(folder = output_folder, folder_kind = "raw")
        with open(SRA_file, "r") as accessions:
            for accession in accessions.readlines():
                result = exec_utils.fetchSRA(output_folder, accession)
                results.append(result)
                
    def textSearchToFasta(self, outputdir: str, email: str):
        hit_dir = os.path.join(outputdir, "hit_fastas")
        self.addFolder("ref", hit_dir)
        for file in self.getFiles("csv", ".textsearch.csv"):
            filename = os.path.splitext(os.path.basename(file))[0]
            sample_name = "_".join(filename.split("_")[:-3])
            accessions = csvHandler.getCSVAccesions(file)
            fasta_filename = f"{os.path.splitext(sample_name)[0]}_accessions.fasta"
            exec_utils.fetchEntrezFastas(id_list = all_accessions, 
                                         outputdir = hit_dir, 
                                         email = email, 
                                         filename = fasta_filename)
            self.addFasta(fasta_filename)
        #Might be worth naming these based on the samples they come from rather than the full csv.
    
    #This only works for a specific format atm. How to go about this?  
    def getSampleName(file: str):
        sample = os.path.splitext(os.path.basename(file))[0]
        sample = sample.split(".")[0]
        return sample
        
    def makeTempFastas(self, raw_dir: str, sample_name: str, fasta_file: str) -> dict:
        self.addFolder("raw", raw_dir)
        tmp_dir = os.path.join(self._folders["ref"], "tmp")
        self.addFolder("tmp", tmp_dir)
        tmp_seqs_dict = {}
        for n, raw_reads in enumerate(self.getFiles("raw", ".fasta")):
            if sample_name in raw_reads:
                tmp_seqs_dict[sample_name] = {"raw" : raw_reads,
                                              "seq_to_tmp" : {}}
            self.addFasta(filename = fasta_file, seqs = SeqIO.parse(fasta_file, "fasta"))
            seq_names, tmp = self.toolBelt.process(fasta_file, "fasta", makeTempFasta, n = n, tmp_folder = self.getFolder["tmp"])
            tmp_seqs_dict[sample_name]["seq_to_tmp"].update(zip(seq_names, tmp))
        return tmp_seqs_dict
    
    def runBwaTS(self, raw_dir: str, bwa_dir: str) -> list:
        if not os.path.exists(bwa_dir):
            os.mkdir(bwa_dir)
        self.addFolder("bwa", bwa_dir)
        tsv_files = []
        all_samples_dict = {}
        for fasta in self.getFiles("ref", ".fasta"):
            sample_name = "_".join(os.path.splitext(os.path.basename(fasta))[0].split("_")[:-1])
            all_samples_dict.update(self.makeTempFastas(raw_dir, sample_name, fasta))
        for sample_name, details in all_samples_dict.items():        
            #bwa_reads = os.path.join(self._folders["raw"], input_reads)
            bwa_reads = details["raw"]
            for i in range(len(details["tmp"])):
                tmp_fa = details["tmp"][i]
                seq_name = details["seq_names"][i]
                subprocess.run(["bwa", "index", tmp_fa])
                out_file = f"{sample_name}_{seq_name}.sorted.bam"
                subprocess.call(f"bwa mem -v 0 -t 12 -v 3 {tmp_fa} {bwa_reads} | samtools view -S -b | samtools sort > {bwa_dir}/{out_file}", shell = True)
                #Just remove the folder.
                os.remove(tmp_fa)
                self.getMappedReads(out_file, sample_name, seq_name)
        #Presumably needs the .clean or whatever?
        Cleanup([self._folders["ref"]], [".amb", ".ann", ".bwt", ".pac", ".sa"])
        return tsv_files
    
    #Some of these static methods might be better moved to their own class.
    @staticmethod
    def splitBbduk(trimmed_dir):
        for file in os.scandir(trimmed_dir):
            if file.name.endswith("_R1.fastq.gz"):
                filename = "".join(file.name.split(".")[-3])
                R2_file = f"{filename.replace('_R1', '_R2')}.fastq.gz"
                os.mkdir(os.path.join(trimmed_dir, filename[:-3]))
                shmove(os.path.join(trimmed_dir, file.name), 
                       os.path.join(trimmed_dir, filename[:-3], 
                       file.name))
                shmove(os.path.join(trimmed_dir, R2_file), 
                       os.path.join(trimmed_dir, filename[:-3], 
                       R2_file))
    
    def getMappedReads(self, bwa_file: str, sample_name: str, seq_name: str) -> str:
        num_mapped_reads = bwaHandler.getNumMappedReads(bwa_file, sample_name, seq_name, self.getFolder("bwa")
        tsv_file = os.join(self.getFolder('bwa'), f"{sample_name}_{seq_name}.tsv")
        csvHandler.mappedReadsTSV(tsv_file, sample_name)
        return tsv_file
        
class csvHandler(seqHandler):
    def __init__(self, header_df: list):
        self.df_all = []
        self.header_df = header_df
    
    def getSampleName(file: str) -> str:
       sample = super().getSampleName(file)
       return sample
    
    @staticmethod
    def mappedReadsTSV(tsv_file: str, sample_name: str) -> None:
        with open(tsv_file, "w+", newline='') as tsv:
            csv_writer = csv.writer(tsv, delimiter = "\t", lineterminator = "\n")
            trunc_sample_name = sample_name.split(".")[0]
            species_name = "_".join(seq_name.split("_")[2:])
            csv_writer.writerow([trunc_sample_name, species_name, num_mapped_reads])
    
    @staticmethod
    def getCSVAccessions(csv: str) -> DataFrame:
        df = pd.read_csv(file, sep = ",")
        all_accessions = df['NCBI accession'].to_list()
        return all_accessions
    
    def appendCSVContents(self, file: str, sample = False):
        tmp_df = pd.read_csv(file, sep = ",")
        if sample:
            sample_name = self.getSampleName(file)
            tmp_df["sample"] = sample_name
        self.df_all.append(tmp_df)
     
    def mergeCSVOutput(self, folder: str) -> str:
        df_merged = pd.concat(self.df_all)
        self.header_df.append("sample")
        df_merged.columns = self.header_df
        all_csv = os.path.join(folder, "all_samples.csv")
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
        tool_kinds = ["fasta", "blast", "pfam", "rma"]
        self.tools = {tool : defaultdict([]) for tool in tool_kinds}
        
    def addFastaTool(self, filename, seqs, frame = 1, ID = "N/A", species = "N/A"):
        #This feels like it repeats some code. Better way maybe?
        if not seqs:
            seqs = SeqIO.parse(filename, "fasta")
            for seq in seqs:
                species = seq.ID
                self.fasta_tools[filename].append(fastaTool(seq, frame, ID, species))
                return
        if type(seqs) == SeqRecord:
            species = seqs.ID
        self.fasta_tools.update({filename : fastaTool(seqs, frame, ID, species)})
            
        
    #Double check how args kwargs works.    
    def process(self, filename: str, tool_kind: str, func, *args, **kwargs):
        try:
            func_to_call = getattr(self.tools[tool_kind][filename], func)
        except AttributeError:
            "Fatal error, no such function."
        return func_to_call(*args, **kwargs)
    
    def outputFasta(self, add_text: str, by_species = False):
        pass

#May or may not be needed as each tool is fairly bespoke.
class Tool:
    def __init__(self):
        pass

class fastaTool(Tool):
    def __init__(self, seqs, frame: int, ID: str, species: str):
        self.filename = filename
        self.ID = ID
        if type(seqs) == str:
            seqs = [SeqRecord(seq = seqs, id = "Unnamed fasta", description = "")
        self.seqs = seqs
        self.frame = frame
        self.species = species
    
    #Check how this works with *args and **kwargs!
    def makeTempFasta(self, tmp_folder: str, n = 1) -> tuple:
        seq_names = []
        tmp = []
        for i, seq in enumerate(seqs):
            seq_name = re.sub("[ :.,()/]", "_", seq.description)
            seq_names.append(seq_name)
            tmp_fa = os.path.join(tmp_folder, f"{sample_name}_seq_{n}{i}_tmp.fasta")
            tmp.append(tmp_fa)
            SeqIO.write(seq, tmp_fa, "fasta")
            return seq_names, tmp
    
class blastTool(Tool):
    def __init__(self):
        pass
    
class rmaTool(Tool):
    def __init__(self):
        pass