# -*- coding: utf-8 -*-
"""
Created on Fri May 19 15:51:33 2023

@author: mwodring
"""

import logging
import os
import re
import subprocess
import pandas as pd
import csv
from shutil import move as shmove
from collections import defaultdict

from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML

from .utils import count_calls, Cleanup
import .exec_utils

LOG = logging.getLogger(__name__).addHandler(logging.NullHandler())

#fileHandler deals with file and folder management and keeps tabs on 
#which fasta reflects which blast, etcetera. It also deals with 
#file renaming. 
class fileHandler:
    def __init__(self, init_dir: str, dir_kind = "misc"):
        self._resetFolders()
        self.addFolder(dir_kind, init_dir)
        self._toolBelt = toolBelt()
    
    #Set folders to an empty dict.
    def _resetFolders(self):
        self._dirs = {}
    
    #Adds a folder to track, with its 
    def addFolder(self, dir_kind: str, add_dir: str):
        if not os.path.exists(add_dir)
            os.mkdir(add_dir)
        self._dirs[dir_kind] = add_dir
    
    def getFolder(self, dir_kind: str):
        if self._dirs[dir_kind]:
            return self._dirs[dir_kind]
        LOG.warn("Attempt to access non-existent directory.")
    
    def getFiles(self, dir_kind: str, file_end = "*"):
        dir_name = self.getFolder(dir_kind)
        for file in os.listdir(dir_name):
            if re.search(file, f".*[{file_end}]$"):
                yield os.path.join(dir_name, file)
                
    def addFastaFile(self, filename: str, frame = 1, 
                     ID = "N/A", species = "N/A"):
        self._toolBelt.addFastaTool(filename, frame, ID, species)
    
    #Simple but it's a function in case it later wants validation.
    def addBlast(self, filename: str, ictv = False):
        self._toolBelt.addBlastTool(filename, ictv)
    
    def findBlastTools(self, ictv = False):
        blasts = self.getFiles("xml", ".xml")
        if not blasts:
            #Raise an error here.
            LOG.critical("No .xml files in input folder.")
        for xml in blasts:
            self.addBlast(file, ictv)
            yield xml
    
    def findFastaTools(self, look_dir = "contigs"):
        fastas = self.getFiles(look_dir, ".fasta")
        if not blasts:
            LOG.critical("No .fasta files bla bla")
        for fasta in fastas:
            self.addFastaFile(fasta)
            yield fasta
        
    def parseAlignments(self, search_params = None):
        self._toolBelt.parseAlignments()
            
    def updateFastaInfo(self):
        info = {filename, info for filename, info in 
                self._toolBelt.process_all("blast", getHitFastaInfo)}
        for filename, fa in fas:
            self._toolBelt.labelFasta(filename, 
                                      frame = info["Frame"], 
                                      to_label = info["contig_id"], 
                                      species = info["Species"])
    
    def mergeCSVOutput(self):
        out_csv = csvHandler(list(blastTool.header))
        for file in self.getFiles(dir_kind = "csv", 
                                  file_end = ".textsearch.csv"):
            out_csv.appendCSVContents(file, sample = True)
        self.merged_csv = out_csv.mergeCSVOutput(self.getFolder("csv"))
    
    def appendMappedToCSV(self, csv_file = None):
        if not csv_file and self.merged_csv:
            csv_file = self.merged_csv
        out_csv = csvHandler(["sample", "species", "read_no"])
        for file in self.getFiles("bwa", ".tsv"):
            out_csv.appendTSVContents(file)
        out_csv.outputMappedReads(dir_name = self.getFolder("csv"), 
                                  csv_file = csv_file)
   
    def renameSRA():
        for file in self.getFiles(dir_kind = "raw", 
                                  file_end = "_1.fastq.gz"):
            sample_num = 1
            #This will need testing!
            whole_sample_name = fileHandler.getSampleName(file)
            samplename = whole_sample_name.split("_")[0]   
            #Uses the fact the SRR reads are stored under the 'raw' folder 
            #type in fileHandler. Extension is based on what Angua anticipates.
            new_filename = f"{samplename}_S{sample_num}_L001_R1_001.fastq.gz"
            full_filename = os.path.join(raw_dir, 
                                        filename)
            os.rename(file, full_filename)
            #Same for second read. This just ensures the reads get 
            #the same sample name. 
            new_filename_R2 = new_filename.replace("R1", "R2")
            file_R2 = file.replace("_1", "_2")
            os.rename(file_R2, new_filename_R2)
            #Last of all increment the sample number.
            sample_num += 1
            #Returns the sample number for logging purposes if so desired.
        return sample_num
    
    def fetchSRAList(self, output_dir: str, SRA_file: str):
        results = []
        self.addFolder(add_dir = out_dir, dir_kind = "raw")
        with open(SRA_file, "r") as accessions:
            for accession in accessions.readlines():
                result = exec_utils.fetchSRA(out_dir, accession)
                results.append(result)
    
    def hitsToCSV(self, add_text: ""):
        self.updateFastaInfo()
        csv_out_folder = os.path.join(self.getFolder("out"), csv)
        self.addFolder("csv", csv_out_folder)
        for filename, info in self._toolBelt.getHitsCSVInfo():
            sample_name = getSampleName(filename)
            out_file = os.path.join(csv_out_file, f"{sample_name}_{add_text}_.csv")
            csvHandler.outputHitsCSV(header = blastTool.header, rows = info)
    
    def hitContigsToFasta(self, add_text: str, by_species = False):
        if not by_species:
            self._toolBelt.outputContigsFastaAll(out_dir = self.getFolder("out"))
        else:
            self._toolBelt.outputContigFastaBySpecies(out_dir = self.getFolder("out"))
    
    def hitAccessionsToFasta(self, email: str):
        hit_dir = os.path.join(self.getFolder("out"), "hit_fastas")
        self.addFolder("accs", hit_dir)
        for file in self.getFiles("csv", ".textsearch.csv"):
            filename = os.path.splitext(os.path.basename(file))[0]
            sample_name = "_".join(filename.split("_")[:-3])
            accessions = csvHandler.getCSVAccesions(file)
            fa_filename = f"{os.path.splitext(sample_name)[0]}_accessions.fasta"
            exec_utils.fetchEntrezFastas(id_list = all_accessions, 
                                         out_dir = hit_dir, 
                                         email = email, 
                                         filename = fa_filename)
            self.addFastaFile(fa_filename)
        #Might be worth naming these based on the samples they come from rather 
        #than the full csv.
    
    #This only works for a specific format atm. How to go about this?  
    def getSampleName(file: str):
        sample = os.path.splitext(os.path.basename(file))[0]
        sample = sample.split(".")[0]
        return sample
        
    def makeTempFastas(self, raw_dir: str, 
                       sample_name: str, fasta_file: str) -> dict:
        self.addFolder("raw", raw_dir)
        tmp_dir = os.path.join(self._dirs["acc"], "tmp")
        self.addFolder("tmp", tmp_dir)
        tmp_seqs_dict = {}
        for n, raw_reads in enumerate(self.getFiles("raw", ".fasta")):
            if sample_name in raw_reads:
                tmp_seqs_dict[sample_name] = {"raw" : raw_reads,
                                              "seq_to_tmp" : {}}
            self.addFastaFile(filename = fasta_file, seqs = SeqIO.parse(fasta_file, 
                          "fasta"))
            seq_names, tmp = self._toolBelt.makeTempFasta(fasta_file, n = n,
                                                         tmp_dir = self.getFolder["tmp"])
            tmp_seqs_dict[sample_name]["seq_to_tmp"].update(zip(seq_names, tmp))
        return tmp_seqs_dict
    
    def runBwaTS(self, raw_dir: str) -> list:
        self.addFolder("raw", raw_dir)
        self.addFolder("bwa", os.path.join(self.getFolder("out"), bwa)
        tsv_files = []
        all_samples_dict = {}
        for fasta in self.getFiles("acc", ".fasta"):
            sample_name = "_".join(
                          os.path.splitext(
                          os.path.basename(fasta))[0]
                          .split("_")[:-1])
            #This need rewriting to be clearer.
            all_samples_dict.update(self.makeTempFastas(self.getFolder("raw"), 
                                                        sample_name, fasta))
        for sample_name, details in all_samples_dict.items():        
            bwa_reads = details["raw"]
            for i in range(len(details["tmp"])):
                tmp_fa = details["tmp"][i]
                seq_name = details["seq_names"][i]
                out_file = os.path.join(bwa_dir, 
                                        f"{sample_name}_{seq_name}.sorted.bam"")
                exec_utils.runBwa(tmp_fa, bwa_reads, out_file)
                #Just remove the folder - write a remove dir function.
                os.remove(tmp_fa)
                self.coverageToTSV(out_file, sample_name, seq_name)
        #Presumably needs the .clean or whatever?
        Cleanup([self._dirs["acc"]], [".amb", ".ann", ".bwt", ".pac", ".sa"])
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
    
    def coverageToTSV(self, bwa_file: str, 
                       sample_name: str, seq_name: str) -> str:
        num_mapped_reads = bwaHandler.getNumMappedReads(bwa_file, 
                                                        sample_name, seq_name, 
                                                        self.getFolder("bwa"))
        tsv_file = os.join(self.getFolder('bwa'), 
                           f"{sample_name}_{seq_name}.tsv")
        csvHandler.mappedReadsTSV(tsv_file, sample_name)
        return tsv_file
        
class csvHandler(fileHandler):
    __slots__ = ("df_all", "header_df")
    def __init__(self, header_df: list):
        self.df_all = []
        self.header_df = header_df
    
    def getSampleName(file: str) -> str:
       sample = super().getSampleName(file)
       return sample
    
    @staticmethod
    def mappedReadsTSV(tsv_file: str, sample_name: str) -> None:
        with open(tsv_file, "w+", newline='') as tsv:
            csv_writer = csv.writer(tsv, delimiter = "\t", 
                                    lineterminator = "\n")
            trunc_sample_name = sample_name.split(".")[0]
            species_name = "_".join(seq_name.split("_")[2:])
            csv_writer.writerow([trunc_sample_name, 
                                 species_name, num_mapped_reads])
    
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
     
    def mergeCSVOutput(self, dir_name: str) -> str:
        df_merged = pd.concat(self.df_all)
        self.header_df.append("sample")
        df_merged.columns = self.header_df
        all_csv = os.path.join(dir_name, "all_samples.csv")
        df_merged.to_csv(all_csv)
        return all_csv
        
    def appendTSVContents(self, file: str):
        tmp_df = pd.read_csv(file, sep = "\t", header = None)
        self.df_all.append(tmp_df)
    
    #This could probably be refactored to be more multi-use and less fragile tbh.
    def outputMappedReads(self, dir_name: str, csv_file: str):
        all_mapped_df = pd.concat(self.df_all)
        all_mapped_df.columns = self.header_df
        merged_df = pd.read_csv(csv_file, sep = ",")
        self.header_df = list(merged_df.columns).append("read_no")
        full_df = all_mapped_df.merge(merged_df, how = "outer")
        full_df.drop("Unnamed: 0", axis=1)
        test_csv_file = os.path.join(dir_name, "test_csv.csv")
        full_df.to_csv(test_csv_file, index = False)
    
    @staticmethod
    def outputHitsCSV(header: list, out_file: str, rows: list):
        with open(out_file, 'w+', encoding=UTF8, newline='') as out_csv:
            csv_writer = csv.writer(out_csv)
            csv_writer.writerow(header) 
            csv_writer.writerows(rows)
            LOG.info(f"csv written to {out_file}.")
   
class toolBelt(fileHandler):
    tool_kinds = ["fasta", "blast", "pfam", "rma"]
    __slots__ = ("tools")
    def __init__(self):
        self.tools = {tool : defaultdict(list) for tool in tool_kinds}
    
    def getSampleName(self, filename: str):
        super().getSampleName(filename)
    
    def addFastaTool(self, filename, seqs = None, 
                     frame = 1, ID = "N/A", species = "N/A"):
        #This feels like it repeats some code. Better way maybe?
        if not seqs:
            seqs = SeqIO.parse(filename, "fasta")
            for seq in seqs:
                contig_ID = seq.ID
                self.tools["fasta"][filename].append(fastaTool(seq, 
                                                     frame, contig_ID, species))
                return
        if type(seqs) == SeqRecord:
            contig_ID = seqs.ID
            self.tools["fasta"].update({filename : fastaTool(seqs, 
                                                             frame, contig_ ID, 
                                                             species)})
        
    def addBlastTool(self, filename: str, ictv: bool):
        self.tools["blast"].update({filename : blastTool(filename, ictv)})
        
    def process_all(tool_kind: str, func, *args, **kwargs):
        for filename in self.tools[tool_kind].keys():
            yield filename, process(filename, tool_kind, func, *args, **kwargs)
            
    #Double check how args kwargs works.    
    def process(self, filename: str, tool_kind: str, func, *args, **kwargs):
        try:
            func_to_call = getattr(self.tools[tool_kind][filename], func)
        except AttributeError:
            "Fatal error, no such function."
        return func_to_call(*args, **kwargs)
    
    def outputContigFastaAll(self, out_dir: str):
        for filename, tools in self.tools["fasta"].items():
            out_file = os.path.join(out_dir,
                                       f"{getSampleName(filename)}_{add_text}.fasta")
            with open(out_file, "w+") as fa:
                for tool in tools:
                    tool.output(fa)
            LOG.info(f"{filename} written with all hits for contig.")
                
    #Note to self, write some functions to make this comprehension less Worse.            
    def outputContigFastaBySpecies(self, out_dir: str):
        for filename in self.tools["fasta"]:
            for species in self.getUniqueSpecies():
                current_tools = [tool for tool in filename 
                                 if tool.seq == species]
                out_file = os.path.join(out_dir, f"{species}_contigs.fasta")
                with open(out_file, 'w+') as fa:
                    for tool in current_tools:
                        tool.output(fa)
            
    def labelFasta(filename: str, frame: int, to_label: str, species: str):
        fa_tools = self.tools["fasta"][filename]
        for tool in fa_tools:
            if tool.contig_id == to_label:
                change_tool.updateSpecies(species)

    def makeTempFasta(self, filename: str, tmp_dir: str, sample_name: str, n = 1) -> tuple:
    seqs = []
    tmp_fa = []
    for i, tool in enumerate(self.tools["fasta"][filename]):
        seq_name = tool.subSeqName(n, i, sample_name)
        tmp_file = os.path.join(tmp_dir, 
                                f"{sample_name}_seq_{n}f{i}_tmp.fasta")
        seq_names.append(seq_name)
        tmp.append(tmp_file)
        with open(tmp_file, "w+") as fa:
            tool.output(fa)
    return seq_names, tmp_file
    
    def getHitsCSVInfo(self):
        for filename in self.tools["blast"]:
            yield filename, self.tools["blast"][filename].getHitCSVInfo()
    
    def getUniqueSpecies(self):
        species = [tool.species for tool in self.getAllTools("fasta")]
        return set(species)
        
    def getAllTools(self, tool_kind: str):
        return [tool for tool in 
                filename for filename in self.tools["tool_kind"]]
    
#May or may not be needed as each tool is fairly bespoke.
class Tool:
    def __init__(self, filename):
        pass

class fastaTool(Tool):
    __slots__ = ("contig_id", "seq", "frame", "species")
    def __init__(self, seq, frame: int, ID: str, species: str):
        self.contig_id = ID
        if type(seq) == str:
            seq = [SeqRecord(seq = seq, 
                             id = "Unnamed_fasta", description = "")
        self.seq = seq
        self.frame = frame
        self.species = subSeqName(species)
    
    def updateSpecies(self):
        self.species, self.seq.description = self.subSeqName(species)
        
    def subSeqName(self, to_sub: None):
        to_sub = self.seq.description if not to_sub else to_sub
        seq_name = re.sub("[ :.,()/]", "_", to_sub)
        return seq_name
        
    def output(output_stream):
        SeqIO.write(self.seq, output_stream, "fasta")

class blastTool(Tool):
    __slots__ = ("header", "_queries", "blast_type", "hits")
    header = ("species", "query % coverage", 
              "% identity", "contig length", 
              "contig name", "NCBI accession",
              "Frame", "alignment", 
              "matched length",
              "bitscore")
    
    def __init__(self, filename: str, ictv: bool):
        with open(filename) as handle:
            self._queries = (SearchIO.parse(handle, "blast-xml") if ictv 
                             else NCBIXML.parse(handle, debug = 0))
            self.blast_type = "BlastN"

    def parseHitData(hit, query, ictv = False):
        aln_info = parseICTVData(hit, query) if ictv else parseNCBIData(hit, query)
        return {title : dict_values[i] for i, title in enumerate(self.header)}
        
    def parseICTVData(hit, query):
        hsp = hit.hsps[0]
        accession = hit.id
        species = re.sub("[ :.,()/]", "_", hit.description)
        ungapped = hsp.hit_span - hsp.gap_num
        coverage = exec_utils.getpercentage(hsp.hit_span,
                                            len(hsp.query.seq))
        identity = exec_utils.getpercentage(hsp.ident_num, 
                                            hsp.aln_span)
        return (species, coverage, identity, 
                len(hsp.query.seq), hit.query_id, 
                accession, "N/A", str(hsp.query.seq), ungapped, hsp.bitscore)
        
    def parseNCBIData(hit, query):
        #Assumes no frame.
        frame = "N.A"
        hsp = hit.hsps[0]
        ungapped = hsp.align_length - hsp.gaps
        coverage = exec_utils.getpercentage(ungapped,
                                            query.query_length)
        identity = exec_utils.getpercentage(hsp.identities, 
                                            hsp.align_length)
        #These go by the formatting outputted by NCBI - 
        #the accession number is in the ID at different places.
        blast_letter = self.blast_type[-1]
        splitnum = 3 if blast_letter == "N" else 3
        if blast_letter == "X":
            frame = hsp.frame
        accession = alignment.hit_id.split("|")[splitnum]  
        #Unused for now - gb is genbank, ref is refseq.
        #db_type = alignment.hit_id.split("|")[0]
        return (alignment.hit_def, coverage, identity, 
                query.query_length, query.query,
                accession, frame[0], hsp.query, ungapped, hsp.bits)
        
    def parseAlignments(self, search_params = None, get_all = True) -> int:
        all_aln = {}
        for n, query in enumerate(self._queries):
            aln = (self._parseICTVAln(query, n) if self.ictv 
                   else self._parseNCBIAln(query, n, get_all))
            all_aln.update(aln)
        self.hits = (checkAlignments(all_aln, search_params) if search_params 
                     else full_aln)
        return n, len(self.hits)
            
    def _parseNCBIAln(self, query, n: int, get_all: bool):
        try:
            #Fetches the application type if it doesn't already know.
            if self.blast_type == "None":
                self.blast_type = query.application
            if len(query.alignments) > 0:
            to_check = len(query.alignments) if get_all == True else 1
            for i in range(to_check):
                alignment = query.alignments[i]
                aln[f"q{n}a{i+1}"] = self.parseHitData(alignment, query, ictv)       
            #Catches empty xml entries and ignores them.
        except ExpatError:
            pass
        return aln
        
    def _parseICTVAln(self, query, n: int) -> dict:
        aln = {}
        if len(query.hits) > 0:
            for i, hit in enumerate(query.hits):
                aln[f"q{n}a{i}"] = self.parseHitData(hit, query, ictv)
        return aln
        
    def checkAlignments(self, aln, search_params: NamedTuple):
        checked_aln = {}
        for key, alignment in aln.items():
            wl_species = any(term in alignment["species"].upper() 
                             for term in search_params.search_terms)
            bl_species = any(term in alignment["species"].upper() 
                             for term in search_params.blacklist)
            species_correct = wl_species and not bl_species
            if (
                species_correct 
                and aln["contig_length"] > search_params.minlen
                and aln["bitscore"] > search_params.bitscore
            ):
                checked_aln[key] = (alignment 
        return checked_aln
    
    def getHitFastaInfo(self) -> list:
        info = []
        for hit in self.hits:
            species = hit["species"]
            contig_name = hit["contig name"]
            if contig_name == "":
                logger.warn(f"No alignment found for {hit}.")
                continue
            contig_id = contig_name.split(" ")[0]
            current_fasta = {"contig_id" : contig_id,
                            "Frame" : hit["Frame"],
                            "species" : species}
            info.append(current_fasta)
        return info 
        
    def getHitCSVInfo(self):
        no_hits = ["No hits found."] * len(header)
        return self.aln.values() if self.aln else no_hits
        
class rmaTool(Tool):
    def __init__(self):
        pass