# -*- coding: utf-8 -*-
"""
Created on Fri May 19 15:51:33 2023

@author: mwodring
"""

import logging, os, re, subprocess
import pandas as pd
#TODO: have pandas do the csv in/out.
import csv
from collections import defaultdict

from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
import dataclasses import dataclass, field
from collections.abc import Generator

from .utils import count_calls, Cleanup, getSampleName
from .exec_utils import *

LOG = logging.getLogger(__name__)
LOG.addHandler(logging.NullHandler())

class fileHandler:
    def __init__(self, dir_kind: str, init_dir: str):
        self._resetFolders()
        self.addFolder(dir_kind, init_dir)
        self._toolBelt = toolBelt()
    
    #Set folders to an empty dict.
    def _resetFolders(self):
        self._dirs = {}
    
    #Might be a good use for @property?
    def addFolder(self, dir_kind: str, add_dir: str):
        if not os.path.exists(add_dir):
            os.mkdir(add_dir)
        self._dirs[dir_kind] = add_dir
    
    def getFolder(self, dir_kind: str):
        try:
            return self._dirs[dir_kind]
        except KeyError:
            LOG.error("Not tracking directory of that kind.")
    
    def getFiles(self, dir_kind: str, file_end = "") -> Generator[str]:
        dir_name = self.getFolder(dir_kind)
        for file in os.listdir(dir_name):
            if file.endswith(file_end):
                yield(os.path.join(dir_name, file))
                
    def addFastaFile(self, filename: str, frame = 1, 
                     ID = "N/A", species = "N/A"):
        self._toolBelt.addFastaTool(filename, frame, ID, species)
    
    #Simple but it's a function in case it later wants validation.
    def addBlast(self, filename: str, ictv = False):
        self._toolBelt.addBlastTool(filename, ictv)
    
    def findBlastFiles(self, ictv = False):
        blasts = self.getFiles("xml", ".xml")
        #if not blasts:
            # TODO Raise an error here.
         #   LOG.critical("No .xml files in input folder.")
        for xml in blasts:
            self.addBlast(xml, ictv)
    
    def findFastaFiles(self, look_dir = "contigs") -> Generator[str]:
        fastas = self.getFiles(look_dir, ".fasta")
        if not blasts:
            LOG.critical("No .fasta files bla bla")
        for fasta in fastas:
            self.addFastaFile(fasta)
            yield fasta
        
class csvHandler():
    __slots__ = ("df_all", "header_df")
    def __init__(self, header_df: list):
        self.df_all = []
        self.header_df = header_df
    
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
    def getCSVAccessions(csv: str) -> pd.DataFrame:
        df = pd.read_csv(csv, sep = ",")
        all_accessions = df['NCBI accession'].to_list()
        return all_accessions
    
    def appendCSVContents(self, csv: str, sample = False):
        tmp_df = pd.read_csv(csv, sep = ",")
        if sample:
            sample_name = getSampleName(csv, extend = 1)
            tmp_df["sample"] = sample_name
        self.df_all.append(tmp_df)
     
    def mergeCSVOutput(self, dir_name: str) -> str:
        df_merged = pd.concat(self.df_all)
        self.header_df.append("sample")
        df_merged.columns = self.header_df
        all_csv = os.path.join(dir_name, "all_samples.csv")
        df_merged.to_csv(all_csv)
        return all_csv
        
    def appendTSVContents(self, tsv: str):
        tmp_df = pd.read_csv(tsv, sep = "\t", header = None)
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
    def outputHitsCSV(header, out_file: str, rows: list):
        with open(out_file, 'w+', encoding='UTF8', newline='') as out_csv:
            csv_writer = csv.writer(out_csv)
            csv_writer.writerow(header) 
            csv_writer.writerows(rows)
            LOG.info(f"csv written to {out_file}.")
    
class toolBelt():
    tool_kinds = ("fasta", "blast", "pfam", "rma")
    __slots__ = ("tools")
    
    def __init__(self):
        self.tools = {tool : defaultdict(list) for tool in toolBelt.tool_kinds}
    
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
                                                             frame, contig_ID, 
                                                             species)})
        
    def addBlastTool(self, filename: str, ictv: bool):
        self.tools["blast"].update({filename : blastTool(filename, ictv)})
    
    #To run a process on all tools of type in all files.
    def process_all(self, tool_kind: str, func: str, 
                    *args, **kwargs) -> Generator[Callable]:
        for filename in self.tools[tool_kind].keys():
            yield filename, self.process(filename, tool_kind, func, 
                                         *args, **kwargs)
            
    #To run a process on all tools of type connected to a file.  
    def process(self, filename: str, tool_kind: str, func: str, 
                *args, **kwargs) -> Generator[Callable]:
        chosen_tools = [self.tools[tool_kind][filename]]
        for tool in chosen_tools:
            yield tool.process(func_to_call, *args, **kwargs)
    
    def outputContigAll(self, out_dir: str):
        for filename, tools in self.tools["fasta"].items():
            out_file = os.path.join(out_dir,
                                    f"{getSampleName(filename)}_{add_text}.fasta")
            with open(out_file, "w+") as fa:
                for tool in tools:
                    tool.output(fa)
            LOG.info(f"{filename} written with all hits for contig.")
                
    #Note to self, write some functions to make this comprehension less Worse.            
    def outputContigBySpecies(self, out_dir: str):
        for filename in self.tools["fasta"]:
            for species in self.getUniqueSpecies():
                current_tools = [tool for tool in filename 
                                 if tool.species == species]
                out_file = os.path.join(out_dir, f"{species}_contigs.fasta")
                with open(out_file, 'w+') as fa:
                    for tool in current_tools:
                        tool.output(fa)
    
    def mapFastaToBlast(self):
        for filename, tool in self.tools["blast"].items():
            all_info = tool.getHitFastaInfo()
            for info in all_info:
                self.labelFasta(filename,
                                frame = info["Frame"],
                                to_label = info["contig_id"],
                                species = info["species"])
        
    def labelFasta(self, filename: str, frame: int, to_label: str, species: str):
        fa_tools = self.tools["fasta"][filename]
        for tool in fa_tools:
            if tool.contig_id == to_label:
                change_tool.updateSpecies(species)
    
    def makeTempFasta(self, filename: str, tmp_dir: str, 
                      sample_name: str, n = 1) -> tuple:
        seqs, tmp_fa = [], []
        for i, tool in enumerate(self.tools["fasta"][filename]):
            seq_name = tool.subSeqName(n, i, sample_name)
            tmp_file = os.path.join(tmp_dir, 
                                    f"{sample_name}_seq_{n}f{i}_tmp.fasta")
            seq_names.append(seq_name)
            tmp.append(tmp_file)
            with open(tmp_file, "w+") as fa:
                tool.output(fa)
        return seq_names, tmp_file
    
    def getHitsCSVInfo(self) -> Generator[str, list]:
        for filename, tool in self.tools["blast"].items():
            yield filename, tool.getHitCSVInfo()
    
    def getUniqueSpecies(self) -> set:
        species = [tool.species for tool in self.getAllTools("fasta")]
        return set(species)
    
    #TODO: CHeck where yield from and generator returns are more appropriate.
    def getAllTools(self, tool_kind: str) -> Generator[Tool]:
        if tool_kind == "fasta":
            return (tool for filename in self.tools["fasta"].values() for tool in filename)
        else:
            return (tool for tool in self.tools[tool_kind].values())
    
    def parseAlignments(self, header, search_params = None, get_all = False):
        all_queries_parsed, all_hits = 0, 0
        for tool in self.getAllTools("blast"):
            queries_parsed, hits = tool.parseAlignments(header, 
                                                        search_params, get_all)
            all_queries_parsed += queries_parsed
            all_hits += hits
        return all_queries_parsed, all_hits
        
@dataclass
class Tool:
    filename: str
    
    def process(self, func: str, *args, **kwargs):
        try:
            func_to_call = getattr(self, func)
        except AttributeError:
            LOG.error("No such function.")
        return func_to_call(*args, *kwargs)

@dataclass
class fastaTool(Tool):
    seq: str | SeqRecord
    frame: int
    contig_id: str
    species: str
    
    def __post_init__(self):
        if type(self.seq) == str:
            seq = [SeqRecord(seq = self.seq, 
                   id = "Unnamed_fasta", description = "")]
        self.seq = seq
        self.species = self.subSeqName(self.species)
    
    def updateSpecies(self, species):
        self.species, self.seq.description = self.subSeqName(species)
        
    def subSeqName(self, to_sub: None):
        to_sub = self.seq.description if not to_sub else to_sub
        seq_name = re.sub("[ :.,()/]", "_", to_sub)
        return seq_name
        
    def output(output_stream):
        SeqIO.write(self.seq, output_stream, "fasta")

@dataclass
class blastTool(Tool):
    ictv: bool
    blast_type: str = "Blastn"
    _queries: Generator = field(init = False)
              
    def __post_init__(self):
        self._queries = (SearchIO.parse(self.filename, "blast-xml") if ictv 
                         else NCBIXML.parse(handle, debug = 0))
    
    #Need to consult documentation to type hint this stuff.
    def parseHitData(self, hit, query, header: list) -> dict:
        aln_info = self.parseICTVData(hit, query) if self.ictv 
                   else self.parseNCBIData(hit, query)
        return {title : aln_info[i] for i, title in enumerate(header)}
        
    @staticmethod
    def parseICTVData(hit, query) -> tuple:
        hsp = hit.hsps[0]
        accession = hit.id
        species = re.sub("[ :.,()/]", "_", hit.description)
        ungapped = hsp.hit_span - hsp.gap_num
        coverage = getpercentage(hsp.hit_span,
                                            len(hsp.query.seq))
        identity = getpercentage(hsp.ident_num, 
                                            hsp.aln_span)
        return (species, coverage, identity, 
                len(hsp.query.seq), hit.query_id, 
                accession, "N/A", str(hsp.query.seq), ungapped, hsp.bitscore)
    
    @staticmethod    
    def parseNCBIData(hit, query) -> tuple:
        #Assumes no frame.
        frame = "N.A"
        hsp = hit.hsps[0]
        ungapped = hsp.align_length - hsp.gaps
        coverage = getpercentage(ungapped,
                                            query.query_length)
        identity = getpercentage(hsp.identities, 
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
        
        #Current implementation of the header feels dorky.
    def parseAlignments(self, header: list, 
                        search_params = None, get_all = True) -> tuple:
        all_aln = {}
        for n, query in enumerate(self._queries):
            aln = (self._parseICTVAln(query, n, header) if self.ictv 
                   else self._parseNCBIAln(query, n, get_all, header))
            all_aln.update(aln)
        self.hits = (self.checkAlignments(all_aln, search_params) if search_params 
                     else all_aln)
        return n, len(self.hits)
            
    def _parseNCBIAln(self, query, n: int, get_all: bool, header: list) -> dict:
        try:
            #Fetches the application type if it doesn't already know.
            if self.blast_type == "None":
                self.blast_type = query.application
            if len(query.alignments) > 0:
                to_check = len(query.alignments) if get_all == True else 1
            for i in range(to_check):
                alignment = query.alignments[i]
                aln[f"q{n}a{i+1}"] = self.parseHitData(alignment, query, 
                                                       header)       
            #Catches empty xml entries and ignores them.
        except ExpatError:
            pass
        return aln
        
    def _parseICTVAln(self, query, n: int, header) -> dict:
        aln = {}
        if len(query.hits) > 0:
            for i, hit in enumerate(query.hits):
                aln[f"q{n}a{i}"] = self.parseHitData(hit, query, header)
        return aln
        
        #Need to check type of argparse options.
    def checkAlignments(self, all_aln: dict, search_params) -> dict:
        if not all_aln:
            return {}
        checked_aln = {}
        for key, alignment in all_aln.items():
            wl_species = any(term in alignment["species"].upper() 
                             for term in search_params.search_term)
            bl_species = any(term in alignment["species"].upper() 
                             for term in search_params.blacklist)
            species_correct = wl_species and not bl_species
            correct = (species_correct 
                       and alignment["contig length"] > search_params.minlen
                       and alignment["bitscore"] > search_params.bitscore)
            if correct:
                checked_aln[key] = alignment
        return checked_aln
        
    def getHitFastaInfo(self) -> list:
        info = []
        for hit in self.hits.values():
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
        
    def getHitCSVInfo(self) -> list:
        #This causes the string to go across the whole .csv but it's so funny.
        return [hit.values() for hit in self.hits.values()] if self.hits else None

@dataclass
class rmaTool(Tool):
    pass