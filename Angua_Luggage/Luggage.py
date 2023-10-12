# -*- coding: utf-8 -*-
"""
Created on Fri May 19 15:51:33 2023

@author: mwodring
"""

import logging, os, re, subprocess, shutil
import pandas as pd
#TODO: have pandas do the csv in/out.
import csv
from collections import defaultdict

from Bio import SeqIO, BiopythonDeprecationWarning
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonDeprecationWarning)
    from Bio import SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
import dataclasses 
from dataclasses import dataclass, field
from collections.abc import Generator

from .utils import Cleanup, getSampleName, getpercentage, subSeqName
from .exec_utils import *

from xml.parsers.expat import ExpatError
from xml.etree.ElementTree import ParseError

import importlib.resources
from . import data

from rpy2 import robjects as r
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
import logging
rpy2_logger.setLevel(logging.ERROR)

LOG = logging.getLogger(__name__)
LOG.addHandler(logging.NullHandler())

class fileHandler:
    def __init__(self, dir_kind: str, init_dir: str, extend = 0):
        self._resetFolders()
        self.addFolder(dir_kind, init_dir)
        self._toolBelt = toolBelt()
        #Extend is the number of _ after the sample name.
        self.extend = extend
    
    #Set folders to an empty dict.
    def _resetFolders(self):
        self._dirs = {}
    
    #Might be a good use for @property?
    def addFolder(self, dir_kind: str, add_dir: str):
        if not os.path.exists(add_dir):
            os.mkdir(add_dir)
        self._dirs[dir_kind] = add_dir
        return self._dirs[dir_kind]
    
    def extendFolder(self, dir_kind: str, new_dir_kind: str, extension: str):
        old_dir = self.getFolder(dir_kind)
        new_dir = os.path.join(old_dir, extension)
        return self.addFolder(new_dir_kind, new_dir)
	
    def extendFolderMultiple(self, orig_dir_kind: str, 
								   dir_kinds: list, dir_names: list):
		#Zip?
        for dir_kind, dir_name in zip(dir_kinds, dir_names):
            self.extendFolder(orig_dir_kind, dir_kind, dir_name)
    
    def removeFolder(self, dir_kind: str):
        to_remove = self.getFolder(dir_kind)
        if to_remove:
            shutil.rmtree(to_remove)
            
    def getFolder(self, dir_kind: str):
        try:
            return self._dirs[dir_kind]
        except KeyError:
            LOG.error(f"Not tracking directory of kind: {dir_kind}.")
            return
    
    def getFiles(self, dir_kind: str, file_end = "") -> Generator[str]:
        dir_name = self.getFolder(dir_kind)
        if type(file_end) == list:
            for end in file_end:
                for file in os.listdir(dir_name):
                    if file.endswith(end):
                        yield(os.path.join(dir_name, file))
        else:
            for file in os.listdir(dir_name):
                if file.endswith(file_end):
                    yield(os.path.join(dir_name, file))
                
    def addFastaFile(self, filename: str, frame = 1, 
                     ID = "N/A", species = "N/A"):
        self._toolBelt.addFastaTool(filename = filename,
                                    frame = frame,
                                    ID = ID, 
                                    species = species)
    
    def findFastaBySample(self, sample_name: str, dir_kind: str):
        files = []
        for file in self.getFiles(dir_kind, [".fasta", ".fq", ".fq.gz",
                                             ".fastq", ".fastq.gz"]):
            if sample_name in file and not ".log" in file:
                files.append(file)
        return files
    
    #Simple but it's a function in case it later wants validation.
    def addBlast(self, filename: str, ictv = False, blast_type = "Blastn"):
        self._toolBelt.addBlastTool(filename, blast_type, ictv)
    
    def findBlastFiles(self, ictv = False, blast_type = "Blastn"):
        blasts = [file for file in self.getFiles("xml", ".xml")]
        samples = [getSampleName(file, self.extend) for file in blasts]
        for i, xml in enumerate(blasts):
            self.addBlast(xml, blast_type, ictv)
        LOG.info(f"Found {i} xml files.")
        return samples
    
    def findFastaFiles(self, look_dir = "contigs") -> Generator[str]:
        fastas = self.getFiles(look_dir, [".fasta", ".fq", ".fastq"])
        if not fastas:
            LOG.critical("No .fasta files in folder!")
        for fasta in fastas:
            self.addFastaFile(fasta)
    
class csvHandler():
    __slots__ = ("df_all", "header_df")
    def __init__(self, header_df: list):
        self.df_all = []
        self.header_df = header_df
    
    @staticmethod
    def mappedReadsTSV(tsv_file: str, sample_name: str,
                       species: str, num_mapped_reads: int) -> None:
        with open(tsv_file, "w+", newline='') as tsv:
            csv_writer = csv.writer(tsv, delimiter = "\t", 
                                    lineterminator = "\n")
            trunc_sample_name = sample_name.split(".")[0]
            csv_writer.writerow([trunc_sample_name, 
                                 species, num_mapped_reads])
    
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
        try:
            df_merged = pd.concat(self.df_all)
        except ValueError:
            LOG.warn("No suitable hits found in any file; no csv to output.")
            return
        if not "sample" in self.header_df:
            self.header_df.append("sample")
        df_merged.columns = self.header_df
        all_csv = os.path.join(dir_name, "all_samples_textsearch.csv")
        df_merged.to_csv(all_csv, index = False)
        return all_csv
        
    def appendTSVContents(self, tsv: str):
        tmp_df = pd.read_csv(tsv, sep = "\t", header = None)
        tmp_df.columns = self.header_df
        self.df_all.append(tmp_df)
    
    def makeAllFromTSV(self):
        return pd.concat(self.df_all)
        
    #This could probably be refactored to be more multi-use and less fragile tbh.
    def outputMappedReads(self, dir_name: str, csv_file: str):
        merged_df = pd.read_csv(csv_file, sep = ",")
        reads_full_df = self.makeAllFromTSV()
        reads_full_df = reads_full_df.rename(columns={"species" : "_species"})
        merged_df["_species"] = merged_df["species"].map(subSeqName)
        with_mapped_reads = merged_df.merge(reads_full_df,
                                            on=["_species", "sample"])
        with_mapped_reads = with_mapped_reads.drop_duplicates()
        with_mapped_reads = with_mapped_reads.drop("_species", axis=1)
        out_csv_file = os.path.join(dir_name, "all_samples_mapped.csv")
        with_mapped_reads.to_csv(out_csv_file, index = False, mode = "w+")
    
    @staticmethod
    def outputHitsCSV(header, out_file: str, rows: list):
        with open(out_file, 'w+', encoding='UTF8', newline='') as out_csv:
            csv_writer = csv.writer(out_csv)
            csv_writer.writerow(header) 
            csv_writer.writerows(rows)
            LOG.info(f"csv written to {out_file}.")
    
class toolBelt():
    tool_kinds = ("fasta", "blast", "orf", "rma", "pfam")
    __slots__ = ("tools")
    
    def __init__(self):
        self.tools = {tool : defaultdict(list) for tool in toolBelt.tool_kinds}
    
    def addFastaTool(self, filename, seqs = None, 
                     frame = 1, ID = "N/A", species = "N/A"):
        if not seqs:
            seq_type = "fastq" if filename.endswith(".fastq") or filename.endswith(".fq") else "fasta"
            seqs = SeqIO.parse(filename, seq_type)
            for seq in seqs:
                self.tools["fasta"][filename].append(fastaTool(
                                                     filename = filename,
                                                     seq = seq, 
                                                     frame = frame, 
                                                     contig_id = seq.id, 
                                                     species = species))
            return
        if type(seqs) == SeqRecord:
            contig_ID = seqs.ID
            self.tools["fasta"].update({filename : fastaTool(
                                                filename = filename,
                                                seq = seqs, 
                                                frame = frame, 
                                                contig_id = contig_ID, 
                                                species = species)})
        
    def addBlastTool(self, filename: str, ictv: bool, blast_type = "Blastn"):
        self.tools["blast"].update({filename : blastTool(filename, blast_type, ictv)})
    
    def addorfTool(self, contig_dir: str):
        new_tool = orfTool(contig_dir)
        self.tools["orf"] = new_tool
        return new_tool
    
    def addpfamTool(self, filename: str, outfile: str):
        new_tool = pfamTool(filename, outfile)
        self.tools["pfam"].update({filename : new_tool})
        return new_tool
    
    def addRmaTool(self, file: str, output: str, db: str, contigs: str, 
                         sample_name: str, blast_kind = ""):
        new_tool = rmaTool(file, output, db, contigs, sample_name, 
                           blast_kind = blast_kind)
        self.tools["rma"].update({file : new_tool})
        return new_tool
    
    def findHitContigNames(self):
        for tool in self.getAllTools("blast"):
            print(tool.hit)

    #To run a process on all tools of type in all files.
    def process_all(self, tool_kind: str, func: str, 
                    *args, **kwargs) -> Generator[callable]:
        for filename in self.tools[tool_kind].keys():
            yield filename, self.process(filename, tool_kind, func, 
                                         *args, **kwargs)
            
    #To run a process on all tools of type connected to a file.  
    def process(self, filename: str, tool_kind: str, func: str, 
                *args, **kwargs) -> Generator[callable]:
        chosen_tools = [self.tools[tool_kind][filename]]
        for tool in chosen_tools:
            yield tool.process(func_to_call, *args, **kwargs)
    
    def outputContigAll(self, out_dir: str, add_txt = "_hits"):
        for filename, tools in self.tools["fasta"].items():
            out_file = os.path.join(out_dir,
                                    f"{getSampleName(filename)}_{add_text}.fasta")
            with open(out_file, "w+") as fa:
                for tool in tools:
                    tool.output(fa)
            LOG.info(f"{filename} written.")
                
    #Note to self, write some functions to make this comprehension less Worse.            
    def outputContigBySpecies(self, out_dir: str):
        #This might need windows of some kind to avoid memory issues.
        for species in self.getUniqueSpecies():
            if species != "N/A" and species != "N_A":
                species_tools = [tool for tool in self.getAllTools("fasta") 
                                if tool.species == species]
                underscore_species = subSeqName(species)
                while len(underscore_species) > 150:
                    underscore_species = "_".join(underscore_species.split()[:-1])
                for tool in species_tools:
                    sample_name = getSampleName(tool.filename)
                    out_file = os.path.join(out_dir, 
                                            f"{sample_name}_{underscore_species}_contigs.fasta")
                    with open(out_file, "a") as fa:
                        tool.output(fa)
    
    def migrateFasta(self, in_file: str, out_file: str):
        tools = self.tools["fasta"][in_file]
        with open(out_file, "a") as fa:
            for tool in tools:
                tool.output(fa)
    
    def getToolsByName(self, kind: str, name: str):
        tools = []
        for tool in self.getAllTools(kind):
            if name in tool.filename:
                tools.append(tool)
        return tools
                
    def mapFastaToBlast(self):
        for filename, tool in self.tools["blast"].items():
            all_info = tool.getHitFastaInfo()
            sample_name = getSampleName(filename)
            for info in all_info:
                these_tools = self.getToolsByName("fasta", sample_name)
                if not these_tools:
                    continue
                else:
                    filename = these_tools[0].filename
                    self.labelFasta(filename,
                                    frame = info["Frame"],
                                    to_label = info["contig_id"],
                                    species = info["species"])
        
    def labelFasta(self, filename: str, frame: int, to_label: str, species: str):
            for tool in self.tools["fasta"][filename]:
                if tool.contig_id == to_label:
                    tool.updateSpecies(species)
    
    def makeTempFastas(self, filename: str, tmp_dir: str, 
                       sample_name: str) -> dict:
        seq_names, tmp_fas = [], []
        for i, tool in enumerate(self.tools["fasta"][filename]):
            seq_name = subSeqName(tool.seq.description)
            species = " ".join(tool.seq.description.split(" ")[1:])
            tool.updateSpecies(species)
            bef_path_seq_name = seq_name.split("path=")[0]
            tmp_file = os.path.join(tmp_dir, 
                                    f"{sample_name}_{bef_path_seq_name}_tmp.fasta")
            seq_names.append(species)
            tmp_fas.append(tmp_file)
            with open(tmp_file, "w+") as fa:
                tool.output(fa)
        return seq_names, tmp_fas
    
    def getHitsCSVInfo(self, tool_kind = "blast") -> Generator[str, list]:
        for filename, tool in self.tools[tool_kind].items():
            yield filename, tool.getHitCSVInfo()
    
    def getUniqueSpecies(self, tool_kind = "fasta") -> set:
        species = [tool.species for tool in self.getAllTools(tool_kind)]
        return set(species)
    
    #TODO: CHeck where yield from and generator returns are more appropriate.
    def getAllTools(self, tool_kind: str):
        if tool_kind == "fasta":
            return (tool for filename in self.tools["fasta"].values() for tool in filename)
        else:
            return (tool for tool in self.tools[tool_kind].values())
    
    def parseAlignments(self, header, search_params = None, get_all = False):
        all_queries_parsed, all_hits = 0, 0
        for tool in self.getAllTools("blast"):
            LOG.info(f"Parsing {tool.filename}.")
            queries_parsed, hits = tool.parseAlignments(header, 
                                                        search_params, get_all)
            all_queries_parsed += queries_parsed
            all_hits += hits
        return all_queries_parsed, all_hits
    
    def getORFs(self, contig_dir: str, aa_dir: str, nt_dir: str):
        this_orfTool = self.addorfTool(contig_dir)
        this_orfTool.getORFs(aa_dir, nt_dir)
    
    def runPfam(self, db_dir: str, fasta_file: str, outfile: str):
        this_pfamTool = self.addpfamTool(fasta_file, outfile)
        this_pfamTool.runPfam(db_dir)
    
    def getAnnotations(self, pfam_dir: str, ORFs_file: str, gff3 = True):
        pfamTool.parseJson(pfam_dir, ORFs_file, gff3)
    
    def plotAnnotations(self, pfam_grl_file: str, pfam_df_file: str, 
                              plot_img_dir: str, backmap_dir: str):
        self.tools["orf"].plotAnnotations(pfam_grl_file, pfam_df_file,
                                          plot_img_dir, backmap_dir)
    
    def sortFastaOnLen(self, in_fa: str, out_fa: str, min_len: int):
        fas = self.getToolsByName("fasta", in_fa)
        sorted_fas = [fa for fa in fas if fa.longerThan(min_len)]
        for fa in sorted_fas:
            with open(out_fa, "a") as file:
                fa.output(file)
                    
    def filterFasta(self, in_fa: str, out_fa: str, on: str, filt: any):
        self.addFastaTool(in_fa)
        if on == "len":
            self.sortFastaOnLen(in_fa, out_fa, filt)
    
    def blast2Rma(self, file: str, output: str, db: str, contigs: str, 
                        blast_kind: str, sample_name: str):
        self.addRmaTool(file, output, db, contigs, sample_name,
                        blast_kind = blast_kind)
    
    def mapFastaToRma(self):
        for tool in self.getAllTools("rma"):
            sample = getSampleName(tool.filename)
            fasta_filename = self.getToolsByName("fasta", sample)[0].filename
            for c, taxon in tool.contig_to_taxon.items():
                self.labelFasta(fasta_filename, 1, c, taxon[1])
            
    def getMeganReports(self, out_dir: str, sortby = "virus"):
        for tool in self.getAllTools("rma"):
            tool.Rma2Info(sortby, out_dir)

@dataclass
class Tool:
    filename: str
    
    def process(self, func: str, *args, **kwargs):
        try:
            func_to_call = getattr(self, func)
        except AttributeError:
            LOG.error("No such function.")
        return func_to_call(*args, *kwargs)

@dataclass(slots=True)
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
        self.species = subSeqName(self.species)
    
    def updateSpecies(self, species):
        self.species = species
        
    def output(self, output_stream):
        SeqIO.write(self.seq, output_stream, "fasta")
    
    def longerThan(self, min_len: int):
        if(len(self.seq) >= int(min_len)):
            return self
            
@dataclass(slots=True)
class blastTool(Tool):
    blast_type: str
    ictv: bool
    _queries: Generator = field(init = False)
              
    def __post_init__(self):
        self._queries = (SearchIO.parse(self.filename, "blast-xml") if self.ictv 
                         else self.parseNCBI())
    
    def parseNCBI(self):
        return NCBIXML.parse(open(self.filename, "rU"), debug = 0)

    #Need to consult documentation to type hint this stuff.
    def parseHitData(self, hit, query, header: list) -> dict:
        aln_info = self.parseICTVData(hit, query) if self.ictv else self.parseNCBIData(hit, query)
        return {title : aln_info[i] for i, title in enumerate(header)}
        
    @staticmethod
    def parseICTVData(hit, query) -> tuple:
        hsp = hit.hsps[0]
        accession = hit.id
        species = hit.description
        ungapped = hsp.hit_span - hsp.gap_num
        coverage = getpercentage(hsp.hit_span,
                                 len(hsp.query.seq))
        identity = getpercentage(hsp.ident_num, 
                                 hsp.aln_span)
        return (species, coverage, identity, 
                len(hsp.query.seq), hit.query_id, 
                accession, "N/A", str(hsp.query.seq), ungapped, hsp.bitscore)
      
    def parseNCBIData(self, hit, query) -> tuple:
        #Assumes no frame.
        alignment = query.alignments[0]
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
        splitnum = 3 if blast_letter.upper() == "N" else 1
        if blast_letter == "X":
            frame = hsp.frame
        accession = alignment.hit_id.split("|")[splitnum]        
        #Unused for now - gb is genbank, ref is refseq.
        #db_type = alignment.hit_id.split("|")[0]
        LOG.debug(alignment.hit_def)
        return (alignment.hit_def, coverage, identity, 
                query.query_length, query.query,
                accession, frame[0], hsp.query, ungapped, hsp.bits)
    
        #Current implementation of the header feels dorky.
    def parseAlignments(self, header: list, 
                        search_params = None, get_all = True) -> tuple:
        all_aln = {}
        try:
            for n, query in enumerate(self._queries):
                aln = (self._parseICTVAln(query, n, header) if self.ictv else self._parseNCBIAln(query, n, get_all, header))
                all_aln.update(aln)
        except ParseError as err:
            LOG.error(err)
        self.hits = (self.checkAlignments(all_aln, search_params) if search_params 
                     else all_aln)
        return n, len(self.hits)
            
    def _parseNCBIAln(self, query, n: int, get_all: bool, header: list) -> dict:
        aln = {}
        to_check = 0
        try:
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
            wl_species = any(term.upper() in alignment["species"].upper() 
                             for term in search_params.search_term)
            bl_species = any(term.upper() in alignment["species"].upper() 
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
            contig_name = hit["contig_name"]
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
        return [hit.values() for hit in self.hits.values()] if self.hits else None
    
class orfTool():
    def __init__(self, contig_dir: str):
        self.contig_dir = contig_dir
    
    @staticmethod
    def init_r():
        r_path = importlib.resources.path(data, "annotatr_cli.r")
        r_script = os.path.abspath(r_path)
        r.r.source(r_script)
        
    def getORFs(self, aa_dir: str, nt_dir: str):
        self.init_r()
        r_output = r.r['ORF_from_fasta'](self.contig_dir, aa_dir, nt_dir, 150)
        log = r_output.rx2("log")
        self.grl_file = os.path.join(self.contig_dir, "grl.rdata")
        LOG.info(log)
    
    def plotAnnotations(self, pfam_grl_file: str, pfam_df_file: str, 
                              plot_img_dir: str, backmap_dir: str):
        self.init_r()
        r.r['generate_orf_plots'](self.grl_file, 
                                  self.contig_dir, plot_img_dir, 
                                  pfam_grl_file, pfam_df_file, backmap_dir)
        
@dataclass
class pfamTool(Tool):
    outfile: str
       
    def runPfam(self, db_dir: str):
        runPfam(self.filename, self.outfile, db_dir)

    def pfam_to_gff3(self):
        print(self.pfam_df)
    
    @staticmethod
    def parseJson(pfam_dir: str, ORF_file: str, gff3 = True):
        pfam_output = r.r['parse_pfam_json'](pfam_dir, ORF_file)
        if gff3:
            pfamTool.pfam_to_gff3()

class rmaTool():
    def __init__(self, file: str, output: str, db: str, 
                 contigs: str, sample_name: str,
                 blast_kind = ""):
        rma_file = os.path.splitext(os.path.basename(file))[0] + ".rma6"
        self.filename = os.path.join(output, rma_file)
        self.sample = sample_name
        self.contig_to_taxon = {}
        if not os.path.exists(self.filename):
            runBlast2Rma(file, output, db, contigs, blast_kind)
        
    def Rma2Info(self, sortby: str, out_dir: str) -> dict:
        file_no_extension = os.path.splitext(self.filename)[0]
        self.rma_txt = os.path.join(out_dir, 
                                    f"{os.path.basename(file_no_extension)}_info.txt")
        if not os.path.exists(self.rma_txt):
            runRma2Info(self.filename, self.rma_txt)
        self.info_dict = self.findRmaInfo()
        return self.info_dict
            
    def findRmaInfo(self) -> dict:
        with open(self.rma_txt, 'r') as info:
            for line in info.readlines():
                contig_name, rank, *virus_name = line.split("\t")
                virus_name = " ".join(virus_name).strip()
                self.contig_to_taxon[contig_name] = [rank, virus_name]
        return self.contig_to_taxon
        
    def getHitCSVInfo(self):
        rows = []
        for c, taxon in self.contig_to_taxon.items():
            rows.append([self.sample, c, taxon[0], taxon[1]])
        return rows