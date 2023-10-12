from .Luggage import fileHandler, csvHandler
from .utils import getSampleName, Cleanup, subSeqName
from .exec_utils import *
import json, importlib.resources
from . import data
import os, traceback, sys, shutil
from pathlib import Path

LOG = logging.getLogger(__name__)
LOG.addHandler(logging.NullHandler())

header_json = importlib.resources.open_binary(data, "header.json")
header = json.load(header_json)

#Could probably split off the parsing and csv functions?
class blastParser(fileHandler):
    def parseAlignments(self, search_params = None, 
                        header = header, get_all = False):
        return self._toolBelt.parseAlignments(header = header,
                                              search_params = search_params,
                                              get_all = get_all)
    
    def updateFastaInfo(self):
        self._toolBelt.mapFastaToBlast()
    
    def mergeCSVOutput(self, alt_header = ""):
        if alt_header:
            this_header = alt_header
        else:
            this_header = header
        out_csv = csvHandler(this_header)
        for file in self.getFiles(dir_kind = "csv", 
                                  file_end = ".textsearch.csv"):
            out_csv.appendCSVContents(file, sample = True)
        self.merged_csv = out_csv.mergeCSVOutput(self.getFolder("csv"))
    
    def appendMappedToCSV(self, csv_file = None):
        if not csv_file and hasattr(self, "merged_csv"):
            csv_file = self.merged_csv
        out_csv = csvHandler(["sample", "species", "read_no"])
        tsvs = list(Path(self.getFolder("bwa")).rglob("*.[tT][sS][vV]"))
        for tsv in tsvs:
            out_csv.appendTSVContents(tsv)
        out_csv.outputMappedReads(dir_name = self.getFolder("csv"), 
                                  csv_file = csv_file)
                                  
    def hitsToCSV(self, add_text = "", tool_kind = "blast",
                  alt_header = ""):
        if alt_header:
            this_header = alt_header
        else:
            this_header = header
        csv_out_folder = os.path.join(self.getFolder("out"), "csv")
        self.addFolder("csv", csv_out_folder)
        for filename, info in self._toolBelt.getHitsCSVInfo(tool_kind):
            if info:
                sample_name = getSampleName(filename, self.extend)
                out_file = os.path.join(csv_out_folder, 
                                        f"{sample_name}_{add_text}.textsearch.csv")
                csvHandler.outputHitsCSV(header = this_header, 
                                         rows = info, out_file = out_file)
            else:
                LOG.info(f"No suitable hits for {filename}.")
    
    def hitContigsToFasta(self):
        self.updateFastaInfo()
        out_dir = os.path.join(self.getFolder("out"), "contigs")
        self.addFolder("parsed_contigs", out_dir)
        self._toolBelt.outputContigBySpecies(out_dir)
    
    def hitAccessionsToFasta(self, email: str, db_type="N"):
        db_type = "nuc" if db_type == "N" else "prot"
        out_dir = self.extendFolder("out", "acc", "hit_fastas")
        for file in self.getFiles("csv", ".textsearch.csv"):
            sample_name = getSampleName(file)
            accessions = csvHandler.getCSVAccessions(file)
            fa_filename = f"{sample_name}_accessions.fasta"
            out = os.path.join(out_dir, fa_filename)
            handle = fetchEntrez(id_list = accessions, 
                                 email = email, 
                                 db_type = db_type)
            if handle:
                sequences = SeqIO.parse(handle, "fasta")
                with open(out, "w+") as fa:
                    #SeqIO returns the count when it works, which is handy.
                    count = SeqIO.write(sequences, fa, "fasta")
                    LOG.info(f"{count} sequences found and written for {sample_name}.")
                    self.addFastaFile(out)
        
    def makeTempFastas(self, sample_name: str, fasta_file: str) -> dict:
        seq_names, tmp_fas = self._toolBelt.makeTempFastas(
                                fasta_file,
                                tmp_dir = self.getFolder("tmp"),
                                sample_name = sample_name)
        return seq_names, tmp_fas 
    
    def runBwaTS(self, raw_dir: str, in_dir_type = "acc", extend = 0,
                 threads = 23, mapq = 0, flag = 2304, text_search = True) -> list:
        self.addFolder("raw", raw_dir)
        self.findFastaFiles("raw")
        self.extendFolder("out", "bwa", "bwa")
        self.extendFolder("out", "hist", "Histograms")
        all_trimmed = (file for file in self.getFiles("raw") if "R1" in file)
        tsv_files = []
        all_samples_dict = {}
        tmp_dir = self.extendFolder(in_dir_type, "tmp", "tmp")
        for fasta in self.getFiles(in_dir_type, ".fasta"):
            dir_name = os.path.splitext(os.path.basename(fasta))[0]
            self.addFastaFile(fasta)
            sample_name = getSampleName(fasta, extend = extend) if text_search else getSampleName(next(all_trimmed), extend = 1)
            seq_names, tmp_fas = self.makeTempFastas(sample_name, fasta)
            all_samples_dict[sample_name] = dict(zip(seq_names, tmp_fas))
        for sample_name, seq_to_tmp in all_samples_dict.items():
            bwa_reads = self.findFastaBySample(sample_name, dir_kind = "raw")
            for seq_name, tmp_fa in seq_to_tmp.items():
                underscore_seq_name = subSeqName(seq_name)
                out_dir = self.extendFolder("bwa", dir_name, dir_name)
                out_file = os.path.join(out_dir, 
                                        f"{sample_name}_{underscore_seq_name}.sam")
                sorted_file = os.path.join(out_dir, 
                                           f"{sample_name}_{underscore_seq_name}.sorted.bam")
                if not os.path.exists(sorted_file):
                    runBwa(tmp_fa, bwa_reads, out_file, threads = threads)
                    samSort(sorted_file, out_file, mapq, flag)
                    hist_dir = self.extendFolder("hist", dir_name, dir_name)
                    hist_file = os.path.join(hist_dir, 
                                             f"{sample_name}_{underscore_seq_name}_hist.txt")
                    outputSamHist(sorted_file, hist_file)
                    self.coverageToTSV(out_file, sample_name, seq_name)
        self.removeFolder("tmp")
        Cleanup(self.getFolder(in_dir_type), [".amb", ".ann", ".bwt", ".pac", ".sa"])
        Cleanup(self.getFolder("bwa"), [".sam"])
        return tsv_files
            
    @staticmethod
    def coverageToTSV(bwa_file: str, 
                      sample_name: str, seq_name: str) -> str:
        bwa_dir = os.path.dirname(bwa_file)
        num_mapped_reads = getNumMappedReads(bwa_file)
        if num_mapped_reads == 0:
            LOG.info(f"No reads mapped for {seq_name} to {sample_name}.")
        seq_name = subSeqName(seq_name)
        tsv_file = os.path.join(bwa_dir, 
                                f"{sample_name}_{seq_name}.tsv")
        csvHandler.mappedReadsTSV(tsv_file, sample_name, seq_name, num_mapped_reads)
        return tsv_file
        
class SRA(fileHandler):
    def fetchSRAList(self, SRA_file: str):
        count = 0
        with open(SRA_file, "r") as accessions:
            to_fetch = accessions.readlines()
        for accession in to_fetch:
            accession = accession.strip()
            complete_file = os.path.join(self.getFolder("raw"),
                                         f"{accession}_SX_L001_R1_001.fastq.gz")
            if not os.path.exists(complete_file):
                proc_info = fetchSRA(self.getFolder("raw"), 
                                     accession)
                LOG.info(proc_info.stdout)
            else:
                LOG.info(f"{accession} already exists.")
            count += self.renameSRA(accession)
        return count 
                       
    def renameSRA(self, sample_name: str) -> int:
        samples = (file for file in self.getFiles("raw", ".fastq") if sample_name in file)
        for i, sample in enumerate(samples):
            new_filename = f"{sample_name}_SX_L001_R{i+1}_001.fastq"
            full_filename = os.path.join(self.getFolder("raw"), 
                                         new_filename)
            os.rename(sample, full_filename)
        return 1
    
    def pigzAll(self):
        for file in self.getFiles("raw", ".fastq"):
            runGzip(file)

class Annotatr(fileHandler):
    def generateorfTools(self):
        self._toolBelt.getORFs(self.getFolder("contigs"), 
                                 self.getFolder("aa"), 
                                 self.getFolder("ORF_nt"))
    
    def getORFs(self, out_dir: str, contig_dir = None):
        self.addFolder("ORFs", out_dir)
        if contig_dir:
            self.addFolder("contigs", contig_dir)
        self.extendFolder("ORFs", "aa", "aa")
        self.extendFolder("ORFs", "ORF_nt", "nt")
        self.generateorfTools()
        self.ORF_file = os.path.join(self.getFolder("contigs"), "ORFs.rdata")
        self.grl_file = os.path.join(self.getFolder("contigs"), "grl.rdata")
            
    def runPfam(self, db_dir: str, add_text = "_virus"):
        pfam_dir = self.extendFolder("ORFs", "pfam", "pfam_json")
        for file in self.getFiles("aa", ".fasta"):
            fasta_filename = os.path.basename(file)
            sample_name = "_".join(fasta_filename.split("_")[:-3])
            outfile = os.path.join(pfam_dir, f"{sample_name}{add_text}.json")
            if not os.path.exists(outfile):
                self._toolBelt.runPfam(db_dir, file, outfile)
            self.pfam_grl_file = os.path.join(pfam_dir, "pfam_grl.rdata")
            self.pfam_df_file = os.path.join(pfam_dir, "pfam_df.rdata")
    
    def getAnnotations(self, trimmed_dir: str, no_plot = False, gff3 = True):
        self._toolBelt.getAnnotations(self.getFolder("pfam"), self.ORF_file, gff3)
        if not no_plot:
            self.addFolder("trimmed", trimmed_dir)
            self.backMap()
            plot_dir = self.extendFolder("ORFs", "plots", "ORF_plots")
            self._toolBelt.plotAnnotations(self.pfam_grl_file, self.pfam_df_file, 
                                           plot_dir, self.getFolder("bedgraph"))

    def backMap(self, threads = 23, mapq = 0, flag = 2304, out_dir = ""):
        backmap_dir = self.extendFolder("contigs", "backmap", "backmap") if not out_dir else self.getFolder("out")
        sorted_files = {}
        for file in self.getFiles("contigs", ".fasta"):
            sample_name = getSampleName(file, extend = self.extend)
            out_file = os.path.join(backmap_dir, f"{sample_name}.bam")
            if not os.path.exists(out_file):
                current_trimmed = self.findTrimmed(sample_name)
                sam_file = os.path.splitext(out_file)[0] + ".sam"
                runBwa(file, current_trimmed, sam_file, threads)
                bam_file = os.path.splitext(out_file)[0] + ".bam"
                samSort(sam_file, bam_file)
                Cleanup(backmap_dir, ".sam")
            sorted_files.update({sample_name : out_file})
        #Probably set flags to keep or not keep them.
        Cleanup([self.getFolder("contigs")], [".64", ".pac", ".fai", 
                                              ".ann", ".amb", ".0123"])    
        out_dir = self.extendFolder("pfam", "bedgraph", "bedGraph")
        for sample, file in sorted_files.items():
            out_file = os.path.join(out_dir, f"{sample}.bedGraph")
            if not os.path.exists(f"{out_file}.gz"):
                self.bamToBG(out_file, file)
        
    @staticmethod
    def bamToBG(out_file: str, bam: str):
        runBedtools(out_file, bam)
        
    def findTrimmed(self, sample: str) -> list[str]:
        return [file for file in 
                self.getFiles("trimmed") 
                if os.path.basename(file).startswith(sample)]

class rmaHandler(blastParser):
    def blast2Rma(self, db: str, blast_kind = "BlastN"):
        output = self.addFolder("megan", self.getFolder("out"))
        for file in self.getFiles("xml", ".xml"):
            sample_name = getSampleName(file, extend = self.extend)
            contig_tools = self._toolBelt.getToolsByName("fasta", sample_name)
            if contig_tools:
                contig_filename = contig_tools[0].filename
            else:
                LOG.error("No fasta/q file in contig directory.")
                sys.exit(1)
            self._toolBelt.blast2Rma(file, 
                                     self.getFolder("out"),
                                     db, 
                                     contig_filename, 
                                     blast_kind,
                                     sample_name)
    
    def getMeganReport(self):
        report_dir = self.extendFolder("out", "reports", "Reports")
        self._toolBelt.getMeganReports(out_dir = report_dir)
    
    def updateFastaInfo(self):
        self._toolBelt.mapFastaToRma()
    
    def hitsToCSV(self, header: list):
        super().hitsToCSV(tool_kind = "rma",
                          alt_header = header)
                          
    def findRmas(self, db: str, blast_kind: str):
        for file in self.getFiles("megan", ".rma6"):
            contig_tools = self._toolBelt.getToolsByName("fasta", sample_name)
            contig_filename = contig_tools[0].filename
            sample_name = getSampleName(file, extend = self.extend)
            self._toolBelt.addRmaTool(file,
                                      self.getFolder("out"),
                                      db,
                                      contig_filename,
                                      blast_kind,
                                      sample_name)
        
class spadesTidy(fileHandler):
    def spadesToDir(self, out_dir: str, cleanup = False):
        in_dir = self.getFolder("in")
        self.addFolder("out", out_dir)
        for item in os.listdir(in_dir):
            item_abs = os.path.join(in_dir, item)
            if os.path.isdir(item_abs):
                LOG.debug(f"Found {item}")
                scaffolds = os.path.join(item_abs, "scaffolds.fasta")
                if os.path.exists(scaffolds):
                    LOG.debug(f"Found {scaffolds}.")
                    new_fasta_name = f"{os.path.basename(os.path.dirname(scaffolds))}.fasta"
                    new_fasta_file = os.path.join(self.getFolder("out"), 
                                                  new_fasta_name)
                    self.addFastaFile(scaffolds)
                    self._toolBelt.migrateFasta(scaffolds, new_fasta_file)
                if cleanup:
                    shutil.rmtree(item_abs)
                    
class dbMaker(fileHandler):
    def fetchEntrezFastas(self, id_list: list[str], email: str, api):
        out_dir = self.getFolder("fastas")
        self.fasta = os.path.join(out_dir, "ICTV_db_fastas.fasta" )
        handle = fetchEntrez(id_list, email = email, api = api)
        if handle:
            sequences = SeqIO.parse(handle, "fasta")
            with open(self.fasta, "w+") as fa:
                #SeqIO returns the count when it works, which is handy.
                count = SeqIO.write(sequences, fa, "fasta")
                LOG.info(f"{count} sequences found and written from ICTV db.")
                self.addFastaFile(self.fasta)
    
    def makeBlastDb(self, db_name = "vir"):
        runMakeblastdb(self.fasta, db_name)