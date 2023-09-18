from .Luggage import fileHandler, csvHandler
from .utils import getSampleName, Cleanup, subSeqName
from .exec_utils import *
import json, importlib.resources
from . import data
import os

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
    
    def mergeCSVOutput(self):
        out_csv = csvHandler(header)
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
                                  
    def hitsToCSV(self, add_text: ""):
        csv_out_folder = os.path.join(self.getFolder("out"), "csv")
        self.addFolder("csv", csv_out_folder)
        for filename, info in self._toolBelt.getHitsCSVInfo():
            if info:
                sample_name = getSampleName(filename, self.extend)
                out_file = os.path.join(csv_out_folder, 
                                        f"{sample_name}_{add_text}.textsearch.csv")
                csvHandler.outputHitsCSV(header = header, 
                                         rows = info, out_file = out_file)
            else:
                LOG.info(f"No suitable hits for {filename}.")
    
    def hitContigsToFasta(self, by_species = False):
        out_dir = os.path.join(self.getFolder("out"), "contigs")
        self.updateFastaInfo()
        self.addFolder("parsed_contigs", out_dir)
        if not by_species:
            self._toolBelt.outputContigsAll(out_dir = out_dir)
        else:
            self._toolBelt.outputContigBySpecies(out_dir = out_dir)
    
    def hitAccessionsToFasta(self, email: str, db_type="N"):
        db_type = "nuc" if db_type == "N" else "prot"
        out_dir = self.extendFolder("out", "acc", "hit_fastas")
        for file in self.getFiles("csv", ".textsearch.csv"):
            sample_name = getSampleName(file, self.extend)
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
    
    def runBwaTS(self, raw_dir: str, in_dir_type: str, extend = 0) -> list:
        self.addFolder("raw", raw_dir)
        self.findFastaFiles("raw")
        self.addFolder("bwa", os.path.join(self.getFolder("out"), "bwa"))
        tsv_files = []
        all_samples_dict = {}
        tmp_dir = os.path.join(self.getFolder("acc"), "tmp")
        self.addFolder("tmp", tmp_dir)
        for fasta in self.getFiles(in_dir_type, ".fasta"):
            self.addFastaFile(fasta)
            sample_name = getSampleName(fasta, extend = extend)
            seq_names, tmp_fas = self.makeTempFastas(sample_name, fasta)
            all_samples_dict[sample_name] = dict(zip(seq_names, tmp_fas))
        for sample_name, seq_to_tmp in all_samples_dict.items():
            bwa_reads = self.findRawBySample(sample_name)
            for seq_name, tmp_fa in seq_to_tmp.items():
                underscore_seq_name = subSeqName(seq_name)
                out_file = os.path.join(self.getFolder("bwa"), 
                                        f"{sample_name}_{underscore_seq_name}.sam")
                sorted_file = os.path.join(self.getFolder("bwa"), 
                                           f"{sample_name}_{underscore_seq_name}.sorted.bam")
                if not os.path.exists(sorted_file):
                    runBwa(tmp_fa, bwa_reads, out_file)
                    samSort(out_file, sorted_file)
                    self.extendFolder("out", "hist", "Histograms")
                    hist_file = os.path.join(self.getFolder("hist"), 
                                                            f"{sample_name}_{underscore_seq_name}_hist.txt")
                    outputSamHist(sorted_file, hist_file)
                    self.coverageToTSV(out_file, sample_name, seq_name)
                #Just remove the folder - write a remove dir function.
                os.remove(tmp_fa)
        Cleanup(self.getFolder("acc"), [".amb", ".ann", ".bwt", ".pac", ".sa"])
        Cleanup(self.getFolder("bwa"), [".sam"])
        return tsv_files
        
    def findRawBySample(self, sample_name: str):
        files = []
        for file in self.getFiles("raw", [".fasta", ".fq"]):
            if sample_name in file:
                files.append(file)
        return files
            
    @staticmethod
    def coverageToTSV(bwa_file: str, 
                      sample_name: str, seq_name: str) -> str:
        bwa_dir = os.path.dirname(bwa_file)
        num_mapped_reads = getNumMappedReads(bwa_file)
        tsv_file = os.path.join(bwa_dir, 
                                f"{sample_name}_{seq_name}.tsv")
        csvHandler.mappedReadsTSV(tsv_file, sample_name, seq_name, num_mapped_reads)
        return tsv_file
        
class SRA(fileHandler):
    def fetchSRAList(self, output_dir: str, SRA_file: str):
        results = []
        self.addFolder(add_dir = out_dir, dir_kind = "raw")
        with open(SRA_file, "r") as accessions:
            for accession in accessions.readlines():
                result = fetchSRA(out_dir, accession)
                results.append(result)
                
    def renameSRA():
        for file in self.getFiles(dir_kind = "raw", 
                                  file_end = "_1.fastq.gz"):
            sample_num = 1
            #This will need testing!
            whole_sample_name = getSampleName(file)
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
    
    def backMap(self):
        backmap_dir = self.extendFolder("contigs", "backmap", "backmap")
        sorted_files = {}
        for file in self.getFiles("contigs", ".fasta"):
            sample_name = getSampleName(file, extend = self.extend)
            out_file = os.path.join(backmap_dir, f"{sample_name}.bam")
            if not os.path.exists(out_file):
                current_trimmed = self.findTrimmed(sample_name)
                sam_file = os.path.splitext(out_file)[0] + ".sam"
                bwa_proc = runBwa(file, current_trimmed, sam_file)
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

class rmaHandler(fileHandler):
    def blast2Rma(self, outdir, db, reads, blast_kind = "BlastN"):
        output = self.addFolder("megan", outdir)
        for file in self.getFolder("xml", ".xml"):
            self.ToolBelt.runBlast2Rma(file, output, db, reads, blast_kind = blast_kind)
    
    def getMeganReport(self):
        self.ToolBelt.getMeganReports()

class spadesTidy(fileHandler):
    def spadesToDir(self, out_dir: str):
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