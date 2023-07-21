from .Luggage import fileHandler, csvHandler
from .utils import getSampleName
from .exec_utils import *
import json, importlib.resources
from . import data

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
        self.updateFastaInfo()
        csv_out_folder = os.path.join(self.getFolder("out"), "csv")
        self.addFolder("csv", csv_out_folder)
        for filename, info in self._toolBelt.getHitsCSVInfo():
            if info:
                sample_name = getSampleName(filename)
                out_file = os.path.join(csv_out_folder, 
                                        f"{sample_name}_{add_text}.textsearch.csv")
                csvHandler.outputHitsCSV(header = header, 
                                         rows = info, out_file = out_file)
            else:
                LOG.info(f"No suitable hits for {filename}.")
    
    def hitContigsToFasta(self, by_species = False):
        out_dir = os.path.join(self.getFolder("out"), "contigs")
        self.addFolder("parsed_contigs", out_dir)
        if not by_species:
            self._toolBelt.outputContigsAll(out_dir = out_dir)
        else:
            self._toolBelt.outputContigBySpecies(out_dir = out_dir)
    
    def hitAccessionsToFasta(self, email: str):
        hit_dir = os.path.join(self.getFolder("out"), "hit_fastas")
        self.addFolder("accs", hit_dir)
        #TODO: This should use the merged csv or better yet the toolBelt.
        for file in self.getFiles("csv", ".textsearch.csv"):
            filename = os.path.splitext(os.path.basename(file))[0]
            sample_name = "_".join(filename.split("_")[:-3])
            accessions = csvHandler.getCSVAccessions(file)
            fa_filename = f"{os.path.splitext(sample_name)[0]}_accessions.fasta"
            fetchEntrezFastas(id_list = accessions, 
                                        out_dir = hit_dir, 
                                        email = email, 
                                        filename = fa_filename)
            self.addFastaFile(fa_filename)
        #Might be worth naming these based on the samples they come from rather 
        #than the full csv.
        
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
        self.addFolder("bwa", os.path.join(self.getFolder("out"), "bwa"))
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
                                        f"{sample_name}_{seq_name}.sorted.bam")
                runBwa(tmp_fa, bwa_reads, out_file)
                #Just remove the folder - write a remove dir function.
                os.remove(tmp_fa)
                self.coverageToTSV(out_file, sample_name, seq_name)
        #Presumably needs the .clean or whatever?
        Cleanup([self._dirs["acc"]], [".amb", ".ann", ".bwt", ".pac", ".sa"])
        return tsv_files
    
    def coverageToTSV(self, bwa_file: str, 
                       sample_name: str, seq_name: str) -> str:
        num_mapped_reads = getNumMappedReads(bwa_file, 
                                                        sample_name, seq_name, 
                                                        self.getFolder("bwa"))
        tsv_file = os.join(self.getFolder('bwa'), 
                           f"{sample_name}_{seq_name}.tsv")
        csvHandler.mappedReadsTSV(tsv_file, sample_name)
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
    