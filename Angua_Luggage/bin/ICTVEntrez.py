# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 17:35:43 2022

Modified from get_virus and others by Sam McGreig.
"""

import argparse
import openpyxl as opxl
from FastaKit import seqHandler
            
def getVir(filename):
    def checkOptions(checklist):
        try:
            for test in checklist:
                if test["opt"] == currentVirusDict[test["key"]]:
                    pass
                else:
                    return
            allVirusDict[values[0]] = currentVirusDict
        except TypeError:
            print(f"Missing taxa on {values[0]}.")
            
    checklist = [{"opt" : options.host,
                  "key" : "Host"}]
    if options.nuc:
        checklist.append({"opt" : options.nuc,
                         "key" : "Genome composition"})
    if options.family:
        checklist.append({"opt" : options.family,
                         "key" : "Family"})
    if options.genus:
        checklist.append({"opt" : options.genus,
                         "key" : "Genus"})
         
    allVirusDict = {}
    wb = opxl.load_workbook(filename)
    sheet = wb.active
    for row in sheet.iter_rows(min_row = 2):
        values = [data.value for data in row]
        currentVirusDict = {"Realm" : values[2],
                            "Kingdom" : values[4],
                            "Class" : values[8],
                            "Order" : values[10],
                            "Family" : values[12],
                            "Genus" : values[14],
                            "Species" : values[16],
                            "Exemplar" : values[17],
                            "Virus name" : values[18],
                            "Abbreviations" : values[19],
                            "Isolate designation" : values[20],
                            "GENBANK accession" : values[21],
                            "Genome coverage" : values[22],
                            "Genome composition" : values[23],
                            "Host" : values[24]}
        checkOptions(checklist)
    print(f"Searching for {len(allVirusDict)} entries.")
    return allVirusDict

def getFastas(virdict, handler):      
    id_list = [virus["GENBANK accession"] for virus in virdict.values()]
    handler.fetchEntrezFastas(id_list = id_list, email = options.email, api = options.api, output = options.output)
        
def parseArguments():
    parser = argparse.ArgumentParser(description = "Fetches a list of viruses from the ICTV formatted file.")
    parser.add_argument("input",
                        help = "Input folder containing .xlsx files. Required.")
    parser.add_argument("output",
                        help = "Output folder for the db. Required.")
    parser.add_argument("email", 
                        help = "Entrez email.")
    parser.add_argument("-a", "--api",
                        help = "api_key for Entrez email. Allows 10 queries per second instead of 3")
    parser.add_argument("-g", "--genus",
                        help = "Restricts db to a genus.")
    parser.add_argument("-f", "--family",
                        help = "Restricts db to a family.")
    parser.add_argument("-n", "--nuc",
                        help = "Restricts db to a nucleotide type, Baltimore classification.",
                        choices = ["dsdna", 
                                   "ssrna+", "ssrnam", "ssrna", 
                                   "ssdna+", "ssdnm", "ssdna", "ssdna+m)"])
    parser.add_argument("-ho", "--host",
                        help = "Restricts db to a host type. Default plant.",
                        choices = ["plants", 
                                   "algae", 
                                   "fungi", 
                                   "archaea",
                                   "vertebrates",
                                   "bacteria"],
                        #Finish filling this out.
                        default = "plants")
    parser.add_argument("-db", "--blastdb",
                        help = "Construct Blastdb from nucleotide fasta.",
                        action = "store_true")
    parser.add_argument("--dbname",
                        help = "Name of the resulting database.",
                        default = "db")
    #Add toggle for exemplar or not. Store_true and exmplar = E etc.
    return parser.parse_args()

options = parseArguments()
if not options.api:
    options.api = False

nucdict = {"dsdna" : "dsDNA",
           "ssrna+" : "ssRNA(+)",
           "ssrnam" : "ssRNA(-)",
           "ssrna" : "ssRNA",
           "ssdna+" : "ssDNA(+)",
           "ssdnam" : "ssDNA(-)",
           "ssdna" : "ssDNA",
           "ssdna+m" : "ssDNA(+/-)"}

handler = seqHandler(folder = options.input, folder_type = "ICTV_db")
toparse = handler.getFiles(file_end = ".xlsx")

virus_Dict = {}

for file in toparse:
    virus_Dict.update(getVir(file))
    
getFastas(virus_Dict, handler)

if options.blastdb:
    handler.makeBlastDb(options.output, options.dbname)