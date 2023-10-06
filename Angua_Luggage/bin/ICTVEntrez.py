# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 17:35:43 2022

Modified from get_virus and others by Sam McGreig.
"""

import argparse
import openpyxl as opxl
from ..LuggageInterface import dbMaker
            
def getVir(filename, options):
    def checkOptions():
        for test in checklist:
            opt = test["opt"]
            to_search = [opt] if isinstance(test["opt"], str) else opt
            to_match = currentVirusDict[test["key"]]
            if not to_match:
                return
            if any(term in to_match for term in to_search):
                pass
            else:
                return
        return currentVirusDict
            
    checklist = [{"opt" : [options.host],
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
    
    allVirusDicts = []
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
                            "Exemplar" : values[17],
                            "Isolate designation" : values[20],
                            "GENBANK accession" : values[21],
                            "Genome coverage" : values[22],
                            "Genome composition" : values[24],
                            "Host" : values[25]}
        has_terms = checkOptions()
        if checkOptions():
            allVirusDicts.append(has_terms)
    print(f"Searching for {len(allVirusDicts)} entries.")
    return allVirusDicts 
        
def parseArguments():
    parser = argparse.ArgumentParser(description = "Fetches a list of viruses from the ICTV VMR file.")
    parser.add_argument("input",
                        help = "Input folder containing .xlsx files. Required.")
    parser.add_argument("email", 
                        help = "Entrez email.")
    parser.add_argument("-a", "--api",
                        help = "api_key for Entrez email. Allows 10 queries per second instead of 3")
    parser.add_argument("-g", "--genus",
                        nargs = "+",
                        help = "Restricts db to a genus or genera.")
    parser.add_argument("-f", "--family",
                        nargs = "+",
                        help = "Restricts db to a family or families.")
    parser.add_argument("-n", "--nuc",
                        help = "Restricts db to a nucleotide type, or types. Baltimore classification.",
                        nargs = "+",
                        choices = ["dsdna", 
                                   "ssrna+", "ssrnam", "ssrna", 
                                   "ssdna+", "ssdnm", "ssdna", "ssdna+m"])
    parser.add_argument("-ho", "--host",
                        help = "Restricts db to a host type or types. Default plant.",
                        nargs = "+",
                        choices = ["plants", 
                                   "algae", 
                                   "fungi", 
                                   "archaea",
                                   "vertebrates",
                                   "bacteria"],
                        #Finish filling this out.
                        default = "plants")
    parser.add_argument("-nodb", "--noblastdb",
                        help = "Do not construct Blastdb from nucleotide fasta.",
                        action = "store_true")
    parser.add_argument("--dbname",
                        help = "Name of the resulting database.",
                        default = "vir")
    #Add toggle for exemplar or not. Store_true and exmplar = E etc.
    return parser.parse_args()

def main():
    options = parseArguments()

    nucdict = {"dsdna" : "dsDNA",
               "ssrna+" : "ssRNA(+)",
               "ssrnam" : "ssRNA(-)",
               "ssrna" : "ssRNA",
               "ssdna+" : "ssDNA(+)",
               "ssdnam" : "ssDNA(-)",
               "ssdna" : "ssDNA",
               "ssdna+m" : "ssDNA(+/-)"}

    dbm = dbMaker("ICTV_db", options.input)
    dbs = []
    
    for file in dbm.getFiles("ICTV_db", ".xlsx"):
        dbs.append(getVir(file, options))
    
    viruses = [virus for viruses in dbs for virus in viruses]
    id_list = [virus["GENBANK accession"] for virus in viruses]
    dbm.extendFolder("ICTV_db", "fastas", "ICTV_db")
    dbm.fetchEntrezFastas(id_list = id_list, email = options.email, 
                          api = options.api)

    if not options.noblastdb:
        dbm.makeBlastDb(options.dbname)
        
if __name__ == '__main__':
	sys.exit(main())