# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 12:52:31 2022

@author: mwodring
"""

import os, pathlib, re
from collections import namedtuple

SearchParams = namedtuple("Search_Params", 
                          "search_term,minlen,bitscore,blacklist")

def subSeqName(to_sub: str):
    seq_name = re.sub("[ :.,;()/\>|]", "_", to_sub)
    return seq_name

def Cleanup(folders: list | str, filetypes: list):
    def cleanup(folder):
        for file in os.scandir(folder):
            last_suffix = (lambda suffixes : suffixes[-1] if len(suffixes) > 0 else "None")(pathlib.Path(file.name).suffixes)
            if last_suffix in filetypes:
                os.remove(file)
    if type(folders) == list:
        for folder in folders:
            cleanup(folder)
    else:
        cleanup(folders)
        
def getpercentage(a, b):
    return(a / b) * 100
    
#This only works for a specific format atm. How to go about this?  
def getSampleName(file: str, extend = None):
    sample = os.path.splitext(os.path.basename(file))[0]
    sample = sample.split(".")[0]
    if extend:
        sample = "_".join(sample.split("_")[:-extend])
    return sample
    
def count_calls():
    pass

