# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 12:52:31 2022

@author: mwodring
"""

import os
import functools
import pathlib
from collections import namedtuple

SearchParams = namedtuple("Search_Params", 
                          "search_term,minlen,bitscore,blacklist")

#This is just magpied from the 'net. Decorators are magic. https://realpython.com/primer-on-python-decorators/
def count_calls(func):
    @functools.wraps(func)
    def wrapper_count_calls(*args, **kwargs):
        wrapper_count_calls.num_calls += 1
        return func(*args, **kwargs)
    wrapper_count_calls.num_calls = 0
    return wrapper_count_calls

def Cleanup(folders: list, filetypes: list):
    for folder in folders:
        for file in os.scandir(folder):
            last_suffix = (lambda suffixes : suffixes[-1] if len(suffixes) > 0 else "None")(pathlib.Path(file.name).suffixes)
            if last_suffix in filetypes:
                os.remove(file)

#This only works for a specific format atm. How to go about this?  
def getSampleName(file: str):
    sample = os.path.splitext(os.path.basename(file))[0]
    sample = sample.split(".")[0]
    return sample