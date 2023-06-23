# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 12:52:31 2022

@author: mwodring
"""

import logging 
import sys
import os
import functools
import pathlib

#This is just magpied from the 'net. Decorators are magic. https://realpython.com/primer-on-python-decorators/
def count_calls(func):
    @functools.wraps(func)
    def wrapper_count_calls(*args, **kwargs):
        wrapper_count_calls.num_calls += 1
        return func(*args, **kwargs)
    wrapper_count_calls.num_calls = 0
    return wrapper_count_calls

#Consider a logger class?
msg = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

def new_file_handler(output_dir, file):
    log_file = os.path.join(output_dir, f"{file}.log")
    file_handler = logging.FileHandler(log_file)
    file_format = logging.Formatter(msg)
    file_handler.setFormatter(file_format)
    return file_handler

def new_logger(name = __name__, output_dir = "", file = "logs"):
    logger = logging.getLogger(name)
    
    #If I don't do this, it repeats itself for some reason.
    if not logger.handlers:
        #logger.propagate = False
        logger.setLevel(logging.DEBUG)
        
        console_handler = logging.StreamHandler(sys.stdout)
        console_format = logging.Formatter(msg)
        console_handler.setFormatter(console_format)
        logger.addHandler(console_handler)
        
        if output_dir != "":
            logger.addHandler(new_file_handler(output_dir, file))
            
    return logger

def add_logger_outdir(output_dir, file, logger):
    logger.addHandler(new_file_handler(output_dir, file))

def Cleanup(folders: list, filetypes: list):
    for folder in folders:
        for file in os.scandir(folder):
            last_suffix = (lambda suffixes : suffixes[-1] if len(suffixes) > 0 else "None")(pathlib.Path(file.name).suffixes)
            if last_suffix in filetypes:
                os.remove(file)