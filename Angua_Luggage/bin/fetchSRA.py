# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 18:40:02 2022

@author: mwodring
"""
from ..LuggageInterface import SRA
import sys, os, logging

logging.basicConfig(stream = sys.stdout, level=logging.DEBUG)
LOG = logging.getLogger(__name__)

def main():
    out_file = sys.argv[1]
    sra_file = sys.argv[2]

    handler = SRA("raw", os.path.basename(out_file), 0)
    count = handler.fetchSRAList(sra_file)
    
    handler.pigzAll()

    LOG.info(f"Sample (or pairs) downloaded: {count}.")
    
if __name__ == "__main__":
    sys.exit(main())