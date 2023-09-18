# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 18:40:02 2022

@author: mwodring
"""

from .FastaKit import seqHandler
import sys

out_file = sys.argv[1]
sra_file = sys.argv[2]

#For now, just init an empty seqHandler and run the SRA. I could have a folder type it checks for and then runs fetchSRA by default but that would be less elegant imo. 
#This way 'raw' reads can be same across other uses, such as tracking something through Angua if desired.
#That's why it's in FastaKit and not on its own. 
handler = seqHandler()
handler.fetchSRA(output_folder = out_file, SRA_file = sra_file)
sample_number = handler.renameSRA()

print("Samples pairs downloaded:", sample_number)