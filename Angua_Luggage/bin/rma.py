# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 15:22:54 2022

@author: mwodring
"""

from .FastaKit import seqHandler
import sys

def processRma(fasta: str, rma: str, out: str, sortby: str):
    handler = seqHandler(fasta = fasta)
    #Build this into seqHandler's init if.
    handler.addRma(rma, sortby)
    handler.rmaToFasta(out)

#Probably best to keep track and mkdir if need be to hold the outputted fastas?
if __name__ == '__main__':
    _, fasta, rma, out, sortby = sys.argv
    sys.exit(processRma(fasta, rma, out, sortby))