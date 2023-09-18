# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 10:19:35 2022

@author: mwodring
"""

from .FastaKit import SeqHandler

orfs_test = seqHandler(folder = "/home/mwodring/testkit/Novelstest", folder_type = "contigs")
orfs_test.runPfam()
