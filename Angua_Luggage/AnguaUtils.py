# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 19:08:04 2023

@author: smcgreig
"""

from pathlib import Path
import os
import subprocess
import pysam

################################################################################
def looper(input_dir, qualifier = "N.A"):
    # Loop through files in a directory and return a list of contents

    for file in os.listdir(input_dir):
        if qualifier != "N.A":
            if(qualifier in file):
                yield file
        elif qualifier == "N.A":
            yield file 


################################################################################
def back_mapper(input_dir, output_dir, reference, threads, delim, mapq, flag, coverage, min_alignment_length, softclip_perc, logger):

    with open(reference) as ref_fasta:
        seq_count = 0
        for line in ref_fasta:
            if(">" in line):
                seq_count += 1

        if(seq_count > 1):
            logger.warning('''Multiple sequences detected in the reference fasta file. Please ensure that each reference file represents ONE DISTINCT genome.
            If multiple genomes are used in a single reference file, then reads that align equally well will ONLY be assigned to a single reference.
            This can and will skew the mapping results.
            Consider using the --flag 0 option to correct this.''')
            


    input_files = looper(input_dir, qualifier = "_R1")

    # Index reference
    # subprocess.call(f"bwa-mem2 index {reference}", shell = True)

    for file in input_files:
        sample_ID = file.split(delim)[0]
        Path(output_dir, sample_ID).mkdir(exist_ok = True, parents = True)
        fileR2 = file.replace("_R1", "_R2")

        # Index reference
        subprocess.run(["bwa-mem2", "index", reference])

        if os.path.isfile(f"{input_dir}/{file}"):
            r1_in = f"{input_dir}/{file}"
            r2_in = f"{input_dir}/{fileR2}"
            sample_out = f"{output_dir}/{sample_ID}/"

            # Map raw reads
            # subprocess.call(f"bwa-mem2 mem -t {threads} {reference} {r1_in} {r2_in} > {sample_out}{sample_ID}.sam", shell = True)
            # subprocess.call(f"samtools view -q {mapq} -F {flag} -bS {sample_out}/{sample_ID}.sam > {sample_out}{sample_ID}.bam", shell = True)
            # subprocess.call(f"samtools sort {sample_out}/{sample_ID}.bam > {sample_out}{sample_ID}_sort.bam", shell = True)
            
            #Another one I'm unsure how to work out so I'll leave it be for now!
            subprocess.call(f"bwa-mem2 mem -t {threads} {reference} {r1_in} {r2_in} | samtools view -q {mapq} -F {flag} -bS - | samtools sort - > {sample_out}{sample_ID}_sort.bam", shell = True)
            subprocess.run(["samtools", "index", f"{sample_out}{sample_ID}_sort.bam"])

            # Pysam filtering
            samfile = pysam.AlignmentFile(f"{sample_out}{sample_ID}_sort.bam", "rb")
            samfile_out = pysam.AlignmentFile(f"{sample_out}{sample_ID}_sort_pysam.bam", "wb", template = samfile)

            sam_iterator = samfile.fetch()
            for alignment in sam_iterator:
                if int(alignment.query_alignment_length) >= int(min_alignment_length):
                    if (float(alignment.query_length) * float(softclip_perc)) >= (float(alignment.query_length) - float(alignment.query_alignment_length)):
                        samfile_out.write(alignment)

            subprocess.call(f"samtools index {sample_out}{sample_ID}_sort_pysam.bam", shell = True)
            subprocess.call(f"samtools idxstats {sample_out}{sample_ID}_sort_pysam.bam > {sample_out}{sample_ID}_stats.txt", shell = True)
            # Calculate Coverage
            if coverage == "Y":
                subprocess.run(["average-coverage.py", f"{sample_out}{sample_ID}_sort_pysam.bam", "-o", f"{sample_out}{sample_ID}_coverage.tsv"])
                logger.info(f"average-coverage.py {sample_out}{sample_ID}_sort_pysam.bam -o {sample_out}{sample_ID}_coverage.tsv")

    logger.info("Back mapping completed.")
