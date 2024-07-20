#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
java -jar trimmomatic-0.39.jar PE \
    input_forward.fq.gz input_reverse.fq.gz \
    output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
    output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

java -jar trimmomatic-0.35.jar SE \
    -phred33 input.fq.gz output.fq.gz \
    ILLUMINACLIP:TruSeq3-SE:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
"""

import os
from subprocess import check_call

from tqdm import tqdm

from db.models import RawData


TRIMMOMATIC = f"java -jar trimmomatic-0.39.jar"

PARAMS = {
    "PE": "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36",
    "SE": "ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
}


def main():
    run_accessions = {}
    for rec in RawData.select(RawData.run_accession, RawData.path):
        if os.path.exists(RawData.path):
            if os.path.exists(RawData.path.replace(".fastq.gz", ".clean.fastq.gz")):
                if rec.run_accession not in run_accessions:
                    run_accessions[rec.run_accession] = []
                run_accessions[rec.run_accession].append(rec.path)

    for accession, files in tqdm(run_accessions.items()):
        files = sorted(files)

        fs = []
        for f in files:
            fs.append(f.replace(".fastq.gz", ".clean.fastq.gz"))
            fs.append(f.replace(".fastq.gz", ".unpaired.fastq.gz"))
        if len(files) > 1:
            cmd = f"{TRIMMOMATIC} {' '.join(files)} {' '.join(fs)} {PARAMS['PE']}"
        else:
            cmd = f"{TRIMMOMATIC} {' '.join(files)} {' '.join(fs)} {PARAMS['SE']}"

        check_call(cmd, shell=True)



if __name__ == '__main__':
    pass
