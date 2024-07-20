#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import gzip, os, time
from glob import glob
from multiprocessing import Pool, Queue, Process
from subprocess import check_call

import pandas as pd
from tqdm import tqdm

from db.models import RMats, Ref, Meta


"""
# rename id
python /mnt/raid62/zhouran/project/2024-database/splice_test/external/addchr.py Mus_musculus.GRCm38.dna.primary_assembly.fa Mus_musculus.GRCm38.dna.primary_assembly.chr.fa fasta
python /mnt/raid62/zhouran/project/2024-database/splice_test/external/addchr.py Mus_musculus.GRCm38.101.filtered.gtf Mus_musculus.GRCm38.101.filtered.chr.gtf gtf
python /mnt/raid62/zhouran/project/2024-database/splice_test/external/gtf2bed.py Mus_musculus.GRCm38.101.filtered.chr.gtf > Mus_musculus.GRCm38.101.filtered.chr.bed12

# bed12 must start with 'chr'

tool=/mnt/raid62/zhouran/project/2024-database/splice_test/tool/SpliceTools-main/bin

## intron retention

### intron exon size
# 这部分拿到的结果都是一个整体的统计分布图，必须要把所有的结构加上是哪个剪接事件可能才好用

最后拿到的结果应该是1_intron_exon_size_SUMMARY.tsv，得到上调和下调内含子的长度分布
然后根据2_intron_exon_sizes_neg_IncDiff.tsv和2_intron_exon_sizes_neg_IncDiff.tsv的里面的图来画boxplot
考虑改成pvalue不？

input=/mnt/raid62/zhouran/project/2024-database/splice_test/Perturb_1932/RI.MATS.JC.txt
bed12=/mnt/raid61/Ref/MusMus/release101/Mus_musculus.GRCm38.101.filtered.chr.bed12
genome=/mnt/raid61/Ref/MusMus/release101/Mus_musculus.GRCm38.dna.primary_assembly.chr.fa
perl $tool/RIIntronExonSizes.pl \
 -r $input \
 -a $bed12 -f 2

### splice site score

最后拿到的结果是
1_neg_IncDiff/out.downstream_acceptor_neg.scores.txt, 1_neg_IncDiff/out.upstream_donor_neg.scores.txt
2_pos_IncDiff/out.downstream_acceptor_pos.scores.txt, 2_pos_IncDiff/out.upstream_donor_pos.scores.txt

perl $tool/RISpliceSiteScoring.pl \
 -r $input \
 -a $bed12 \
 -g $genome -f 2

## SE
"""


__DIR__ = os.path.dirname(os.path.abspath(__file__))


def modify_fasta(input_file, output_file):
    infile = open(input_file, "r") if input_file.endswith("fa") else gzip.open(input_file, "rt")
    with open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                if line.startswith('>MT'):
                    line = '>chrM' + line[3:]
                elif not line.startswith('>chr'):
                    line = '>chr' + line[1:]
            outfile.write(line)
    infile.close()


def modify_gtf(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            if len(fields) >= 9:
                if fields[0] == 'MT':
                    fields[0] = 'chrM'
                elif not fields[0].startswith('chr'):
                    fields[0] = 'chr' + fields[0]
                outfile.write('\t'.join(fields) + '\n')
            else:
                outfile.write(line)


def parse_gtf_to_bed12(gtf_file, bed_file):
    transcript_info = {}

    file = open(gtf_file, 'r') if gtf_file.endswith(".gtf") else gzip.open(gtf_file, 'rt')

    for line in file:
        if line.startswith('#'):
            continue  # Skip header lines
        fields = line.strip().split('\t')
        feature_type = fields[2]
        chrom = fields[0]

        if chrom == "MT":
            chrom = "chrM"
        elif not chrom.startswith('chr'):
            chrom = "chr" + chrom

        start = int(fields[3]) - 1  # BED format is 0-based
        end = int(fields[4])
        strand = fields[6]

        # Parse attributes
        attributes = fields[8]
        gene_id = ''
        transcript_id = ''
        gene_name = ''

        for attribute in attributes.split(';'):
            if 'gene_id' in attribute:
                gene_id = attribute.split('"')[1]
            if 'transcript_id' in attribute:
                transcript_id = attribute.split('"')[1]
            if 'gene_name' in attribute:
                gene_name = attribute.split('"')[1]

        # Initialize transcript_info entry if new
        if transcript_id not in transcript_info:
            transcript_info[transcript_id] = {
                'chrom': chrom,
                'start': start,
                'end': end,
                'name': f"{transcript_id}_{gene_name}",
                'score': '1000',
                'strand': strand,
                'thickStart': float('inf'),  # Initialize to infinities for comparison
                'thickEnd': -1,
                'itemRgb': '0',
                'blockCount': 0,
                'blockSizes': [],
                'blockStarts': [],
            }

        entry = transcript_info[transcript_id]
        if feature_type == "exon":
            entry['start'] = min(entry['start'], start)
            entry['end'] = max(entry['end'], end)
            entry['blockSizes'].append(end - start)
            entry['blockStarts'].append(start - entry['start'])
            entry['blockCount'] += 1
        elif feature_type == "start_codon" and (strand == '+' or strand == '-'):
            entry['thickStart'] = min(entry['thickStart'], start) if strand == '+' else max(entry['thickStart'], start)
        elif feature_type == "stop_codon" and (strand == '+' or strand == '-'):
            entry['thickEnd'] = max(entry['thickEnd'], end) if strand == '+' else min(entry['thickEnd'], end)

    file.close()

    # Adjust thickStart and thickEnd if no codon was found
    for transcript_id, info in transcript_info.items():
        if info['thickStart'] == float('inf'):  # If no start_codon found
            info['thickStart'] = info['end']
        if info['thickEnd'] == -1:  # If no stop_codon found
            info['thickEnd'] = info['end']

    # Write the BED12 file
    with open(bed_file, 'w') as outfile:
        for transcript_id, info in transcript_info.items():
            bed12_fields = [
                info['chrom'],
                info['start'],
                info['end'],
                info['name'],
                info['score'],
                info['strand'],
                info['thickStart'],
                info['thickEnd'],
                info['itemRgb'],
                info['blockCount'],
                ','.join(map(str, info['blockSizes'])) + ',',
                ','.join(map(str, info['blockStarts'])) + ','
            ]
            outfile.write('\t'.join(map(str, bed12_fields)) + '\n')


def generate_expr(path):
    conditions = ["b1.txt", "b2.txt"]
    condition_group = []

    data = None
    for group, condition in enumerate(conditions):
        with open(os.path.join(path, condition)) as r:
            for line in r:
                line = line.strip().split(",")

                for l in line:
                    condition_group.append(group)

                    for f in glob(l.replace("Aligned.sortedByCoord.out.bam", "ReadsPerGene.out.tab*")):
                        try:
                            df = pd.read_csv(f, sep="\t", skiprows=4, header=None)
                        except Exception:
                            continue
                        df = df.iloc[:, 0:2]
                        df.columns = ["gene", os.path.basename(f).split(".")[0]]

                        if data is None:
                            data = df
                        else:
                            data = data.merge(df, how="inner", on="gene")
    if data is not None:
        data.to_csv(os.path.join(path, "expression.tsv"), sep="\t", index=False)
        return condition_group


def decode(line):
    res = {}
    for row in line.split(";"):
        row = row.split()
        if len(row) > 1:
            res[row[0].strip()] = row[1].replace("\"", "").replace(";", "").strip()
    return res


def load_ref():
    path = os.path.join(__DIR__, "bed12")
    os.makedirs(path, exist_ok=True)

    references = {}
    genomes = {}
    for ref in Ref.select():
        if not os.path.exists(os.path.join(str(path), f"{ref.sci_name}.bed12")):
            print("generate bed12 for %s" % ref.sci_name)

            for f in glob(os.path.join(ref.path, "*.gtf")) + glob(os.path.join(ref.path, "*.gtf.gz")):
                if "filter" not in f and "chr" not in f and "sorted" not in f:
                    parse_gtf_to_bed12(f, os.path.join(str(path), f"{ref.sci_name}.bed12"))
                    break

        if not os.path.exists(os.path.join(str(path), f"{ref.sci_name}.fa")):
            print("generate fasta for %s" % ref.sci_name)

            for f in glob(os.path.join(ref.path, "*.fa")) + glob(os.path.join(ref.path, "*.fa.gz")):
                modify_fasta(f, os.path.join(str(path), f"{ref.sci_name}.fa"))
                break

        references[ref.sci_name] = os.path.join(str(path), f"{ref.sci_name}.bed12")
        genomes[ref.sci_name] = os.path.join(str(path), f"{ref.sci_name}.fa")

    return references, genomes


__HOME__ = "/mnt/raid62/zhouran/project/2024-database/splice_test/tool/SpliceTools-main/bin"


def generate_cmds(path: str, ref: str, genome:str):
    # groups = generate_expr(path)
    # if groups is None:
    #     return []
    #
    # tpm = ",".join(["2" for _ in groups])
    #
    # temp_groups = {}
    # for i in groups:
    #     temp_groups[i] = 1 + temp_groups.get(i, 0)
    # sn = ",".join([str(temp_groups[x]) for x in groups])

    os.makedirs(os.path.join(path, "temp"), exist_ok=True)
    cmds = []
    # for i in glob(os.path.join(path, "*.JCEC.txt")):
    #     if not os.path.exists(os.path.join(path, "temp", os.path.basename(i))):
    #         os.symlink(i, os.path.join(path, "temp", os.path.basename(i)))

    # if not glob(os.path.join(path, "SpliceCompare_out_FDR_0.05")):
    for j in ["SE"]: # "RI",
        i = os.path.join(path, f"{j}.MATS.JCEC.txt")
        if os.path.exists(os.path.join(path, f"{j}.MATS.JCEC.txt")):
            param = "-r" if j == "RI" else "-s"

            cmds += [
                # f"perl {__HOME__}/{j}FractionExpressed.pl {param} {i} -e {os.path.join(path, 'expression.tsv')} -TPM {tpm} -SN {sn} -f 0.05 > {path}/{j}SpliceSiteScoring.log",
                # f"perl {__HOME__}/{j}IntronExonSizes.pl {param} {i} -a {ref} -f 0.05 > {path}/{j}SpliceSiteScoring.log",
                # f"perl {__HOME__}/{j}SpliceSiteScoring.pl {param} {i} -g {genome} -f 0.05 -a {ref} > {path}/{j}SpliceSiteScoring.log",
                f"perl {__HOME__}/SETranslateNMD.pl {param} {i} -g {genome} -f 0.05 -a {ref} > {path}/{j}TranslateNMD.log"
            ]

        # cmds.append(f"perl {__HOME__}/SpliceCompare.pl -i {path}/temp -o {path} -m 280 -p 1 -f 0.05 > {path}/{j}SpliceSiteScoring.log")
    return cmds


def call(cmd):
    with open(os.devnull) as w:
        try:
            check_call(cmd, shell=True,)
        except Exception:
            pass


# def calls(inputQueue: Queue, outputQueue: Queue):
#
#     while True:
#         perturb_id, cmds = inputQueue.get()
#         for cmd in cmds:
#             call(cmd)
#         RMats.update(**{"splice_tools": True}).where(RMats.perturb_id == perturb_id).execute()
#         outputQueue.put(perturb_id)


def calls(args):
    perturb_id, cmds = args
    for cmd in cmds:
        call(cmd)
    RMats.update(**{"splice_tools": True}).where(RMats.perturb_id == perturb_id).execute()


def generate_cmds___(args):
    query, sci_names, reference, genomes = args
    query["sci_name"] = sci_names[query["project_accession"]]

    path = os.path.join(
        os.path.dirname(query["path"]),
        str(os.path.basename(query["path"])).replace("|", "_")
    )

    key = query["sci_name"]
    if not os.path.exists(os.path.join(path, os.path.basename(reference[key]))):
        os.system(f"ln -sf {reference[key]} {os.path.join(path, os.path.basename(reference[key]))}")

    return [query['perturb_id'], generate_cmds(path, ref=os.path.join(path, os.path.basename(reference[key])), genome=genomes[key])]


def main(server: str=1140, n_jobs: int = 2):
    inqueue, outqueue = Queue(), Queue()
    reference, genomes = load_ref()

    # Meta.sci_name,

    queries = RMats.select().where(
        (RMats.server == str(server)) # & (RMats.splice_tools == False) & (RMats.properly == "")
    )

    sci_names = {}
    for rec in tqdm(
            Meta.select(Meta.secondary_study_accession, Meta.sci_name).distinct(),
            total=Meta.select(Meta.secondary_study_accession, Meta.sci_name).distinct().count()):
        sci_names[rec.secondary_study_accession] = rec.sci_name

    cmds0 = []
    for query in tqdm(queries.dicts(), total=queries.count()):
        summ = os.path.join(query["path"], "summary.txt")
        if not os.path.exists(summ) or os.path.getsize(summ) < 10:
            continue

        if query["project_accession"] not in sci_names:
            continue

        cmds0.append([query, sci_names, reference, genomes])

    with Pool(n_jobs) as p:
        cmds = list(tqdm(p.imap(generate_cmds___, cmds0), total=len(cmds0), desc="generating cmds"))

    with Pool(n_jobs) as p:
        list(tqdm(p.imap(calls, cmds), total=len(cmds), desc="Splice Tools"))


if __name__ == '__main__':
    # perl /mnt/raid62/zym/env/SpliceTools/bin/RIFractionExpressed.pl -r /mnt/data6/zym/rmats/PRJNA801846/OE/RI.MATS.JC.txt
    """
    perl /mnt/raid62/zym/env/SpliceTools/bin/RIMedley.pl \
        -r /mnt/data6/zym/rmats/PRJNA433339/OE/RI.MATS.JCEC.txt \
        -a /mnt/raid62/zym/project/dataexplorer/bed12/HomSap.bed12 \
        -g /mnt/raid61/Ref/HomSap/release101/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -e /mnt/data6/zym/rmats/PRJNA433339/OE/expression.tsv \
        -TPM 2,2 -SN 1,1 -f 0.05

    perl RIFractionExpressed.pl \
        -r /mnt/data6/zym/rmats/PRJNA801846/OE/RI.MATS.JCEC.txt \
        -e /mnt/data6/zym/rmats/PRJNA801846/OE/expression.tsv \
        -TPM 2,2 -SN 1,1 -f 0.05
        
    perl RIIntronExonSizes.pl \
	    -r /mnt/data6/zym/rmats/PRJNA801846/OE/RI.MATS.JCEC.txt \
	    -a /mnt/raid62/zym/project/dataexplorer/bed12/HomSap.bed12 \
	    -f 0.05
	    
    perl RISpliceSiteScoring.pl \
        -r /mnt/data6/zym/rmats/PRJNA801846/OE/RI.MATS.JCEC.txt \
        -g /mnt/raid61/Ref/HomSap/release101/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -f 0.05 \
        -a /mnt/raid62/zym/project/dataexplorer/bed12/HomSap.bed12
        
    perl /mnt/raid62/zym/env/SpliceTools/bin/SEMedley.pl \
        -r /mnt/data6/zym/rmats/PRJNA433339/OE/SE.MATS.JCEC.txt \
        -a /mnt/raid62/zym/project/dataexplorer/bed12/HomSap.bed12 \
        -g /mnt/raid61/Ref/HomSap/release101/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -e /mnt/data6/zym/rmats/PRJNA433339/OE/expression.tsv \
        -TPM 2,2 -SN 1,1 -f 0.05
       
    
    perl SEFractionExpressed.pl \
        -s /mnt/data6/zym/rmats/PRJNA801846/OE/SE.MATS.JCEC.txt \
        -e /mnt/data6/zym/rmats/PRJNA801846/OE/expression.tsv \
        -TPM 2,2 -SN 1,1 -f 0.05
        
    perl SEIntronExonSizes.pl \
	    -s /mnt/data6/zym/rmats/PRJNA801846/OE/SE.MATS.JCEC.txt \
	    -a /mnt/raid62/zym/project/dataexplorer/bed12/HomSap.bed12 \
	    -f 0.05
	    
    perl SESpliceSiteScoring.pl \
        -s /mnt/data6/zym/rmats/PRJNA801846/OE/SE.MATS.JCEC.txt \
        -g /mnt/raid61/Ref/HomSap/release101/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -f 0.05 \
        -a /mnt/raid62/zym/project/dataexplorer/bed12/HomSap.bed12
    
       
	perl /mnt/raid62/zym/env/SpliceTools/bin/SpliceCompare.pl \
	    -i /mnt/data6/zym/rmats/PRJNA801846/OE/ \
	    -o /mnt/data6/zym/rmats/PRJNA801846/OE/spliceTools -m 280 -p 1 -f 0.05
    """
    #
    from fire import Fire
    Fire(main)
