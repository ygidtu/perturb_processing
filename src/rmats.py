#!/usr/bin/env python3
# -*-coding:utf-8 -*-

u"""
批量运行rMATs

python rmats.py --b1 /path/to/prep2.txt --gtf /path/to/the.gtf -t paired --readLength 50 --nthread 4 --od /path/to/output --tmp /path/to/tmp_output_prep_2 --task prep
/path/to/pair_1_a.bam,/path/to/pair_2_a.bam,/path/to/pair_3_a.bam

Using the paired stats model
The default statistical model considers the samples to be unpaired. The --paired-stats flag can be used if each entry in --b1 is matched with its pair in --b2. As an example, if there are three replicates where each replicate has paired "a" and "b" data, then b1.txt and b2.txt should look like:

/path/to/b1.txt
/path/to/pair_1_a.bam,/path/to/pair_2_a.bam,/path/to/pair_3_a.bam
/path/to/b2.txt
/path/to/pair_1_b.bam,/path/to/pair_2_b.bam,/path/to/pair_3_b.bam
"""
import os, json
from datetime import datetime
from glob import glob
from db.models import Meta, StarData, Ref, RMats
from tqdm import tqdm
from subprocess import check_call
from multiprocessing import Pool, Lock
from loguru import logger

import pandas as pd
import pysam
import numpy as np

# from multiprocessing import set_start_method
# set_start_method("spawn")

pysam.set_verbosity(0)

logger.remove()
logger.add(lambda msg: tqdm.write(msg, end=""), colorize=True)


# df = pd.read_excel("merge_0523.xlsx")
# conditions = {}
# for _, row in df.iterrows():
#     conditions[row["ID"]] = row["Condition"]
#
# for id, cdt in tqdm(conditions.items()):
#     if Meta.get_or_none(Meta.experiment_accession == id):
#         Meta.update(condition=cdt).where(Meta.experiment_accession == id).execute()


__DIR__ = os.path.dirname(os.path.abspath(__file__))
__CMDS__ = os.path.join(__DIR__, 'rmats_cmds')


def get_references():
    refs = {}
    for ref in Ref.select():
        fs = glob(os.path.join(ref.path, "*.gtf"))
        if fs:
            for f in fs:
                if "chr" not in f and "sorted" not in f and "filter" not in f:
                    refs[ref.sci_name] = f
                    break
    return refs

REFERENCES = {}


# def paired_data(server: str = "1130", species: str = None):
#     project_accession = pw.TextField(index=True)
#     perturb_id = pw.TextField(index=True)
#     path = pw.TextField()
#     server = pw.CharField(index=True, null=True)
#     cmd = pw.TextField()
#     properly = pw.TextField(index=True)
#
#     query = (Meta.
#              select(Meta.study_accession, StarData.experiment_accession,
#                     StarData.path, Meta.sci_name, Ref.path.alias("ref"), StarData.cmd).
#              join(StarData, on=(Meta.experiment_accession == StarData.experiment_accession)).
#              join(Ref, on=(Meta.sci_name == Ref.sci_name)).
#              where(StarData.server == str(server)))
#
#     if species:
#         query = query.where(Meta.sci_name == species)
#
#     projs = {}
#     for data in query.dicts():
#         data["condition"] = data["condition"].replace("|", "_")
#
#         data["path"] = os.path.join(data["path"], data["experiment_accession"])
#         if data["study_accession"] not in projs.keys():
#             projs[data["study_accession"]] = {}
#
#         if data["condition"] not in projs[data["study_accession"]].keys():
#             projs[data["study_accession"]][data["condition"]] = []
#
#         projs[data["study_accession"]][data["condition"]].append(data)
#
#     return projs


def guess_strandness(path) -> str:
    for p in path:
        try:
            df = pd.read_csv(p, sep='\t', names=['gene', 'unstrand', 'strand', 'reverse']).iloc[4:]
            df = df.loc[df['strand'] + df['reverse'] != 0]
            align_rate = (df['strand'] / (df['strand'] + df['reverse'])).median()
            if align_rate >= 0.8:
                strand_inf = "fr-secondstrand"
            elif align_rate <= 0.2:
                strand_inf = "fr-firststrand"
            else:
                strand_inf = "fr-unstranded"
            return strand_inf
        except Exception:
            continue
    return "fr-unstranded"


def parse_cmd(ctrls, bams, gtf, libtype, output, n_jobs, paired, length, **kwargs):
    summ = os.path.join(output, "summary.txt")
    if os.path.exists(summ) and datetime.fromtimestamp(os.path.getctime(summ)) >= datetime(2024, 5, 20):
        try:
            df = pd.read_csv(summ, sep="\t")
            if df.shape == (5, 10):
                logger.info(f"skip: {output}")
                return ""
        except Exception:
            pass

    if ctrls and bams:
        os.makedirs(os.path.join(output, "tmp"), exist_ok=True)

        b1 = os.path.join(output, "b1.txt")
        with open(b1, "w+") as w:
            w.write(",".join(sorted(bams)))

        b2 = os.path.join(output, "b2.txt")
        with open(b2, "w+") as w:
            w.write(",".join(sorted(ctrls)))

        paired = 'paired' if paired else 'single'

        gtf = REFERENCES.get(gtf)

        cmd = f"python /mnt/raid62/zym/env/rmats-turbo/rmats.py --b1 {b1} --b2 {b2} --gtf {gtf} --variable-read-length -t {paired} --libType {libtype} --readLength {length} --nthread {n_jobs} --od {output} --tmp {output}/tmp"
        return cmd


def guess_seq_length(files):
    seqLen = []
    for file in files:
        try:
            with pysam.AlignmentFile(file) as r:
                for idx, rec in enumerate(r):
                    seqLen.append(len(rec.query_sequence))

                    if idx > 10:
                        break
        except (ValueError, OSError) as err:
            logger.error(f"{file} {err}")
            return 0
    try:
        return int(np.mean(seqLen))
    except ValueError as err:
        return 50


def genounzip(files: str):
    for file in glob(files):
        check_call(f"genounzip -f -^ {file}", shell = True)
        if os.path.exists(file):
            os.remove(file)


def update_records(proj: str, perturb_id: str, cmd: str, error, path, server):
    if RMats.get_or_none(RMats.perturb_id == perturb_id):
        RMats.update(**{
            "project_accession": proj,
            "cmd": cmd, "path": path, "properly": error,
            "server": str(server),
        }).where(RMats.perturb_id == perturb_id).execute()
    else:
        RMats.create(**{
            "project_accession": proj, "perturb_id": perturb_id,
            "cmd": cmd, "path": path, "properly": error,
            "server": str(server),
        })


def run(data):
    cmd = parse_cmd(**data)

    if cmd:
        error = ""
        try:
            w = open(os.devnull)
            check_call(cmd, shell=True, stdout=w, stderr=w)
        except Exception as e:
            error = str(e)
        # update_records(data["proj"], data["perturb_id"], cmd, error, data["output"], data["output_server"])


def parse_pairs_from_meta_single(args):
    row = args["row"]
    star_datas = args["star_datas"]
    n_jobs = args["n_jobs"]
    output = args["output"]
    output_server = args["output_server"]
    CMDS = args["CMDS"]

    servers = set()
    bams, ctrls = [], []
    strandness, paired = None, None
    for query in row["perturb_label"].split(","):
        if query not in star_datas:
            bams = []
            break
        query = star_datas[query]
        if not paired and "_2.fastq.gz" in query['cmd']:
            paired = True
        if strandness is None:
            strandness = guess_strandness(glob(os.path.join(query['path'], "*.ReadsPerGene.out.tab*")))

        servers.add(query['server'])
        bams += glob(os.path.join(query['path'], query["experiment_accession"] + "*.bam"))

    ctrls = []
    for query in row["control_label"].split(","):
        if query not in star_datas:
            ctrls = []
            break
        query = star_datas[query]
        if not paired and "_2.fastq.gz" in query['cmd']:
            paired = True
        if strandness is None:
            strandness = guess_strandness(glob(os.path.join(query['path'], "*.ReadsPerGene.out.tab*")))

        servers.add(query['server'])
        ctrls += glob(os.path.join(query['path'], query["experiment_accession"] + "*.bam"))

    if ctrls and bams and len(servers) <= 2:
        org = row['org']
        if org == 'AraAth':
            org = "AraTha"

        seq_length = guess_seq_length(bams)

        if seq_length < 1:
            return

        lock.acquire()
        cmds = {}
        __CMDS__ = CMDS + f"_{output_server}.json"
        if os.path.exists(__CMDS__) and os.path.getsize(__CMDS__) > 10:
            with open(__CMDS__) as r:
                cmds = json.load(r)

        cmds[row["Perturb_ID"]] = {
            "ctrls": ctrls, "bams": bams, "gtf": org, "libtype": strandness,
            "output": os.path.join(output, row["Perturb_ID"]),
            "n_jobs": n_jobs, "paired": paired, "length": seq_length,
            "proj": row["study_group"], "perturb_id": row["Perturb_ID"], "output_server": output_server
        }

        with open(__CMDS__, "w+") as w:
            json.dump(cmds, w, indent=4)
        lock.release()


def init(l):
    global lock
    lock = l


def parse_pairs_from_meta(path: str, server: str, output: str, n_jobs: int, output_server: str):
    df = pd.read_csv(path, sep="\t")

    # star_datas = {}
    # for query in StarData.select().where(StarData.server.in_([str(server), "1130"])).dicts():
    #     star_datas[query['experiment_accession']] = query
    #
    # refs = {}
    # for query in Ref.select():
    #     fs = [x for x in glob(os.path.join(query.path, "*.gtf")) if "chr" not in x and "sorted" not in x and "filter" not in x]
    #     refs[query.sci_name] = fs[0]
    #
    # cmds = {}
    # CMDS = __CMDS__ + f"_{output_server}.json"
    # if os.path.exists(CMDS) and os.path.getsize(CMDS) > 10:
    #     with open(CMDS) as r:
    #         cmds = json.load(r)
    #
    # args = []
    # for _, row in tqdm(df.iterrows(), total=df.shape[0], desc="Load meta"):
    #
    #     # if RMats.get_or_none(RMats.perturb_id == row["Perturb_ID"]) or row["Perturb_ID"] in cmds.keys():
    #     #     continue
    #
    #     args.append({
    #         "row": row, "star_datas": star_datas,
    #         "n_jobs": n_jobs, "output": output,
    #         "output_server": output_server,
    #         "CMDS": __CMDS__
    #     })
    #
    # l = Lock()
    # with Pool(n_jobs, initializer=init, initargs=(l, )) as p:
    #     list(tqdm(p.imap(parse_pairs_from_meta_single, args), total=len(args), desc="Generating cmds..."))
    CMDS = __CMDS__ + f"_{output_server}.json"
    if os.path.exists(CMDS) and os.path.getsize(CMDS) > 10:
        with open(CMDS) as r:
            cmds = json.load(r)
    return cmds


def main(meta, output, output_server, server: str = "1130", species: str = None, n_jobs=30, parallel = 10):
    cmds = parse_pairs_from_meta(meta, server, output, n_jobs, output_server)
    if not RMats.table_exists():
        RMats.create_table()

    global REFERENCES
    if not REFERENCES:
        REFERENCES = get_references()

    print("run")
    with Pool(parallel) as p:
        list(tqdm(p.imap_unordered(run, cmds.values()), total=len(cmds), desc="Running..."))


if __name__ == "__main__":
    from fire import Fire
    Fire(main)
