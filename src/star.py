#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Find and update current fastq files
"""

import os
from glob import glob
from multiprocessing import Pool
from shutil import rmtree
from subprocess import check_call, CalledProcessError

import pandas as pd
import pysam
from tqdm import tqdm
from peewee import fn, SQL

from db.models import RawData, StarData, Meta, Ref, ext_db
from db.utils import calculate_md5


def move_file(infile: str, outfile: str):
    size = 1

    fsize = os.path.getsize(infile)
    while (fsize / (1024 * size)) > 10:
        size += 1

    pbar = tqdm(total=round(fsize / (1024 * size), ndigits=2), desc=os.path.basename(infile))
    with open(outfile, "wb+") as w:
        with open(infile, "rb") as r:
            while buf := r.read(4096):
                pbar.update(len(buf) / (1024 * size))
                w.write(buf)

    os.remove(infile)


def move_star(indir: str, outdir: str):
    for parent, _, files in os.walk(indir):
        for f in files:
            f = os.path.join(parent, f)
            o = f.replace(indir, outdir)

            if (os.path.exists(o) and os.path.getsize(o) < os.path.getsize(f)) or not os.path.exists(o):
                os.makedirs(os.path.dirname(o), exist_ok=True)

                try:
                    os.rename(f, o)
                except Exception:
                    move_file(f, o)
            elif os.path.exists(o) and os.path.getsize(o) == os.path.getsize(f):
                os.remove(f)
            elif os.path.exists(o) and os.path.getsize(o) > os.path.getsize(f):
                print(f, o, os.path.getsize(f), os.path.getsize(o))


def __extract_cmd_from_bam__(path: str):
    with pysam.AlignmentFile(path) as r:
        for k, v in r.header.items():
            if k == "CO":
                v = v[0]
                return v.replace("user command line: ", "").strip()
    return ""


def generate_path(indir: str, server: int = 1130):
    if not StarData.table_exists():
        StarData.create_table()
    for parent, _, files in tqdm(os.walk(indir)):
        for f in tqdm(files):
            if f.endswith("Aligned.sortedByCoord.out.bam"):
                failed = False
                if not os.path.exists(os.path.join(parent, f) + ".bai"):
                    tqdm.write(f)
                    try:
                        pysam.index(os.path.join(parent, f), "-@ 20")
                    except Exception:
                        failed = True
                key = os.path.basename(f).split(".")[0]

                if not StarData.get_or_none(experiment_accession=key) and Meta.get_or_none(experiment_accession=key):
                    md5 = calculate_md5(os.path.join(parent, f))
                    cmd = __extract_cmd_from_bam__(os.path.join(parent, f)) if not failed else ""

                    if StarData.get_or_none(StarData.experiment_accession == key):
                        StarData.update(**{
                            "experiment_accession": key,
                            "cmd": cmd, "bam_md5": md5,
                            "properly": not failed, "path": parent,
                            "server": str(server)
                        }).execute()
                    else:
                        StarData.create(**{
                            "experiment_accession": key,
                            "cmd": cmd, "bam_md5": md5,
                            "properly": not failed, "path": parent,
                            "server": str(server)
                        })
                # w.write(f"{key}\t{parent}\t{server}\n")

                # rec = Meta.get_or_none(experiment_accession=key)
                # if rec and rec.secondary_study_accession != os.path.basename(parent):
                #     os.makedirs(os.path.join(parent, rec.secondary_study_accession), exist_ok=True)
                #     for g in glob(os.path.join(parent, key + "*")):
                #         os.rename(g, os.path.join(parent, rec.secondary_study_accession, os.path.basename(g)))
                #     # print(key, rec.secondary_study_accession, parent)
                #     # exit(0)


def __get_exist_experiment_ids__(species: str = None):
    queries = Meta.select(
        StarData.experiment_accession.distinct(),
    ). \
        join(StarData, on=(Meta.experiment_accession == StarData.experiment_accession))\
    .where(StarData.properly == True)

    if species:
        queries = queries.where(Meta.sci_name.in_(species.split(",")))

    return set([x.experiment_accession for x in queries])


def dump_(path, data: str):
    with open(path, "a+") as w:
        w.write(data)


class StarRequiredData:

    def __init__(self, data: dict, writer = None):
        self.experiment = data["experiment_accession"]
        self.run = {data["run_accession"]: []}
        self.sci_name = data["sci_name"]
        self.path = data["path"]
        self.server = data["server"]
        self.output_server = data["output_server"]
        self.ref = None
        self.__get_fqs__()
        self.writer = writer

    def __get_fqs__(self):
        for k in self.run.keys():
            paths = []
            for i in RawData.select().where((RawData.run_accession == k) & (RawData.server.in_(self.server))):
                if i.properly:
                    paths.append(i.path)
                else:
                    paths = []
                    break

            self.run[k] = paths

    def __add__(self, other):
        if self.experiment == other.experiment:
            self.run.update(other.run)
        return self

    def __str__(self):
        u"""
        STAR --runThreadN 20
        --genomeDir /mnt/raid61/Ref/HomSap/release101/2.7.10b
        --outSAMtype BAM SortedByCoordinate --outBAMcompression 9
        --limitBAMsortRAM 100000000000 --readFilesCommand zcat
        --quantMode GeneCounts
        --readFilesIn /mnt/raid61/public_data/Perturbation_data/rawdata/Homo_sapiens/SRX8148326/SRR11580626.fastq.gz
        --outFileNamePrefix /mnt/data3/zym/STAR/HomSap/SRP257714/SRX8148326.
        --genomeLoad LoadAndKeep
        :return:
        """
        STAR = "STAR --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outBAMcompression 9 " \
               "--readFilesCommand zcat --quantMode GeneCounts --limitBAMsortRAM 119590295963 --genomeLoad LoadAndKeep"
        fq1, fq2, fq = [], [], []
        for fqs in self.run.values():
            if not fqs or any([x is None for x in fqs]):
                return ""
            for f in fqs:
                if "_1" in os.path.basename(f):
                    fq1.append(f)
                elif "_2" in os.path.basename(f):
                    fq2.append(f)
                else:
                    fq.append(f)

        if len(fq1) == len(fq2) > 0:
            return f"{STAR} --readFilesIn {','.join(fq1)} {','.join(fq2)} --genomeDir {self.ref} --outFileNamePrefix {self.path}/{self.experiment}."
        elif fq:
            return f"{STAR} --readFilesIn {','.join(fq)} --genomeDir {self.ref} --outFileNamePrefix {self.path}/{self.experiment}."
        return ""


def __run_star__(data: StarRequiredData):
    if data is not None and str(data):
        for f in glob(os.path.join(data.path, data.experiment + "._STARtmp")):
            if os.path.isdir(f):
                rmtree(f)

        finished = False
        logs = os.path.join(data.path, data.experiment + ".Log.progress.out")
        if os.path.exists(logs):
            with open(logs) as r:
                for line in r:
                    if "ALL DONE" in line:
                        finished = True

        logs = os.path.join(data.path, data.experiment + ".Log.progress.out.gz")
        if os.path.exists(logs):
            import gzip
            with gzip.open(logs, "rt") as r:
                for line in r:
                    if "ALL DONE" in line:
                        finished = True

        if not finished:
            if data.writer is not None:
                dump_(data.writer, str(data) + "\n")

        bam = glob(os.path.join(data.path, data.experiment + "*.bam"))

        rec = StarData.get_or_none(StarData.experiment_accession == data.experiment)
        if rec is not None:
            return "update", {
                "experiment_accession": data.experiment,
                "cmd": str(data), "bam_md5": None,  # calculate_md5(bam),
                "properly": len(bam) > 0, "path": data.path,
                "server": data.output_server
            }
        else:
            return "create", {
                "experiment_accession": data.experiment,
                "cmd": str(data), "bam_md5": None,  # calculate_md5(bam),
                "properly": len(bam) > 0, "path": data.path,
                "server": data.output_server
            }


def __selected_raw_data__(server, species):
    run_accessions = set()
    queries = RawData.select(). \
        where((RawData.properly == True) & (RawData.server == server))

    if species:
        queries = queries.where(Meta.sci_name.in_(species.split(",")))

    for query in queries:
        run_accessions.add(query.run_accession)

    return Meta.select().where(Meta.run_accession.in_(run_accessions))


def run_star(output: str = None, species: str = None, server: str = "1130", n_jobs=1, output_server: str = None, bash = None):
    server = str(server)
    if output_server is None:
        output_server = server

    # 先生成没跑完STAR的样本records
    # print("获取failed records")
    queries = Meta.select(
        Meta.sci_name, Meta.secondary_study_accession,
        Meta.run_accession, StarData.path,
        StarData.experiment_accession,
    ). \
        join(StarData, on=(Meta.experiment_accession == StarData.experiment_accession)). \
        where((StarData.properly != True) & (StarData.server == server))

    if species:
        queries = queries.where(Meta.sci_name.in_(species.split(",")))

    data = [x for x in queries.dicts()]

    if not output_server:
        output_server = str(server)

    # data = []
    if output:
        newdata = []
        for d in data:
            d["path"] = os.path.join(output, d["sci_name"], d["secondary_study_accession"])
            newdata.append(d)
        data = newdata
        del newdata

        # print("获取not properly records")
        # 如果指定了输出目录，把直接没有star的都给提出来
        os.makedirs(output, exist_ok=True)

        # 获取有记录的experiment_accession
        ids = __get_exist_experiment_ids__(species)
        queries = __selected_raw_data__(server, species)

        if server == "1130":
            server = [server]
        else:
            server = [server, "1130"]

        # 获取已经下好，但没跑过的记录
        required_ids = set()
        for i in queries:
            if i.experiment_accession not in ids:
                required_ids.add(i.experiment_accession)

        queries = Meta.select().where(Meta.experiment_accession.in_(required_ids))
        if species:
            queries = queries.where(Meta.sci_name.in_(species.split(",")))

        for query in queries.dicts():
            query["path"] = os.path.join(output, query["sci_name"], query["secondary_study_accession"])
            data.append(query)

    if bash:
        with open(bash, "w+") as w:
            w.write("#!/bin/bash\n")

    # print("配置fastq文件")
    formatted_data = {}  # species: experiment_accession: StarRequiredData
    for d in data:
        key = d["sci_name"]
        experiment = d["experiment_accession"]
        d["server"] = server
        d["output_server"] = output_server

        val = StarRequiredData(d, writer=bash)

        ref = Ref.get(Ref.sci_name == key)
        val.ref = os.path.join(ref.path, ref.star)
        if key not in formatted_data:
            formatted_data[key] = {experiment: val}
        elif experiment not in formatted_data[key]:
            formatted_data[key][experiment] = val
        else:
            formatted_data[key][experiment] = formatted_data[key][experiment] + val
            # print(formatted_data[key][experiment])

    # print("Run STAR")
    for specie, data in tqdm(formatted_data.items()):
        if len(data) > 0:
            # for k, v in data.items():
            #     print(k, str(v))

            data = list(data.values())
            cmd = f"STAR --genomeLoad LoadAndExit --genomeDir {data[0].ref}"

            if bash is not None:
                dump_(bash, cmd + "\n")
            else:
                check_call(cmd, shell=True)
            with Pool(n_jobs) as p:
                configs = list(tqdm(p.imap_unordered(__run_star__, data), total=len(data), desc=specie))

            cmd = f"STAR --genomeLoad Remove --genomeDir {data[0].ref}"

            if bash is not None:
                dump_(bash, cmd + "\n")
            else:
                check_call(cmd, shell=True)

            with ext_db.atomic():
                for conf in tqdm(configs):
                    if conf is not None:
                        process, conf = conf
                        if process == "create":
                            StarData.create(**conf)
                        if process == "update":
                            StarData.update(**conf).where(
                                StarData.experiment_accession == conf["experiment_accession"]).execute()


def dump(output: str):
    data = []
    for query in StarData.select(
            StarData.experiment_accession, StarData.path, StarData.server,
            Meta.sci_name, Meta.tax_id
    ). \
            join(Meta, on=(Meta.experiment_accession == StarData.experiment_accession)). \
            where(StarData.properly == True).dicts():
        data.append(query)
        if len(data) > 10:
            break

    data = pd.DataFrame(data)
    data.to_excel(data, output)


def __check__(experiment):
    ok, failed = None, []
    for rec in StarData.select().where(StarData.experiment_accession == experiment):
        if rec.path is not None and os.path.exists(rec.path) and glob(
                os.path.join(rec.path, rec.experiment_accession + ".*.bam*")) and rec.properly and ok is None:
            ok = rec.id
        else:
            failed.append(rec.id)

    return ok, failed


def check():
    experiments = StarData.select(StarData.experiment_accession, fn.COUNT(StarData.id).alias("count")). \
        group_by(StarData.experiment_accession).order_by(SQL("count"), SQL("count").desc())

    for e in experiments:
        if e.count > 1:
            print(e.experiment_accession, e.count)
            ok, failed = __check__(e.experiment_accession)
            StarData.delete().where(StarData.id.in_(failed)).execute()


def discover_failed_rawdata(server: str = "1130"):
    query = StarData.select(
        StarData.experiment_accession, StarData.server, RawData.path, Meta.fastq_aspera
    ).join(
        Meta, on=(StarData.experiment_accession == Meta.experiment_accession)
    ).join(
        RawData, on=(Meta.run_accession == RawData.run_accession)
    ).where((StarData.properly == False) & (RawData.server == str(server))).dicts()

    for q in query:
        print(q)


if __name__ == '__main__':
    from fire import Fire

    Fire({
        "move": move_star,
        "find": generate_path,
        "run": run_star,
        "dump": dump,
        "check": check,
        "discover": discover_failed_rawdata
    })
