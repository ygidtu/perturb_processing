#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IU (an unstranded paired-end library where the reads face each other)
SF (a stranded single-end protocol where the reads come from the forward strand)
OSR (a stranded paired-end protocol where the reads face away from each other,
     read1 comes from reverse strand and read2 comes from the forward strand)
"""

# salmon quant -i transcripts_index -l A -1 <(zcat lib_1_1.fq.gz lib_2_1.fq.gz) -2 <(zcat lib_1_2.fq.gz lib_2_2.fq.gz) --validateMappings -o transcripts_quant -p
# salmon quant -i transcripts_index -l A -r <(zcat reads.fq.gz) --validateMappings -o transcripts_quant -p
from shutil import rmtree

import os
from multiprocessing import Pool
from subprocess import check_call, CalledProcessError
from typing import Dict, List

from tqdm import tqdm
from peewee import fn, SQL

from db.models import SalmonData, Meta, RawData, StarData
from db.utils import reference


def call(data: Dict[str, any]):
    if data["cmd"]:
        try:
            check_call(f"bash -c '{data['cmd']}'", shell=True)
            data["properly"] = True
        except CalledProcessError as err:
            print(err)
            data["properly"] = False

        try:
            if SalmonData.select().where(SalmonData.experiment_accession == data["experiment_accession"]).count() > 0:
                SalmonData.update(**data).where(SalmonData.experiment_accession == data["experiment_accession"]).execute()
            else:
                SalmonData.create(**data)
        except Exception as err:
            print(err)


def generate_cmd(output: str, files: List[str], sci_name: str, n_jobs: int = 10):
    REFS = reference("salmon")

    if sci_name.strip() not in REFS:
        return

    fqs1, fqs2, fqs = set(), set(), set()
    for f in files:
        if os.path.basename(f).endswith("_1.fastq.gz"):
            fqs1.add(f)
        elif os.path.basename(f).endswith("_2.fastq.gz"):
            fqs2.add(f)
        else:
            fqs.add(f)

    if (len(fqs1) < 1 and len(fqs2) < 1) and not fqs:
        return

    os.makedirs(output, exist_ok=True)

    cmd = f"salmon quant -i {REFS[sci_name.strip()]} -l A --validateMappings --numBootstraps 100 -o {output} -p {n_jobs}"
    if fqs1 and fqs2:
        cmd = f"{cmd} -1 <(zcat {' '.join(sorted(fqs1))}) -2 <(zcat {' '.join(sorted(fqs2))})"
    elif fqs:
        cmd = f"{cmd} -r <(zcat {' '.join(sorted(fqs))})"

    return f"{cmd} 2> {output}/salmon.log"


def __experiments__(server: str):
    for query in StarData.select(StarData.experiment_accession, Meta.sci_name)\
            .join(Meta, on=(StarData.experiment_accession == Meta.experiment_accession))\
            .where(StarData.server == server).dicts():
        yield query["sci_name"], query["experiment_accession"]


def __raw_data__(server: str, experiment: str):
    files = set()

    if server == "1130":
        server = [server]
    else:
        server = [server, "1130"]

    for query in RawData.select()\
            .join(Meta, on=(Meta.run_accession == RawData.run_accession))\
            .where(
                (Meta.experiment_accession == experiment) & (RawData.path != None) & (RawData.server.in_(server))
            ):
        files.add(query.path)
    return files


def main(output: str, server: str = "1130", n_jobs: int = 10, parallel: int = 3):
    server = str(server)
    if not SalmonData.table_exists():
        SalmonData.create_table()

    data = {}

    for sci_name, experiment in tqdm(__experiments__(server)):
        query = SalmonData.get_or_none(experiment_accession=experiment)
        if query and query.properly:
            continue

        files = list(__raw_data__(server, experiment))

        data[experiment] = {
            "experiment_accession": experiment,
            "path": os.path.join(output, sci_name, experiment),
            "server": server,
            "cmd": generate_cmd(
                os.path.join(output, sci_name.strip(), experiment),
                files, sci_name.strip(), n_jobs=n_jobs
            )
        }

    with Pool(parallel) as p:
        list(tqdm(p.imap(call, data.values()), total=len(data)))


def __check__(experiment):
    ok, failed = None, []
    for rec in SalmonData.select().where(SalmonData.experiment_accession == experiment):
        # print(rec.path, os.path.exists(rec.path), rec.properly)
        if rec.path is not None and rec.properly and ok is None:
            ok = rec.id
        else:
            failed.append(rec.id)

    return ok, failed


def check(server = "1130"):
    experiments = SalmonData.select(SalmonData.experiment_accession, fn.COUNT(SalmonData.id).alias("count")).\
            group_by(SalmonData.experiment_accession).order_by(SQL("count"), SQL("count").desc())

    for e in experiments:
        if e.count > 1:
            print(e.experiment_accession, e.count)
            ok, failed = __check__(e.experiment_accession)
            SalmonData.update(**{"deleted": True}).where(SalmonData.id.in_(failed)).execute()

    # for i in SalmonData.select().where(SalmonData.deleted == True):
    #     if os.path.exists(i.path):
    #         rmtree(i.path)
    # SalmonData.delete().where(SalmonData.deleted == True).execute()


def find(indir, server = "1140"):
    for parent, _, files in os.walk(indir):
        for f in files:
            if f == "salmon.log":
                with open(os.path.join(parent, f)) as r:
                    for line in r:
                        if "writing output" in line:
                            key = os.path.basename(parent)

                            if not SalmonData.get_or_none(SalmonData.experiment_accession == key):
                                SalmonData.create(**{
                                    "experiment_accession": key,
                                    "path": parent,
                                    "server": str(server),
                                    "cmd": "", "properly": True,
                                })


__ASCP__ = "/mnt/raid62/zym/env/aspera/cli"


def call1(cmd, depth=0):
    if depth > 5:
        return
    try:
        check_call(cmd, shell=True)
    except CalledProcessError:
        call1(cmd, depth=depth + 1)


def discover_failed_rawdata(server: str="1130"):
    query = SalmonData.select(
        SalmonData.experiment_accession, SalmonData.server, RawData.path, Meta.fastq_aspera
    ).join(
        Meta, on=(SalmonData.experiment_accession == Meta.experiment_accession)
    ).join(
        RawData, on=(Meta.run_accession == RawData.run_accession)
    ).where((SalmonData.properly == False) & (RawData.server == str(server))).dicts()

    for q in query:
        print(q)

        if os.path.exists(q["path"]):
            os.remove(q["path"])

        for u in q["fastq_aspera"].split(";"):
            if os.path.basename(q["path"]) == os.path.basename(u):
                cmd = f"{__ASCP__}/bin/ascp -QT -k1 -l 300m -P33001 --overwrite=diff+older -k 3 -i {__ASCP__}/etc/asperaweb_id_dsa.openssh era-fasp@{u} {q['path']}"
                call1(cmd)


if __name__ == '__main__':
    from fire import Fire

    Fire({
        "run": main, "check": check,
        "find": find, "discover": discover_failed_rawdata
    })
