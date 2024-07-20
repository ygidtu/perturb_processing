#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Find and update current fastq files
"""

import os
import psutil
from multiprocessing import Pool, Lock
from subprocess import check_call
from shutil import move, rmtree
from typing import Optional

from tqdm import tqdm

from db.models import Meta, RawData
from db.utils import bulk_insert, calculate_md5


def call(args):
    f, o = args

    tqdm.write(f"{f} -> {o}")
    # pbar = tqdm(total=os.path.getsize(f))
    with open(o, "wb+") as w:
        with open(f, "rb") as r:
            while buf := r.read(4096):
                # pbar.update(4096)
                w.write(buf)

    try:
        os.remove(f)
    except Exception:
        pass


def move_fastqs(
        indir: str, output: str, n_jobs: int = 10,
        left_threshold: int = 1024 * 1024 * 100,
        server: str = "1130", clean: bool = False
):
    server = str(server)
    hdd = psutil.disk_usage(output)
    total_disk = hdd.free

    fs = {}
    for parent, _, files in tqdm(os.walk(indir), desc="Walking..."):
        for f in files:
            if f.endswith("fastq.gz") and not os.path.exists(os.path.join(parent, f + ".aria2")):
                run_accession = os.path.basename(f).replace(".fastq.gz", "").split("_")[0]
                if run_accession not in fs.keys():
                    fs[run_accession] = []
                fs[run_accession].append(os.path.join(parent, f))
    print(len(fs))
    cmds = []
    for run_accession, files in tqdm(fs.items()):
        query = Meta.get_or_none(run_accession=run_accession)
        if query:
            if query.sci_name.strip() in ["HomSap", "MusMus"] and server == "1140":
                continue
            elif server == "1130" and query.sci_name.strip() != "HomSap":
                continue
            elif server == "1150" and query.sci_name.strip() != "MusMus":
                continue
            else:
                if os.path.basename(output) != query.sci_name.strip():
                    out = os.path.join(
                        output, query.sci_name.strip(),
                        query.secondary_study_accession,
                        query.experiment_accession
                    )
                else:
                    out = os.path.join(
                        output,
                        query.secondary_study_accession,
                        query.experiment_accession
                    )

                for f in files:
                    o = os.path.join(out, os.path.basename(f))
                    if f == o:
                        continue

                    try:
                        os.makedirs(out, exist_ok=True)
                    except FileExistsError as err:
                        print(o)
                        print(err)
                        exit(0)
                    if not clean:
                        if os.path.exists(f) and (not os.path.exists(o) or os.path.getsize(f) > os.path.getsize(o)):
                            total_disk -= os.path.getsize(f)
                            if total_disk < left_threshold:
                                break

                            try:
                                cmds.append([f, o])
                                # move(f, o)
                            except Exception as err:
                                print(f, o)
                                print(err)
                                exit(0)
                        elif os.path.exists(o) and os.path.isdir(o) and os.path.dirname(f) == o:
                            total_disk -= os.path.getsize(f)
                            if total_disk < left_threshold:
                                break

                            move(f, o + ".temp")
                            rmtree(o)
                            move(o + ".temp", o)
                    elif clean and os.path.exists(o) and os.path.getsize(o) <= os.path.getsize(f):
                        os.remove(o)
                        # print(f, o)

    if cmds:
        with Pool(n_jobs) as p:
            list(tqdm(p.imap(call, cmds), total=len(cmds), desc="Moving..."))


def update_records(args):
    """
    run_accession: str, files: List[str]
    :param args:
    :return:
    """

    file, server, md5, give_up = args

    # query = RawData.get_or_none(filename=os.path.basename(file))
    # if query and (not query.properly or not query.path):

    properly = True if give_up else calculate_md5(file) == md5
    lock.acquire()
    RawData.update(properly=properly, server=server, path=file) \
        .where(RawData.filename == os.path.basename(file)) \
        .execute()
    lock.release()


def update(indir: str, n_jobs: int = 10, server: str = "1130", give_up = False):
    # check_call(
    #     f"{os.path.join(os.path.dirname(__file__), 'md5/perturb')} -i {indir} --server {server} -p {n_jobs}",
    #     shell=True
    # )

    query = RawData.select()
    filenames = {}
    for q in query:
        filenames[q.filename] = q

    fs = []
    for parent, _, files in tqdm(os.walk(indir), desc="Walking..."):
        for f in files:
            if f.endswith("fastq.gz") and f in filenames:
                info = filenames[f]
                # print(info.path, os.path.join(parent, f))
                if info.path == os.path.join(parent, f) and info.properly:
                    continue

                if "download" in parent:
                    continue
                fs.append([os.path.join(parent, f), server, info.md5, give_up])
    # print(fs)
    lock = Lock()

    def init(lock_):
        global lock
        lock = lock_

    with Pool(n_jobs, initializer=init, initargs=(lock,)) as p:
        list(tqdm(p.imap(update_records, fs), total=len(fs)))


# reduce duplicated records
def reduce():
    res = {}
    for query in RawData.select().dicts():
        query.pop("id")
        res[query["filename"]] = query
    bulk_insert(list(res.values()), RawData, True)


def __download__(args):
    __ASCP__ = "ascp"

    url, md5_, path, depth, aria2c, give_up, server = args.get("url"), args.get("md5"), args.get("out_path"), args.get(
        "depth", 0), args.get("aria2c"), args.get("give_up"), args.get("server")

    if depth > 5:
        return

    if aria2c:
        cmd = f"aria2c --file-allocation=none -s 20 -c -o {os.path.basename(path)} -d {os.path.dirname(path)} http://{url}"
    else:
        cmd = f"ascli -Pera server download --to-folder={os.path.basename(url)} {url.split(':')[-1]}"

    try:
        check_call(cmd, shell=True, cwd=os.path.dirname(path), stderr=open(os.devnull), stdout=open(os.devnull))
        print(cmd)
        if os.path.exists(path):
            md = calculate_md5(path)
            if give_up:
                RawData.update(**{
                    "properly": True, "download": True,
                    "mismatch_md5": md.strip(), "server": server,
                    "path": path
                }).where(RawData.filename == os.path.basename(path)).execute()
            else:
                md = calculate_md5(path)
                tqdm.write(f"{md} {md5_.strip()}")
                RawData.update(**{
                    "properly": md.strip() == md5_.strip(),
                    "download": True, "server": server,
                    "path": path
                }).where(RawData.filename == os.path.basename(path)).execute()
    except FileNotFoundError:
        return
    except Exception:
        args["depth"] = depth + 1
        __download__(args)


def download(
        server: str = "1130", outdir: Optional[str] = None,
        aria2c: bool = False, give_up: bool = False,
        species: Optional[str] = None):
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    queries = RawData.select().where(
        (RawData.server == str(server)) & (RawData.path != None) & (RawData.properly == False)
    ).order_by(RawData.filename)

    if species:
        # print(f"get records in species: {species}")
        species = species.split(",")
        queries = RawData.select() \
            .join(Meta, on=(RawData.run_accession == Meta.run_accession)) \
            .where(
            (RawData.path == None) & (RawData.properly == False) & (Meta.sci_name.in_(species))
        ).order_by(RawData.filename)

    cmds = []
    for query in queries:
        rec = Meta.get(run_accession=query.run_accession)

        for url, md5_ in zip(rec.fastq_aspera.split(";") if not aria2c else rec.fastq_ftp.split(";"),
                             rec.fastq_md5.split(";")):
            if os.path.basename(url) == query.filename:
                if not query.path:
                    out_path = os.path.join(rec.sci_name, rec.secondary_study_accession, rec.experiment_accession,
                                            query.filename)
                    if outdir:
                        out_path = os.path.join(outdir, out_path)

                    if not os.path.exists(os.path.dirname(out_path)):
                        os.makedirs(os.path.dirname(out_path), exist_ok=True)
                else:
                    out_path = query.path if not outdir else os.path.join(outdir, os.path.basename(query.path))

                if os.path.exists(out_path + ".aria2"):
                    cmds.append({
                        "url": url, "md5": md5_, "out_path": out_path,
                        "depth": 0, "aria2c": aria2c, "give_up": give_up, "server": str(server)
                    })

    for cmd in tqdm(cmds):
        __download__(cmd)



if __name__ == '__main__':
    from fire import Fire

    Fire({
        "move": move_fastqs,
        "update": update,
        "reduce": reduce,
        "download": download
    })
