#!/usr/bin/env python3
# -*- coding:utf-8 -*
u"""
generate database data
"""
import math
import os

import pandas as pd
import requests as rq

from tqdm import tqdm

from db.models import *
from db.utils import bulk_insert


def format_sci_name(name: str):
    name = name.split()
    if len(name) < 2:
        return name[0]
    return "".join([x[:3].title() for x in name])


def insert_meta(path: str, meta: str):
    df = pd.read_csv(path, sep="\t")
    df = df.drop_duplicates()
    sci_name = update_sci(meta)
    data = []
    for _, row in tqdm(df.iterrows(), total=df.shape[0]):
        temp = {}
        for key, val in row.to_dict().items():
            if key not in vars(Meta):
                continue

            try:
                if math.isnan(val):
                    continue
            except TypeError:
                pass

            if not val:
                continue

            if key == "scientific_name":
                temp[key] = sci_name.get(row["secondary_study_accession"], val)
                temp["sci_name"] = format_sci_name(temp[key])
                continue
            elif key in ["tax_id", "base_count", "read_count", "sra_bytes"]:
                val = int(val)
            elif key in ["nominal_length"]:
                val = float(val)
            temp[key] = val

        data.append(temp)

    bulk_insert(data, Meta, trunc=True)


def update_sci(path):
    df = pd.read_csv(path)
    res = {}
    for _, row in tqdm(df.iterrows(), total=df.shape[0]):
        res[row["study_accession"]] = row["scientific_name"]
        # Meta.update(
        #     scientific_name=row["scientific_name"],
        #     sci_name=format_sci_name(row["scientific_name"])
        # ).where(Meta.secondary_study_accession == row["study_accession"]).execute()
    return res


def update_meta(logs: str = "./logs"):
    os.makedirs(logs, exist_ok=True)

    fields = "fields=study_accession,secondary_study_accession," \
             "sample_accession,secondary_sample_accession," \
             "experiment_accession,run_accession,submission_accession," \
             "tax_id,scientific_name,instrument_platform,instrument_model," \
             "library_name,nominal_length,library_layout,library_strategy," \
             "library_source,library_selection,read_count,base_count," \
             "center_name,first_public,last_updated,experiment_title," \
             "study_title,study_alias,experiment_alias,run_alias,fastq_bytes," \
             "fastq_md5,fastq_ftp,fastq_aspera,fastq_galaxy,submitted_bytes," \
             "submitted_md5,submitted_ftp,submitted_aspera,submitted_galaxy," \
             "submitted_format,sra_bytes,sra_md5,sra_ftp,sra_aspera,sra_galaxy," \
             "sample_alias,broker_name,sample_title,nominal_sdev,first_created"

    query = Meta.select().where(Meta.fastq_md5.is_null()).dicts()

    for q in tqdm(query):
        try:
            l = os.path.join(logs, q['run_accession'])
            if os.path.exists(l):
                with open(l) as r:
                    data = json.load(r)
            else:
                r = rq.get(f"https://www.ebi.ac.uk/ena/portal/api/filereport?result=read_run&accession={q['experiment_accession']}"
                           f"&offset=0&limit=1000&format=json&{fields}")

                data = r.json()
                with open(l, "w+") as w:
                    json.dump(data, w)

            for row in data:
                if row["run_accession"] == q['run_accession']:
                    for k, v in q.items():
                        if v is None and row.get(k) and k != "id":
                            q[k] = row.get(k)

                    q.pop("id")
                    Meta.update(**q).where(Meta.run_accession == q['run_accession']).execute()
        except Exception as err:
            tqdm.write(f"{q['experiment_accession']}: {err}")


def generate_rawdata():
    data = []
    for query in Meta.select():
        if query.fastq_ftp and query.fastq_md5:
            for url, md5 in zip(query.fastq_ftp.split(";"), query.fastq_md5.split(";")):
                data.append({
                    "filename": os.path.basename(url),
                    "md5": md5, "run_accession": query.run_accession
                })

    bulk_insert(data, RawData, True)


def generate_ref(meta: str):
    df = pd.read_csv(meta, header=None)

    data = []
    for _, row in df.iterrows():
        data.append({
            'scientific_name': row[0],
            'sci_name': row[1],
            'path': os.path.dirname(row[2]),
            'star': os.path.basename(row[2]),
            'salmon': "1.10.2"
        })

    bulk_insert(data, Ref, True)


if __name__ == '__main__':
    # insert_meta(
    #     "/mnt/raid62/zhouran/project/2023-database/01.meta/data_sra_inf.used.tsv.gz",
    #     "/mnt/raid62/zym/project/dataexplorer/meta.csv"
    # )
    # update_meta()
    # generate_rawdata()
    # generate_ref("/mnt/raid61/Ref/ref.tsv")
    pass