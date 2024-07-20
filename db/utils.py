#!/usr/bin/env python3
# -*- coding: utf-8 -*-

u"""
Database related utils
"""
import hashlib
import os
from typing import Dict, List, Optional

from peewee import DataError, Model, FieldAccessor
from tqdm import trange

from db.models import Meta, RawData, Ref, ext_db


def model_attrs(table: Model) -> Dict[str, FieldAccessor]:
    res = {}
    for key, value in table.__dict__.items():
        if isinstance(value, FieldAccessor):
            res[key] = value
    return res


def bulk_insert(data: List, table: Model, trunc: bool = False):
    if table.table_exists() and trunc:
        table.drop_table()

    if not table.table_exists():
        table.create_table()

    with ext_db.atomic():
        for idx in trange(0, len(data), 100, desc="Inserting..."):
            try:
                table.insert_many(data[idx:idx + 100]).execute()
            except DataError as err:
                print(err)
                print(data[idx])
                exit(0)


def calculate_md5(path: str):
    m = hashlib.md5()

    # 讀取檔案內容，計算 MD5 雜湊值
    with open(path, "rb") as f:
        buf = f.read()
        m.update(buf)

    return m.hexdigest()


def reference(source: str = "star"):
    ref = {}
    for rec in Ref.select().dicts():
        ref[rec['sci_name'].strip()] = os.path.join(rec['path'].strip(), rec[source].strip())

    return ref


def raw_data(
        server: Optional[str] = None,
        experiment_accession: Optional[List[str]] = None,
        run_accession: Optional[List[str]] = None
) -> Dict[str, Dict[str, List[str]]]:
    query = RawData.select(
        Meta.sci_name,
        Meta.experiment_accession,
        RawData.run_accession,
        RawData.path
    )\
        .join(Meta, on=(Meta.run_accession == RawData.run_accession))\
        .where(RawData.properly==True)

    if experiment_accession:
        query = query.where(Meta.experiment_accession.in_(experiment_accession))

    if run_accession:
        query = query.where(RawData.run_accession.in_(run_accession))

    if server:
        if server == "1130":
            server = [server]
        else:
            server = [server, "1130"]
        query = query.where(RawData.server.in_(server))

    fqs = {}
    for row in query.dicts():
        check = RawData.select().\
            join(Meta, on=(Meta.run_accession == RawData.run_accession)).\
            where(Meta.experiment_accession == row['experiment_accession'].strip())

        if all([x.properly for x in check]):
            if row["sci_name"].strip() not in fqs.keys():
                fqs[row["sci_name"].strip()] = {}
            if row['experiment_accession'].strip() not in fqs[row["sci_name"].strip()].keys():
                fqs[row["sci_name"].strip()][row['experiment_accession'].strip()] = []
            fqs[row["sci_name"].strip()][row['experiment_accession'].strip()].append(row['path'].strip())

    return fqs




if __name__ == '__main__':
    pass
