#!/usr/bin/env python
# -*- coding:utf-8 -*-
u"""
Database handle the meta info
"""
import peewee as pw

from playhouse.postgres_ext import *


ext_db = pw.Proxy()


def init_db(name: str, user: str, password: str, host: str, port: int):
    if user and password:
        print("init postgres")
        db = PostgresqlExtDatabase(
            name, user=user, password=password,
            host=host.strip('"'), port=port
        )
    else:
        print("init sqlite")
        db = SqliteDatabase(name)
    ext_db.initialize(db)


class Model(Model):
    class Meta:
        database = ext_db


class Ref(Model):
    scientific_name = pw.CharField(index=True, unique=True)
    sci_name = pw.CharField(index=True, unique=True)
    path = pw.TextField()
    star = pw.CharField()
    salmon = pw.CharField()

    class Meta:
        table_name = "reference"


class Meta(Model):
    u"""
    Meta info
    """
    study_accession = pw.TextField(index=True)
    secondary_study_accession = pw.TextField(index=True)
    sample_accession = pw.TextField(index=True, null=True)
    secondary_sample_accession = pw.TextField(index=True, null=True)
    experiment_accession = pw.TextField(index=True, null=True)
    run_accession = pw.TextField(index=True)
    submission_accession = pw.TextField(null=True)
    tax_id = pw.IntegerField(index=True)
    scientific_name = pw.TextField(index=True)
    sci_name = pw.TextField(index=True)
    instrument_platform = pw.TextField(null=True)
    instrument_model = pw.TextField(null=True)
    library_name = pw.TextField(null=True)
    nominal_length = pw.FloatField(null=True)
    library_layout = pw.TextField(null=True)
    library_strategy = pw.TextField(null=True)
    library_source = pw.TextField(null=True)
    library_selection = pw.TextField(null=True)
    read_count = pw.BigIntegerField(null=True)
    base_count = pw.BigIntegerField(null=True)
    center_name = pw.TextField(null=True)
    first_public = pw.TextField(null=True)
    last_updated = pw.TextField(null=True)
    experiment_title = pw.TextField(null=True)
    study_title = pw.TextField(null=True)
    study_alias = pw.TextField(null=True)
    experiment_alias = pw.TextField(null=True)
    run_alias = pw.TextField(null=True)
    fastq_bytes = pw.CharField(null=True)
    fastq_md5 = pw.TextField(null=True)
    fastq_ftp = pw.TextField(null=True)
    fastq_aspera = pw.TextField(null=True)
    fastq_galaxy = pw.TextField(null=True)
    submitted_bytes = pw.TextField(null=True)
    submitted_md5 = pw.TextField(null=True)
    submitted_ftp = pw.TextField(null=True)
    submitted_aspera = pw.TextField(null=True)
    submitted_galaxy = pw.TextField(null=True)
    submitted_format = pw.TextField(null=True)
    sra_bytes = pw.BigIntegerField(null=True)
    sra_md5 = pw.TextField(null=True)
    sra_ftp = pw.TextField(null=True)
    sra_aspera = pw.TextField(null=True)
    sra_galaxy = pw.TextField(null=True)
    sample_alias = pw.TextField(null=True)
    broker_name = pw.TextField(null=True)
    sample_title = pw.TextField(null=True)
    nominal_sdev = pw.TextField(null=True)
    first_created = pw.TextField(null=True)
    condition = pw.TextField()
    new_sci_name = pw.TextField()

    class Meta:
        table_name = "meta"


class RawData(Model):
    run_accession = pw.TextField(index=True)
    filename = pw.CharField(index=True)
    path = pw.TextField(null=True, index=True)
    md5 = pw.TextField(null=True)
    server = pw.CharField(null=True)
    properly = pw.BooleanField(index=True, default=False)
    download = pw.BooleanField(index=True, default=False)
    mismatch_md5 = pw.TextField(null=True)

    class Meta:
        table_name = "raw_data"


class StarData(Model):
    experiment_accession = pw.TextField(index=True, unique=True)
    path = pw.TextField(null=True)
    server = pw.CharField(index=True, null=True)
    cmd = pw.TextField(null=True)
    bam_md5 = pw.TextField(null=True)
    properly = pw.BooleanField(index=True, default=False)

    class Meta:
        table_name = "star_data"


class SalmonData(Model):
    experiment_accession = pw.TextField(index=True, unique=True)
    path = pw.TextField(null=True)
    server = pw.CharField(index=True, null=True)
    cmd = pw.TextField()
    properly = pw.BooleanField(index=True, default=False)

    class Meta:
        table_name = "salmon_data"


class RMats(Model):
    project_accession = pw.TextField(index=True)
    perturb_id = pw.TextField(index=True)
    path = pw.TextField()
    server = pw.CharField(index=True, null=True)
    cmd = pw.TextField()
    properly = pw.TextField(index=True)
    splice_tools = pw.BooleanField(index=True, default=False, null=False)

    class Meta:
        table_name = "rmats"


if __name__ == '__main__':
    pass
