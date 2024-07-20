#!/usr/bin/env python
from db.models import init_db
from src.spliceTools import main as spliceTools
from src.rmats import main as rmats
from src.download import download
from src.star import run_star
from src.salmon import main as salmon
from src.sleuth import main as sleuth
from src.trimmomatic import main as trimmomatic
from configparser import ConfigParser



if __name__ == '__main__':
    from fire import Fire

    cfg = ConfigParser()
    cfg.read("config.ini")

    init_db(**{x: y for x, y in cfg["db"].items()})

    Fire({
        "splice": spliceTools,
        "rmats": rmats,
        "download": download,
        "star": run_star,
        "salmon": salmon,
        "sleuth": sleuth,
        "trimmomatic": trimmomatic,
    })
