# perturb_processing

The raw data management and processing scripts used for https://perturbatlas.kratoss.site/ database


## Code structure

- db: contains database structure file for raw data and processing task management
- src: contains separate scripts to execute software for quality control, alignment and quantification software etc.
 - [trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)
 - [star](https://github.com/alexdobin/STAR)
 - [salmon](https://github.com/COMBINE-lab/salmon) 
 - [rmats](https://github.com/Xinglab/rmats-turbo)
 - [splicetools](https://github.com/flemingtonlab/SpliceTools)
 - etc.
- R: contains separate scripts to execute differential analysis, enrichment etc.
- config.ini: contains the basic configuration for database
- adapters: contains PE and SE adapters retrieved from `https://github.com/usadellab/Trimmomatic/tree/main/adapters`
  - PE.fa: merged by `NexteraPE-PE.fa`, `TruSeq2-PE.fa`, `TruSeq3-PE-2.fa`, `TruSeq3-PE.fa`
  - SE.fa: merged by `TruSeq2-SE.fa`, `TruSeq3-SE.fa`

    
## Prepare environments

1. basic python environments

```bash
poetry install

pooetry run python main.py --help
```

2. others tools

please follow the official tutorials to install required software and make sure executables properly pleased in your $PATH.

**Note**: for trimmomatic, please make sure the java is properly installed, and the jar file of trimmomatic placed under adapters directory.


## Usage

```help
python main.py --help

NAME
    main.py

SYNOPSIS
    main.py COMMAND

COMMANDS
    COMMAND is one of the following:

     splice
     rmats
     download
     star
     salmon
     sleuth
     trimmomatic
```

The usage for different subcommands:

```bash
python main.py download --help

NAME
    main.py download

SYNOPSIS
    main.py download <flags>

FLAGS
    --server=SERVER
        Type: str
        Default: '1130'
    -o, --outdir=OUTDIR
        Type: Optional[Optional]
        Default: None
    -a, --aria2c=ARIA2C
        Type: bool
        Default: False
    -g, --give_up=GIVE_UP
        Type: bool
        Default: False
    --species=SPECIES
        Type: Optional[Optional]
        Default: None
```


