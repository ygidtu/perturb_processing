# perturb_processing

The raw data management and processing scripts used for https://perturbatlas.kratoss.site/ database


## Code structure

- db: contains database structure file for raw data and processing task management
- src: contains separate scripts to execute software for quality control, alignment and quantification software etc.
- R: contains separate scripts to execute differential analysis, enrichment etc.
- config.ini: contains the basic configuration for database
- adapters: contains PE and SE adapters retrieved from `https://github.com/usadellab/Trimmomatic/tree/main/adapters`
  - PE.fa: merged by `NexteraPE-PE.fa`, `TruSeq2-PE.fa`, `TruSeq3-PE-2.fa`, `TruSeq3-PE.fa`
  - SE.fa: merged by `TruSeq2-SE.fa`, `TruSeq3-SE.fa`

    
## Prepare environments

### 1. basic python environments

```bash
poetry install

pooetry run python main.py --help
```

### 2. R packages required
   - clusterProfiler
   - data.table
   - DESeq2
   - dplyr
   - DRIMSeq
   - glue
   - GO.db
   - pbapply
   - stageR
   - stringr
   - reshape2
   - tidyr
   - rjson
   - tximport

### 3. the annotation database for multi-species required by enrichment analysis

  The content of `03.generate_pair/annotation.csv` required by `R/enrichment_split_run.R`. 
  - First column is the short scientific name for species.
  - Second column is the R packages or absolute path to generated `Gene_go_inf.Rds` by `R/prepare_go_inf.R` script.
  - Third column is the key type used by `clusterProfiler`.

  ```csv
  AraAth,org.At.tair.db,TAIR
  CaeEle,org.Ce.eg.db,ENSEMBL
  DanRer,org.Dr.eg.db,ENSEMBL
  DroMel,org.Dm.eg.db,ENSEMBL
  EscCol,org.EcK12.eg.db,SYMBOL
  GalGal,org.Gg.eg.db,ENSEMBL
  HomSap,org.Hs.eg.db,ENSEMBL
  MusMus,org.Mm.eg.db,ENSEMBL
  OrySat,OrySat/release52/Gene_go_inf.Rds,IG
  RatNor,org.Rn.eg.db,ENSEMBL
  SacCer,org.Sc.sgd.db,ENSEMBL
  SolLyc,SolLyc/release52/Gene_go_inf.Rds,IG
  SusScr,org.Ss.eg.db,SYMBOL
  ```

4. others tools

please follow the official tutorials to install required software and make sure executables properly pleased in your $PATH.

 - [trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)
 - [star](https://github.com/alexdobin/STAR)
 - [salmon](https://github.com/COMBINE-lab/salmon) 
 - [rmats](https://github.com/Xinglab/rmats-turbo)
 - [splicetools](https://github.com/flemingtonlab/SpliceTools)
 - etc.

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


