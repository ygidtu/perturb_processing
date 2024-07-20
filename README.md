# perturb_processing

The raw data management and processing scripts used for https://perturbatlas.kratoss.site/ database

## Prepare environments

```bash
poetry install

pooetry run python main.py --help
```


```help
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