# Molecular Docking with Vina: 

**This code is very unsafe**
## Installation:
1. Install all the required libaries from ```environment.yml```. Hopefully you have anaconda installed.

2. Build openbabel Dockerfile from OpenBabel. 
```
$ cd OpenBabel 

$ docker build -t myobabel . 

```

## Running: 
```
$uvicorn main:app --reload
```

## Usage: 
1. Upload pdb file of both protein and ligand at /uploadPDBFiles

For binding pocket:

http://cao.labshare.cn/drugrep/

https://www.playmolecule.com/deepsite/


*Work In Progress*