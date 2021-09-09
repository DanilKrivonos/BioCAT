# BioCAT: Biosynthesis Cluster Analysis Tool
<img src="https://user-images.githubusercontent.com/53526550/132544644-86306499-133d-44e2-8e4c-e2603fb7d0f0.png" width="290" height="200">BioCAT is the tool which find cluster of of biosynthesis of nonribosomal peptides(NRP). The tool using structure in SMILES format and trying to match restored chemiclal structure on condadate NRPS cluster.
# **Dependence:**
- antiSMASH 6.0.0 
### python3 librarys:
- sklearn 0.24.2
- pandas 1.2.5
- bipython 1.79
- numpy 1.21.0
- rdkit 2021.03.4
## **Installation:**
```git install https://github.com/DanilKrivonos/BioCAT```

## **Run:**
Simple run:
```
python3 BioCAT.py -smiles "SMILES" -genome "GENOME.fna" -out 
```
