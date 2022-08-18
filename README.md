# HyperVR: a hybrid deep ensemble learning approach for simultaneously predicting virulence factors and antibiotic resistance genes
To better simultaneous predict VFs and ARGs, we propose a hybrid deep ensemble learning approach called HyperVR. By considering both best hit scores and statistical gene sequence patterns, HyperVR combines classical machine learning and deep learning to simultaneously and accurately predict VFs, ARGs and negative genes (neither VFs nor ARGs). This repository contains codes and datas for HyperVR model.
# Data description

| File name  | Description |
| ------------- | ------------- |
| Uniprot_ARG.fasta  | Antibiotic resistance genes in the UNIPROT database were used for model training and validation  |
| Uniprot_ARG_ind.fasta  | Independent antibiotic resistance genes in the UNIPROT database  |
| Uniprot_VF.fasta  | Virulence factors in the UNIPROT database were used for model training and validation  |
| Uniprot_VF_ind.fasta| Independent virulence factors in the UNIPROT database |
| Uniprot_NS.fasta| Negative genes (neither VFs nor ARGs) in the UNIPROT database were used for model training and validation| 
| Uniprot_NS_ind.fasta|  Independent negative genes (neither VFs nor ARGs) in the UNIPROT database| 
| Uniprot_ARG+VF+NS.fasta|  Total of 3 types of genes in the UNIPROT database were used for model training and validation| 
| Label_ARG+VF+NS.csv|  Label of 3 types of genes in the UNIPROT database were used for model training and validation| 
| Database_GENE.zip| The remaining database genes were used as known ARGs and VFs| 

# Requirements
HyperVR is tested to work under:

Python 3.8

Tensorflow 2.8.0

Keras 2.8.0

numpy 1.21.2

sklearn 1.1.1

# Quick start
To reproduce our results:

1, Download uniref dataset and install the required toolkit
```
cd tools/ncbi-blast && wget -c https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz && tar -zxvf ncbi-blast-2.13.0+-x64-linux.tar.gz

cd tools/diamond && wget -c https://github.com/bbuchfink/diamond/releases/download/v2.0.5/diamond-linux64.tar.gz && tar -zxvf diamond-linux64.tar.gz 

cd tools/uniref50 && wget -c https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz && tar -zxvf uniref50.fasta.gz

/tools/ncbi-blast/ncbi-blast-2.13.0+/bin/makeblastdb -dbtype prot -in uniref50.fasta -input_type fasta -parse_seqids -out uniref50_blast

```
2, Run generate_pssm_profile.py to generate pssm profiles for each gene sequence, the options are:
```
python src/generate_pssm_profile.py --file /data/Uniprot_ARG+VF+NS.fasta --blastpgp /tools/ncbi-blast/ncbi-blast-2.13.0+/bin --db /tools/uniref50/uniref50_blast --outdir /src/pssm_profile

--file: protein sequence file in fasta format
--blastpgp: the path of NCBI psiblast program
--db: the path of unief50 database
--outdir: the path of out dir
```
3, Run generate_bitscore.py to generate bitscore features for each gene sequence, the options are:
```

python src/generate_bitscore.py --file /data/Uniprot_ARG+VF+NS.fasta --db_file /data/Database_GENE.fasta --diamond_path /tools/diamond/diamond --outdir /src/bitscore
--file: protein sequence file in fasta format
--db_file: protein sequence file in fasta format
--diamond_path: the path of diamond program
--outdir: the path of out dir
```


# License
This source code is licensed under the MIT license found in the LICENSE file in the root directory of this source tree.

# Contacts
If you have any questions or comments, please feel free to email: byj@hnu.edu.cn.
