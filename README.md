# HyperVR: a hybrid deep ensemble learning approach for simultaneously predicting virulence factors and antibiotic resistance genes
To better simultaneous predict VFs and ARGs, we propose a hybrid deep ensemble learning approach called HyperVR. By considering both best hit scores and statistical gene sequence patterns, HyperVR combines classical machine learning and deep learning to simultaneously and accurately predict VFs, ARGs and negative genes (neither VFs nor ARGs). This repository contains codes and datas for HyperVR model.
# Data description
Uniprot_ARG.fasta: 2000 antibiotic resistance genes in the UNIPROT database were used for model training and validation.

Uniprot_ARG_ind.fasta: 209 independent antibiotic resistance genes in the UNIPROT database.

Uniprot_VF.fasta: 2000 virulence factors in the UNIPROT database were used for model training and validation.

Uniprot_ARG_ind.fasta: 209 independent virulence factors in the UNIPROT database.

Uniprot_NS.fasta: 2000 negative genes (neither VFs nor ARGs) in the UNIPROT database were used for model training and validation.

Uniprot_NS_ind.fasta: 209 independent negative genes (neither VFs nor ARGs) in the UNIPROT database.

# Requirements
HyperVR is tested to work under:

Python 3.8

Tensorflow 2.8.0

Keras 2.8.0

numpy 1.21.2

sklearn 1.1.1

# License
This source code is licensed under the MIT license found in the LICENSE file in the root directory of this source tree.

# Contacts
If you have any questions or comments, please feel free to email: byj@hnu.edu.cn.
