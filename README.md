# msc_project
A semester long research project into the properties of Mesenchymal Stromal Cells (MSCs).

This project expanded on previous work completed by Dr. Raga Krishnakumar at Sandia National Laboratories.

The purpose of this project was to determine some of the genetic features that determine if a Mesenchymal Stromal Cell is antimicrobial or not.

My group was provided with mRNA data of ~9,000 MSCs, and experimentally determined labels of antimicrobial or non-antimicrobial. We sought to use statistical and computational methods better understand what made some MSCs antimicrobial, and what made some non-antimicrobial.


My contributions included:
- [A Leiden Clustering based Classification Method](leiden_experiment/Ben_seurat_experiment_leiden.md)
  - The purpose of using this cluster based classification method was to try to develop an
  understanding of the cellular features (i.e. gene expression
  signatures) that may underlie antimicrobial/non antimicrobial
  behavior.
  - This particular dataset is heterogeneous and high noise, which is why
  we hypothesize may be effective.
  - This methodology was suggested by Dr. Mansoor Haider.
- [Using Classic Supervised ML Methods on Dataset](Supervised_learning_exps.ipynb)
  - The purpose of this experiment was to recreate some previous work on the dataset.
  - I aimed to investigate the ability of classic ML models to accurately classify our data, using a limited selection of genes.
  - I intially used only a subset of genes as features that were provided to us by Dr. Raga Krishnakumar.
  - I also derived my own subset of genes to be used in the ML models for classification.
- [Gene Ontology Analysis, from Clustering](Cluster_pathwayID_test.ipynb)
  - Since the mRNA sequencing data gives an idea of what a cell is doing at the time it was sequenced, we aimed to better understand the cellular activites of different  cell clusters in the dataset.
  - First, I clustered our data using Leiden clustering with the Scanpy package.
  - Then, for each cluster, I found genes that were overexpressed in that cluster relative to the others.
  - Then, for each cluster, I used the GSEAPY package to query a database to identify pathways that were overexpressed (or enriched), based on the overexpressed genes.
  - I plotted these enriched pathways so we could get a quick sense of what cells in each cluster were "doing."

