# Integrative Transcriptomic and Machine Learning Approaches for Biomarker Discovery and Fairness Assessment in Lung Cancer

## Overview
This project focuses on identifying potential biomarkers in Lung Adenocarcinoma (LUAD) using transcriptomic analysis and evaluating fairness in machine learning-based classification models.

The study integrates:
- Differential gene expression analysis
- Functional enrichment analysis
- Protein-protein interaction (PPI) network analysis
- Hub gene identification
- Survival analysis
- ROC analysis
- Stage-wise expression analysis
- Fairness-aware machine learning evaluation

## Thesis Title
"Integrative Transcriptomic and Machine Learning Approaches for Biomarker Discovery and Fairness Assessment in Lung Cancer"

## Objectives
- Identify biomarkers involved in LUAD development and progression
- Evaluate diagnostic and prognostic significance of key hub genes
- Develop and assess fairness-aware machine learning models

## Key Findings
- Identified "2269 differentially expressed genes (DEGs)"
- Enrichment revealed pathways related to:
  - cell cycle regulation
  - extracellular matrix organization
  - cytokine-mediated signaling
- Hub genes identified:
  - "CCNB1"
  - "CENPA"
  - "TOP2A"
- Survival analysis showed these genes are associated with poor prognosis
- ROC analysis indicated strong diagnostic potential
- Fairness evaluation was performed across gender groups using machine learning

## Tools and Technologies
- R
- DESeq2
- biomaRt
- clusterProfiler
- STRING / Cytoscape
- Kaplan-Meier survival analysis
- ROC analysis
- SMOTE
- Machine Learning (Logistic Regression)
- Fairness metrics

## Repository Structure
scripts/   -> analysis scripts
data/      -> input/processed data
results/   -> figures and tables
