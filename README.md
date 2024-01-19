# Transposable Element Analysis with SQuIRE

## Analysis overview

![Picture9](https://github.com/Gico1941/Transposable-Element-Analysis-with-SQuIRE/assets/127346166/ce406057-40e6-46c7-8675-18238d2c7f29)

## Description
A analysis pipeline that integrates SQuIRE, GSEA and Hommer for TE/gene expression quantification, annotation and regulatory association study.

## Software denpendencies
Ubuntu

R 

SQuIRE (https://github.com/wyang17/SQuIRE)

Hommer (http://homer.ucsd.edu/homer/)

## Results (use GSE211061 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211061 as example)
### Step 1. Data QC with multiQC (Overrepresented reads detected)
<img src="https://github.com/Gico1941/Transposable-Element-Analysis-with-SQuIRE/assets/127346166/7fdc5d41-e0bd-4057-ab7a-af0ee6c5029a" width="500" />

Subsequent Data QC after rRNA removal with sortMeRNA (Overrepresented reads detected)
<img src="https://github.com/Gico1941/Transposable-Element-Analysis-with-SQuIRE/assets/127346166/f9e4888c-c29a-4877-a4a4-037b6655843d" width="500" />

### Step 2. run squire_batch_run.sh, squire_call_batch.sh and squire_count_batch.sh to generate TE/gene expression matrix 

### Step 3. Alignment Quality control 
#### Visualization of STAR mapping summary with Mapped_summary.R
<img src="https://github.com/Gico1941/Transposable-Element-Analysis-with-SQuIRE/assets/127346166/dfd80791-0f77-4271-bfcd-464700e7b5fb" width="500" />

#### Identities of reads 
<img src="https://github.com/Gico1941/Transposable-Element-Analysis-with-SQuIRE/assets/127346166/9623be09-6ddf-44ce-b4c4-74d1e2472f79" width="500" />



