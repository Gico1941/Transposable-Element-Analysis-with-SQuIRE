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
### Step 1. run squire_batch_run.sh, squire_call_batch.sh and squire_count_batch.sh to generate TE/gene expression matrix 

### Step 2. Alignment Quality control 
#### Visualization of STAR mapping summary with Mapped_summary.R
<img src="https://github.com/Gico1941/Transposable-Element-Analysis-with-SQuIRE/assets/127346166/dfd80791-0f77-4271-bfcd-464700e7b5fb" width="500" />

### Step 3. Differential expression analysis with Differential_analysis_plot.R
#### PCA
<img src="https://github.com/Gico1941/Transposable-Element-Analysis-with-SQuIRE/assets/127346166/b3ce81c2-aa68-4e85-a7f5-70aee61a788d" width="500" />

#### Volcano plot of gene 
<img src="https://github.com/Gico1941/Transposable-Element-Analysis-with-SQuIRE/assets/127346166/9aa56c33-750c-446e-b146-f2e2dc7cbebc" width="500" />

#### GSEA (https://www.gsea-msigdb.org/gsea/index.jsp)

#### Step 4. Hommer annotation to study genomic distribution of Genes and TEs to study their regulatory association
##### Step 4.1 Hommer installatiton :
  ```
1. Download homer config : wget file http://homer.ucsd.edu/homer/configureHomer.pl                save it to /root
2. install homer : perl configureHomer.pl -install homer  (may counter require for basic softwares)
3. add homer to system path : PATH=$PATH:/root/bin/   (PATH=$PATH:/root/homer/bin)   (or PATH=$PATH:   + the directory of /bin of homer)

environment preparation:

4. download ghostscript : wget https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs1000/ghostscript-10.0.0.tar.gz  (origin site : http://pages.cs.wisc.edu/~ghost/)
5. tar zxvf ghostcript-10.0.0.tar.gz
6. cd ghostcript-10.0.0
7. run all codes below in a row:
		./configure
            make
            make all
            sudo make install
8. download weblogo : wget http://weblogo.berkeley.edu/release/weblogo.2.8.2.tar.gz  (origin site : http://weblogo.berkeley.edu/)
7. tar zxvf weblogo.2.8.2.tar.gz
9. add weblogo into PATH :  PATH=$PATH:/root/weblogo/  (or PATH=$PATH:   + the directory of /bin of weblogo)


10(optional for genome position based motif finding): wget https://genome-test.gi.ucsc.edu/~kent/exe/linux/blatSuite.37.zip
11. unzip blatSuite.37.zip
12. add into path : PATH=$PATH:/root/blat
  ```
##### Step 4.2 generate BED format files for hommer process with Homer_file_generation_for_DE_genes.R and Homer_file_generation_for_TEs_only.R 
(up_xxx   : up-regulated TE/gene (log2fc > 2, p < 10^-4))
(down_xxx   : down-regulated TE/gene (log2fc < -2, p < 10^-4))

##### Step 4.3











