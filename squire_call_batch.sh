#cd /mnt/h/LYL_RNAseq/GSE211061/TE_analysis

#conda activate squire
#mkdir output/mapped_with_human
#source ~/.profile
#squire Clean -r /mnt/h/LYL_RNAseq/GSE211061/TE_analysis/Ref/mm39_rmsk.txt -b mm39 -o /mnt/h/LYL_RNAseq/GSE211061/TE_analysis/Ref_cleaned
#bash squire_call_batch.sh


/mnt/g/LYL_RNA_Workspace/TEanalysis/SQuIRE/SQuIRE/squire/Call.py -i Counted -1 $(cat decoding/ctrl.txt) -2 $(cat decoding/lps.txt) -A CTRL -B LPS -o Call -p 36 -N GSE211061

	
#I prefer using awk. If there is only one column, use $0, else replace it with the last column.