#cd /mnt/h/LYL_RNAseq/GSE211061/TE_analysis
#conda activate squire
#mkdir output/mapped_with_human
#source ~/.profile
#bash squire_count_batch.sh
#squire Clean -r /mnt/h/LYL_RNAseq/GSE211061/TE_analysis/Ref/mm39_rmsk.txt -b mm39 -o /mnt/h/LYL_RNAseq/GSE211061/TE_analysis/Ref_cleaned



for i in Mapped/*;
	do

	echo processing" "$(basename $i )

	


	mkdir "Counts/"$(basename $i ) 

	/mnt/g/LYL_RNA_Workspace/TEanalysis/SQuIRE/SQuIRE/squire/Count.py -m "Mapped/"$(basename $i) -c ./Ref_cleaned -o Counts -f ./Ref -r 102 -p 36

	done	


#/mnt/g/LYL_RNA_Workspace/TEanalysis/SQuIRE/SQuIRE/squire/Count.py -m Mapped/SRR21017360 -c ./Ref_cleaned -o Counts -f ./Ref -r 102 -p 36
	

