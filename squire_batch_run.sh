#cd /mnt/h/LYL_RNAseq/GSE211061/TE_analysis
#conda activate squire
#mkdir output/mapped_with_human
#source ~/.profile
#bash squire_batch_run.sh

for i in /mnt/h/LYL_RNAseq/GSE211061/SRA2FASTQ/*;
	do
	if [ ! -f "Mapped/"$(basename $i)"/"$(basename $i )"_1SJ.out.tab" ];
	then
	echo processing" "$(basename $i )

	mkdir "Mapped/"$(basename $i ) 
	
	/mnt/g/LYL_RNA_Workspace/TEanalysis/SQuIRE/SQuIRE/squire/Map.py  -1 $i"/"$(basename $i )"_1.fastq" -2 $i"/"$(basename $i )"_2.fastq" -o "Mapped/"$(basename $i ) -f ./Ref -r 102 -p 36    
	fi
	done