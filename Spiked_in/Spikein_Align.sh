#!/bin/bash

module load bowtie2
module load bwa/0.7.4 

for file in $(<SpikeFiles.txt); do
	bowtie2 --maxins 2000 -fr -p 16 -q -x ./bowtieindex/ERCC92 -1 "../Long_insertSize_run/Spikein/$file/$file.rRNAFiltered.UnpairedRemoved_Read1.fastq" -2 "../Long_insertSize_run/Spikein/$file/$file.rRNAFiltered.UnpairedRemoved_Read2.fastq" -S "./Bowtie2_Sam/$file.Spikein.sam"
	#bwa mem -M -t 12 ERCC92 "../Long_insertSize_run/Spikein/$file/$file.rRNAFiltered.UnpairedRemoved_Read1.fastq" "../Long_insertSize_run/Spikein/$file/$file.rRNAFiltered.UnpairedRemoved_Read2.fastq" > "./Mem_Sam/$file.Spikein.sam"
	cat ./Bowtie2_Sam/$file.Spikein.sam | grep -v "^@" | cut -f 3 | grep -v "*" | sort | uniq -c > ./Bowtie2_Counts/$file.spikeins.bow.txt
	cat ./Bowtie2_Counts/$file.spikeins.bow.txt | wc -l >> SpikeinNumbersCounts.bow.txt
done
