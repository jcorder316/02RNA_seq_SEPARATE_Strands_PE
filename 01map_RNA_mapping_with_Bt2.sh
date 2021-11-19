#!/bin/bash
source  ~/.bash_profile

pa=path_to_raw_files

BA=path_toBAMS
tri=/home/epi4/001_Softwares/Trimmomatic-0.36
geno_mm10=/home/epi4/002_Genomes/01mm/mm10/mm10_bt2/bt2/mm10
ADAPTERS=/home/epi4/001_Softwares/Trimmomatic-0.36/adapters/TruSeq3-PE.fa
igv=path_to_IGV
#RE=/media/heart/KI-MA002/018_Heineke/01BGI_Runs/03HEI03_CM_EC/Reports/02EC
######### move report from Fastqc to the Report folders
FA=path_to_QC ## QC means quality control

start="15"
end="90"
min="15"
trim1="10"
trim2="70"
mini="15"
#rm trimmomatic.cmds
################################## this script will take the mapped bam files and separate the Strands!!! and produce a BigWig file, to see in IGV browser  ############

######## IMportant Delete all the files you will not need, Sam files for example  ##########
for R1 in $pa/*_2.fastq
 do
 R1_trim=${R1//.fastq.gz/_trim.fq.gz}
 yolo=$(echo "$R1" | rev | cut -c 7- | rev)
echo "${yolo}_POS.sam"
   #echo " java -jar $tri/trimmomatic-0.36.jar PE -threads 4 -phred33 ${yolo}_1.fastq.gz ${yolo}_2.fastq.gz ${yolo}_R1_paired.fq ${yolo}_R1_unpaired.fq ${yolo}_R2_paired.fq ${yolo}_R2_unpaired.fq ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$mini CROP:$trim2 HEADCROP:$trim1 >>"
java -jar $tri/trimmomatic-0.36.jar PE -threads 4 -phred33 ${yolo}_1.fastq ${yolo}_2.fastq ${yolo}_R1_paired.fq ${yolo}_R1_unpaired.fq.gz ${yolo}_R2_paired.fq ${yolo}_R2_unpaired.fq.gz ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$mini CROP:$trim2 HEADCROP:$trim1 >> trimmomatic.cmds
	yolo=$(echo "$R1" | rev | cut -c 7- | rev)
  bowtie2 -D 15 -R 3 -N 0 -L 20 -i S,1,0.56 -k 1 -p8 -x $geno_mm10 -1 ${yolo}_fw_pair.fq -2 ${yolo}_rw_pair.fq -S ${yolo}.sam
 
	samtools view -Sb -u ${yolo}.sam  |samtools sort - -o ${yolo}_s.bam
  				
		#samtools view -Sb ${yolo}.sam | samtools view -F 1804 -f 2 -u | samtools sort - > ${yolo}_F8f2_s.bam #this works when the samples are PE
			
			#java -jar $du/MarkDuplicates.jar INPUT=${yolo}_F8f2_s.bam OUTPUT=${yolo}_F8f2_s_rmd.bam METRICS_FILE=${yolo}_F8f2_s_rmd_dupl_INFO.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true
				
				samtools index ${yolo}_s.bam
			
		bamCoverage -b ${yolo}_s.bam -p max -bs 20 --smoothLength 40 --normalizeUsing RPKM -e 150 -o ${yolo}_RPKM_bt2_mm10.bw	
		
fastqc $R1
rm ${yolo}.sam
mv ${yolo}_R1.fastqc.html $FA
mv ${yolo}_R1.fastqc.zip $FA 
mv *.bw $igv
mv *.bam.bai $BA
mv ${yolo}_STAR_s.bam $BA
mv ${yolo}_STAR_ucsc $BA

	
done

