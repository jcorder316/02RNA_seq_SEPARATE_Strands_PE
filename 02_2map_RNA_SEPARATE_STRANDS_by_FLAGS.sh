#!/bin/bash
p=/media/heart/KI-MA001/NCBI_Data/004_RNA_seq_ENCO/07_Zf_PRJNA509429/danRer11/00sham
p4=/media/heart/KI-MA001/NCBI_Data/004_RNA_seq_ENCO/07_Zf_PRJNA509429/danRer11/04dpi

source ~/.bash_profile

tri=~/001_sofwares/Trimmomatic-0.36
mm9=~/002_Genomes/01_Mm/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome
geno=~/002_Genomes/01_Mm/mm9/ucsc_mm9/STAR/ 
geno10=~/002_Genomes/01_Mm/mm10/mm10_bt2_Zv9/
zv9=~/002_Genomes/02_ZF_danRer7/Danio_rerio_UCSC_danRer7/Danio_rerio/UCSC/danRer7/Sequence/Bowtie2Index/genome
zv9c=~/002_Genomes/02_ZF_danRer7/Zv9_Chr_sizes 
Eq=/media/heart/KI-MA001/dob26_b60c_sor/03_m


ADAPTERS=~/001_sofwares/Trimmomatic-0.36/adapters/TruSeq3-PE.fa	
p=path/to/bamfiles
trim1="10"
trim2="70"
mini="15"
#rm trimmomatic.cmds
################################## this script will take the mapped bam files and separate the Strands!!! and produce a BigWig file, to see in IGV browser  ############
for R1 in $p/*_danrer11_s.bam
 do
   
   R1_trim=${R1//.fastq.gz/_trim.fq.gz}
   yolo=$(echo "$R1" | rev | cut -c 7- | rev)
echo "${yolo}_POS.sam"
   #echo " java -jar $tri/trimmomatic-0.36.jar PE -threads 4 -phred33 ${yolo}_1.fastq.gz ${yolo}_2.fastq.gz ${yolo}_R1_paired.fq ${yolo}_R1_unpaired.fq ${yolo}_R2_paired.fq ${yolo}_R2_unpaired.fq ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$mini CROP:$trim2 HEADCROP:$trim1 >>"
#java -jar $tri/trimmomatic-0.36.jar PE -threads 4 -phred33 ${yolo}_1.fastq ${yolo}_2.fastq ${yolo}_R1_paired.fq ${yolo}_R1_unpaired.fq.gz ${yolo}_R2_paired.fq ${yolo}_R2_unpaired.fq.gz ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$mini CROP:$trim2 HEADCROP:$trim1 >> trimmomatic.cmds
		samtools view -h $R1 | awk '{if ($2==83||$2==163) print $0}'> ${yolo}_POS.sam
		samtools view -h $R1 | awk '{if ($2==99||$2==147) print $0}'> ${yolo}_NEG.sam
		samtools view -SH  $R1 > header
							cat header ${yolo}_POS.sam > ${yolo}_POS_2.sam
							cat header ${yolo}_NEG.sam > ${yolo}_NEG_2.sam

	samtools view -Sb -u ${yolo}_POS_2.sam | samtools sort - -o ${yolo}_POS_2_s.bam
	samtools view -Sb -u ${yolo}_NEG_2.sam | samtools sort - -o ${yolo}_NEG_2_s.bam

	
	bamToBed -i  ${yolo}_POS_2_s.bam >  ${yolo}_POS_2_s.bed
						
						makeTagDirectory ${yolo}_POS_ucsc ${yolo}_POS_2_s.bed

	bamToBed -i  ${yolo}_NEG_2_s.bam >  ${yolo}_NEG_2_s.bed
						
						makeTagDirectory ${yolo}_NEG_ucsc ${yolo}_NEG_2_s.bed

	samtools index ${yolo}_POS_2_s.bam
	samtools index ${yolo}_NEG_2_s.bam

	bamCoverage -b ${yolo}_POS_2_s.bam -bs 20 --smoothLength 40 -p max --normalizeUsing RPKM -e 150 -o ${yolo}_POS_RPKM_sm20_40.bw
bamCoverage -b ${yolo}_NEG_2_s.bam -bs 20 --smoothLength 40 -p max --normalizeUsing RPKM --scaleFactor -1 -e 150 -o ${yolo}_NEG_RPKM_sm20_40.bw
	
done

