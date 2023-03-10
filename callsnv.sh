#!/bin/sh
#Usage cleaner.sh file_prefix
#1. clean fastq file
cutadapt -a GGGGGGGGGG -A GGGGGGGGGG -o $1.R1.new.fastq.gz -p $1.R2.new.fastq.gz $1.R1.fastq.gz $1.R2.fastq.gz
cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o $1.R1.new2.fastq.gz -p $1.R2.new2.fastq.gz $1.R1.new.fastq.gz $1.R2.new.fastq.gz
#2. mapping reads
bwa mem -t 64 -M -Y -R "@RG\tID:${1}\tPL:ILLUMINA\tPU:NEXTSEQ\tLB:${1}library\tSM:${1}" /data/reference/hg38/gatk/Homo_sapiens_assembly38.fasta $1.R1.fastq.gz $1.R2.fastq.gz | samtools sort -m 4G -@ 16 -o ./$1.bam
samtools index -@ 16 ./$1.bam ./$1.bam.bai
#3. mark duplicates
java -jar /bin/picard.jar MarkDuplicates I=./$1.bam O=./$1.mkdup.bam M=./$1.mtx.txt
#3.5. call CNV --with reference
cnvkit.py batch ./T1000.mkdup.bam -n ./N1000.mkdup.bam \ 
	-m wgs -p 64 -f /data/reference/hg38/gatk/Homo_sapiens_assembly38.fasta \
	--annotate /ref/hg38_refFlat.txt 
cnvkit.py batch *.bam -r /10x_colon_RNA/CRC1/scDNA/sam/reference.cnn -p 8 -d ./scsam
#3.5. call CNV --no reference
cnvkit.py batch *Tumor.bam -n -f hg38.fasta \
	-m wgs -p 64 -f /data/reference/hg38/gatk/Homo_sapiens_assembly38.fasta \
    --annotate refFlat.txt \
    --output-reference ./my_flat_reference.cnn -d example3/
cnvkit.py batch *.bam -r /data/scDNA/sam/reference.cnn -p 8 -d ./scsam
#4. BQSR
gatk BaseRecalibrator -R /data/reference/hg38/gatk/Homo_sapiens_assembly38.fasta --known-sites /ref/dbSNP146_hg38.vcf.gz --known-sites /ref/gnomAD30_hg38.vcf.gz -I ./$1.mkdup.bam -O ./$1.mkdup.metrics
gatk ApplyBQSR -bqsr ./$1.mkdup.metrics -R /data/reference/hg38/gatk/Homo_sapiens_assembly38.fasta -I ./$1.mkdup.bam -O $1.mkdup.BQSR.bam 
#5.
gatk Mutect2 -R reference.fa -I sample.bam --native-pair-hmm-threads 32 --tmp-dir /home/data/tmp -O single_sample.vcf.gz

gatk Mutect2 -R /data/reference/hg38/gatk/Homo_sapiens_assembly38.fasta -I N1000.mkdup.BQSR.bam --max-mnp-distance 0 -O normal.vcf.gz

gatk Mutect2   -R /data/reference/hg38/gatk/Homo_sapiens_assembly38.fasta   \
	-I /10x_colon_RNA/CRC1/scDNA/sam/T1000.mkdup.BQSR.bam  \
	--germline-resource  /ref/gnomAD30_hg38.vcf.gz  \
	--panel-of-normals /10x_colon_RNA/CRC1/scDNA/sam/normal.vcf.gz \
	--native-pair-hmm-threads 64  --tmp-dir /home/data/tmp \
	-O /10x_colon_RNA/CRC1/scDNA/sam/T1000_panel.vcf.gz
#6.
gatk FilterMutectCalls -R ref.fasta -V unfiltered.vcf -O filtered.vcf
