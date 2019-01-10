#!/usr/local/bash
REFERENCE="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
READ_1="data_read1.fq"
READ_2="data_read2.fq"
READS="alignment"
#########################
#     Read ALIGNMENT    #
########################

#### Generate the BWA index #####       This will generate all the files needed by BWA for alignment
bwa index -a bwtsw $REFERENCE          # bwtsw is an Algorithm implemented in BWT-SW. This method works with the whole human genome.

#### Generate the FASTA file index ####
samtools faidx $REFERENCE              # This will create a ".fai" file (FASTA index). This allows efficient random access to the reference bases.

#### Start read alignment ####
bwa mem $REFERENCE $READ_1 $READ_2 > $READS.sam    # This will produce a alignment.sam file

#### Converting SAM to BAM ####
samtools view -Sb $READS.sam > $READS.unsorted.bam

#### Sort the alignment file before calling variants ####
samtools sort $READS.unsorted.bam -o ANYTHING > $READS.bam

#### Index the alignment file [creates a $READS.bam.bai (binary alignment index) file] ####
samtools index $READS.bam



########################
#    CALL VARIANTS     #
########################

# All the processes will be forked in parallel (&) and the script will wait ("wait") until all the processes have completed before moving to the next command
Platypus.py callVariants --bamFiles = $READS.bam --refFile=$REFERENCE --output = platypus_variants.vcf  &       # The & is for forking all the processes
freebayes -f $ $REFERENCE $READS.bam >freeebayes_variants.vcf &
bcftools mpileup -f $REFERENCE $READS.bam | bcftools call -mv -Ob | bcftools view > bcftools_variants.vcf &

wait
echo "All processes complete"
