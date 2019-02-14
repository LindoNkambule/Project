#!/usr/local/bash
REFERENCE="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
READ_1="data_read1.fq"
READ_2="data_read2.fq"
READS="aligned_reads"
AF_THR="0.01"  # minimum allele frequency
BED="aligned"
PICARD="/usr/local/pipeline/Tools/picard.jar"
GATK="/usr/local/pipeline/Tools/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar"


##################################################
#    Prepare a FASTA file to use as reference    #
##################################################

#Generate the FASTA file index
samtools faidx $REFERENCE              # This will create a ".fai" file (FASTA index). This allows efficient random access to the reference bases.

#Generate a FASTA sequence dictionary file
java -jar $PICARD CreateSequenceDictionary R= $REFERENCE O= $REFERENCE.dict          #Create a dictionary of the contig names and sizes

#Generate the BWA index                This will generate all the files needed by BWA for alignment
bwa index -a bwtsw $REFERENCE          # bwtsw is an Algorithm implemented in BWT-SW. This method works with the whole human genome.


#########################
#     Read ALIGNMENT    #
########################

#Start read alignment
bwa mem $REFERENCE $READ_1 $READ_2 > $READS.sam    # This will produce a alignment.sam file

# Sorting SAM by coordinate and Converting SAM to BAM
java -jar $PICARD SortSam \
    INPUT=aligned_reads.sam \
    OUTPUT=sorted_reads.bam \
    SORT_ORDER=coordinate   #The output contains the aligned reads (BAM) sorted by coordinate

# Mark MarkDuplicates
java -jar $PICARD MarkDuplicates \
    INPUT=sorted_reads.bam \
    OUTPUT=dedup_reads.bam \
    METRICS_FILE=metrics.txt

#Index the alignment file [creates a $READS.bam.bai (binary alignment index) file]
java -jar picard.jar BuildBamIndex \
    INPUT=dedup_reads.bam

#Generating BED file to use for VarDict
bedtools bamtobed -i dedup_reads.bam > $BED.bed


########################
#    CALL VARIANTS     #
########################

# All the processes will be forked in parallel (&) and the script will wait ("wait") until all the processes have completed before moving to the next command
Platypus.py callVariants --bamFiles = $READS.bam --refFile=$REFERENCE --output = platypus_variants.vcf  &       # The & is for forking all the processes
freebayes -f $ $REFERENCE $READS.bam >freeebayes_variants.vcf &
bcftools mpileup -f $REFERENCE $READS.bam | bcftools call -mv -Ob | bcftools view > bcftools_variants.vcf &
vardict -G $REFERENCE -f $AF_THR -N sample_name -b $READS.bam -c 1 -S 2 -E 3 -g 4 $BED.bed | teststrandbias.R | var2vcf_valid.pl -N sample_name -E -f $AF_THR &
wait
echo "All processes complete"


#######################
#       STORAGE      #
######################
