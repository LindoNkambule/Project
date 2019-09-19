#!/usr/local/bash
REFERENCE="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
DICTIONARY=${REFERENCE%.*}
READ_1="data_read1.fq"
READ_2="data_read2.fq"
READS="aligned_reads"
PICARD="/usr/local/pipeline/Tools/picard.jar"
gatk="/usr/local/pipeline/Tools/gatk-4.1.0.0/gatk"
java="/usr/lib/jvm/java-11-openjdk-amd64/bin/java"

alias picard="$java -jar $PICARD"

##################################################
#    Prepare a FASTA file to use as reference    #
##################################################

#Generate the FASTA file index
samtools faidx ${REFERENCE}              # This will create a ".fai" file (FASTA index). This allows efficient random access to the reference bases.

#Generate a FASTA sequence dictionary file
picard CreateSequenceDictionary R= ${REFERENCE} O= ${DICTIONARY}.dict          #Create a dictionary of the contig names and sizes

#Generate the BWA index                This will generate all the files needed by BWA for alignment
bwa index -a bwtsw ${REFERENCE}          # bwtsw is an Algorithm implemented in BWT-SW. This method works with the whole human genome.


#########################
#     Read ALIGNMENT    #
########################

#Start read alignment
bwa mem ${REFERENCE} ${READ_1} ${READ_2} > ${READS}.sam    # This will produce a alignment.sam file

# Sorting SAM by coordinate and Converting SAM to BAM
picard SortSam \
    INPUT=aligned_reads.sam \
    OUTPUT=sorted_reads.bam \
    SORT_ORDER=coordinate   #The output contains the aligned reads (BAM) sorted by coordinate

# Mark MarkDuplicates
picard MarkDuplicates \
    INPUT=sorted_reads.bam \
    OUTPUT=dedup_reads.bam \
    METRICS_FILE=metrics.txt

#Add read groups. Assigns all the reads in a file to a single new read-group
picard AddOrReplaceReadGroups \
      I=dedup_reads.bam \
      O=output.bam \
      RGID=4 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=20

#Index the BAM file [creates a output.bam.bai (binary alignment index) file]
picard BuildBamIndex \
    INPUT=output.bam

########################
#    CALL VARIANTS     #
########################

# All the processes will be forked in parallel (&) and the script will wait ("wait") until all the processes have completed before moving to the next command
### Single sample calling ###
freebayes -f ${REFERENCE} output.bam >freeebayes.vcf &   # The & is for forking all the processes
bcftools mpileup -Ou -f ${REFERENCE} output.bam | bcftools call -mv -Ob | bcftools view > bcftools.vcf &
$gatk HaplotypeCaller \
    -R ${REFERENCE} \
    -I output.bam \
    -O gatk.vcf &
wait
echo "All processes complete"
