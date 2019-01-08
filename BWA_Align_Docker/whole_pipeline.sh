#!/usr/local/bash
REFERENCE="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
READ_1="read1.fq"
READ_2="read2.fq"
READS="alignment"
#########################
#     Read ALIGNMENT    #
########################

echo "Indexing the reference genome"
bwa index -a bwtsw $REFERENCE
echo "Done indexing"

echo "Start read alignment"
bwa mem $REFERENCE $READ_1 $READ_2 > $READS.sam
echo "Done aligning"

# Converting SAM to BAM
echo "Converting SAM ===> BAM"
samtools view -Sb $READS.sam > $READS.unsorted.bam
echo "Done converting SAM ===> BAM"

# Sort the alignment file before calling variants
echo "Sorting BAM file"
samtools sort $READS.unsorted.bam -o ANYTHING > $READS.bam
echo "Done sorting BAM file"

# Index the alignment file [creates a $READS.bam.bai (binary alignment index) file]
echo "Start indexing BAM file"
samtools index $READS.bam
echo "Done indexing BAM file"


########################
#    CALL VARIANTS     #
########################

# All the processes will be forked in parallel (&) and the script will wait ("wait") until all the processes have completed before moving to the next command
Platypus.py callVariants --bamFiles = $READS.bam --refFile=$REFERENCE --output = platypus_variants.vcf  &       # The & is for forking all the processes
freebayes -f $ $REFERENCE $READS.bam >freeebayes_variants.vcf &
bcftools mpileup -f $REFERENCE $READS.bam | bcftools call -mv -Ob | bcftools view > bcftools_variants.vcf &

wait
echo "All processes complete"
