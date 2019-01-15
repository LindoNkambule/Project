#!/usr/local/bash
REFERENCE="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
READ_1="data_read1.fq"
READ_2="data_read2.fq"
READS="alignment"
AF_THR="0.01"  # minimum allele frequency
BED="aligned"
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
samtools sort $READS.unsorted.bam -o $READS.bam

#### Index the alignment file [creates a $READS.bam.bai (binary alignment index) file] ####
samtools index $READS.bam

#### Generating BED file to use for VarDict ####
bedtools bamtobed -i $READS.bam > $BED.bed

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
mkdir Reference/; mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna* alignment.* Reference/       # Make a directory called Reference and move all the reference and alignment files to the folder
mkdir Variants/; mv *.vcf Variants/  # Make a directory called Variants and move all the vcf files into it

zip -r Variants.zip Variants/  # Zip the folder for long term storage
zip -r Reference.zip Reference/ # Zip the folder for long term storage 
